## run Simulation pipeline

WD='</path/to/AmpliPiper>'
# e.g. WD=/media/inter/pipelines/AmpliPiper

## install dependencies
mamba create \
    -p ${WD}/simulations/scripts/nanosim_env
conda activate ${WD}/simulations/scripts/nanosim_env

mamba install \
    -y \
    -c conda-forge \
    -c bioconda \
    scikit-learn regex six samtools pysam pybedtools minimap2 joblib htseq genometools-genometools piecewise-regression biopython edlib vsearch pandas

## install NanoSim
cd ${WD}/simulations/scripts
git clone https://github.com/bcgsc/NanoSim.git
cd NanoSim

## run pipleine

mkdir ${WD}/simulations/data

# get Sanger COI fa from here /media/inter/mkapun/projects/HAPLOTYPES-Benchmark/Syrphid/data/Syrphidae_COI_Sanger_red.fasta, append the fwd and rev primers
echo """>Stor-1
GGTCAACAAATCATAAAGATATTGGATATAAAATTGGATCTCCTCCTCCTGCTGGGTCAAAGAATGAAGTATTTAAATTTCGATCTGTTAATAATATAGTAATAGCTCCTGCTAAAACTGGTAATGATAGAAGTAATAATAAAGCTGTAATTACTACTGATCATACAAATAAAGGTATTCGATCATAAGTAAGTCCACTAGAACGTATATTAATTACTGTTGTAATAAAATTTACTGCTCCTAAAATTGAAGATATTCCTGCCAAATGAAGAGAAAAAATAGCTAAATCAACAGAAGCTCCTCCATGTGCAATACTAGCAGAAAGAGGAGGGTAAACTGTTCAACCTGTTCCAGCTCCGTTTTCAACTATTCTTCTTACTAATAAAAGTGTTAAAGAAGGAGGTAATAATCAAAAACTTATATTATTTATTCGAGGAAATGCTATATCAGGAGCTCCTAACATTAAAGGAACTAATCAATTTCCAAATCCTCCAATTATAATTGGTATTACTATAAAAAAAATTATTACAAAAGCATGTGCAGTTACAATAACATTATAAATTTGGTCATCCCCAATTAAAGCACCTGGATGTCCAAGTTCTGCACGAATTAATACACTTAAAGAGGTTCCTACTATACCTGCTCAAGTTCCAAATAAAAAATATAATGTTCAACACCTATTCTGATTTTTTGG""" >${WD}/simulations/data/Sy_tor_31_COI_short.fa

### Here, I am randomly adding mutations according to the divergence proportions defined below in the region WITHIN the primers (but) not in the primers. In addition, I append 50 A's at the 5' and 3' end of the sequence. I call this sequence "long" (corresponding to the reference) and the sequence without A's "short" (which is corresponding to a transcript). See usage below for the training and simulations
python ${WD}/simulations/scripts/mutateRef.py \
    --primers GGTCAACAAATCATAAAGATATTGG,CCAAAAAATCARAATARRTGYTG \
    --input ${WD}/simulations/data/Sy_tor_31_COI_short.fa \
    --divergence 0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5 \
    --HowManyAs 50 \
    --output ${WD}/simulations/data/Sy_tor_31_COI

conda activate /media/inter/pipelines/AmpliPiper/envs/python_dependencies

## demultiplex to get COX-specific reads
mkdir -p ${WD}/simulations/data/demulti
python ${WD}/scripts/DemultFastq.py \
    -i ${WD}/testdata/reads/Sy_tor_31.fastq.gz \
    -th 0.15 \
    -sr 150 \
    -rp 1000 \
    -p ${WD}/testdata/data/primers.csv \
    -o ${WD}/simulations/data/demulti \
    -mr 100

conda deactivate

conda activate ${WD}/simulations/scripts/nanosim_env

mkdir ${WD}/simulations/training_trans

## train error model based on real data
python3 ${WD}/simulations/scripts/NanoSim/src/read_analysis.py transcriptome \
    -rg ${WD}/simulations/data/Sy_tor_31_COI_long.fa \
    -rt ${WD}/simulations/data/Sy_tor_31_COI_short.fa \
    -i ${WD}/simulations/data/Sy_tor_31.fq.gz \
    -o ${WD}/simulations/training_trans/training_COI \
    --no_intron_retention \
    --fastq \
    -t 200

printf """target_id\test_counts\tpm\nStor-1\t1000\t1000\n""" >${WD}/simulations/training_trans/Sy_tor_31_cox1.exp

### now simulate reads for each of the unmutated and all mutated reference sequences
for DIV in 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5; do

    mkdir ${WD}/simulations/simulation_trans_${DIV}

    conda activate ${WD}/simulations/scripts/nanosim_env

    python3 ${WD}/simulations/scripts/NanoSim/src/simulator.py transcriptome \
        -rg ${WD}/simulations/data/Sy_tor_31_COI_long_${DIV}.fa \
        -rt ${WD}/simulations/data/Sy_tor_31_COI_short_${DIV}.fa \
        -c ${WD}/simulations/training_trans/training_COI \
        -o ${WD}/simulations/simulation_trans_${DIV}/COI_simulated_${DIV} \
        -e ${WD}/simulations/training_trans/Sy_tor_31_cox1.exp \
        -n 500000 \
        --no_model_ir \
        -t 200 \
        -min 690 \
        -max 692 \
        --fastq

    ## only keep aligned sequences
    gzip ${WD}/simulations/simulation_trans_${DIV}/COI_simulated_${DIV}_aligned_error_profile
    gzip ${WD}/simulations/simulation_trans_${DIV}/COI_simulated_${DIV}_aligned_reads.fastq
    rm -rf ${WD}/simulations/simulation_trans_${DIV}/COI_simulated_${DIV}_unaligned_reads.fastq

done

## make primers file
echo """ID,FWD,REV,SIZE
nCOX1,GGTCAACAAATCATAAAGATATTGG,CCAAAAAATCARAATARRTGYTG,690
""" >${WD}/simulations/data/Primers.csv

# test the sensitivity with default parameters and diploidy

## do 20 replicates
for i in {1..20}; do

    ## Pick 40k random aligned reads WITHOUT mutations (seed set to ${i} for reproducablity)
    (
        gunzip -c ${WD}/simulations/simulation_trans_0/COI_simulated_0_aligned_reads.fastq.gz |
            paste - - - - |
            shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:${i} -nosalt </dev/zero 2>/dev/null) -n 40000 |
            tr '\t' '\n'
    ) >${WD}/simulations/simulation_trans_0/COI_simulated_0_aligned_reads_40k_${i}.fastq

    ## loop over all divergences
    for DIV in 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5; do

        mkdir -p ${WD}/simulations/results/Divergence/amplipiper_test_0-${DIV}_${i}/data

        ## combined 40k unmutated reads with 40k randomly drawn reads from mutated dataset to get 50:50 ratio
        (
            cat ${WD}/simulations/simulation_trans_0/COI_simulated_0_aligned_reads_40k_${i}.fastq
            ## randomly subsample 40k reads from the input file and use $i as seed
            gunzip -c ${WD}/simulations/simulation_trans_${DIV}/COI_simulated_${DIV}_aligned_reads.fastq.gz |
                paste - - - - |
                shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:${i} -nosalt </dev/zero 2>/dev/null) -n 40000 |
                tr '\t' '\n'
        ) | gzip >${WD}/simulations/results/Divergence/amplipiper_test_0-${DIV}_${i}/data/amplipiper_test_0-${DIV}_${i}.fq.gz

        ## create samples file
        echo "amplipiper_test_0-${DIV},"${WD}/simulations/results/Divergence/amplipiper_test_0-${DIV}_${i}/data/amplipiper_test_0-${DIV}_${i}.fq.gz \
            >${WD}/simulations/results/Divergence/amplipiper_test_0-${DIV}_${i}/data/amplipiper_test_0-${DIV}_${i}.csv

        conda deactivate

        conda activate /media/inter/pipelines/AmpliPiper/envs/python_dependencies

        ## run AmpliPiper
        bash /media/inter/pipelines/AmpliPiper/shell/AmpliPiper.sh \
            --samples ${WD}/simulations/results/Divergence/amplipiper_test_0-${DIV}_${i}/data/amplipiper_test_0-${DIV}_${i}.csv \
            --primers ${WD}/simulations/data/Primers.csv \
            --output ${WD}/simulations/results/Divergence/amplipiper_test_0-${DIV}_${i}/analysis \
            --quality 10 \
            --nreads 10000 \
            --threads 1 \
            --kthreshold 0.05 \
            --minreads 50 \
            --sizerange 100 \
            --force &

        conda deactivate

    done
done

## summarize the output and count from Amplicon_Sorter the number of haplotypes and groups (i.e. differernt detected loci)
mkdir ${WD}/simulations/results/summary
printf "Divergence Replicate Groups Haplotypes\n" >${WD}/simulations/results/summary/HapCount.txt
for i in {1..20}; do
    for DIV in 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5; do
        C=0
        for j in ${WD}/simulations/results/Divergence/amplipiper_test_0-${DIV}_${i}/analysis/results/consensus_seqs/amplipiper_test_0-${DIV}/COX1/COX1_?_consensussequences.fasta; do
            HAP=$(grep '^>' ${j} | wc -l)
            let C++
            echo $DIV $i $C $HAP >>${WD}/simulations/results/summary/HapCount.txt
        done
    done

done

## plot with R
Rscript -e """

library(tidyverse)
library(ggpubr)
DATA<-read.table('${WD}/simulations/results/summary/HapCount.txt',header=T)

long <- DATA %>% 
  pivot_longer(
    cols = Groups:Haplotypes, 
    names_to = 'Type',
    values_to = 'Value'
)

DATA.mean <- long %>%
    group_by(Divergence,Type) %>%
    summarise(Mean=mean(Value),SD=sd(Value))

PLOT<-ggplot(DATA.mean,aes(x=Divergence,y=Mean))+
    geom_bar( stat='identity', position='dodge',alpha=0.7,fill='#21614A') +
    geom_errorbar( stat='identity', position='dodge', aes(
            ymin=Mean-SD, 
            ymax=Mean+SD),  
        colour='black', 
        alpha=0.9)+
    ylab('Count')+
    facet_grid(Type~.,scales = 'free_y')+
    theme_bw() 
PLOT

ggsave('${WD}/simulations/results/summary/HapCount.png',
    PLOT,
    width=10,
    height=5)

ggsave('${WD}/simulations/results/summary/HapCount.pdf',
    PLOT,
    width=10,
    height=5)
"""

## Now test for the effect of the --similar_consensus parameter in Amplicon_Sorter (newly added tp AmpliPiper)

## make 20 replicates
for i in {1..20}; do

    ## loop through all divergences
    for DIV in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1; do

        ## test range of --similar_consensus parameters
        for SC in 96 96.5 97 97.5 97.5 97.5.5 99; do

            conda activate /media/inter/pipelines/AmpliPiper/envs/python_dependencies

            ## run AmpliPiper
            bash /media/inter/pipelines/AmpliPiper/shell/AmpliPiper.sh \
                --samples ${WD}/simulations/results/Divergence/amplipiper_test_0-${DIV}_${i}/data/amplipiper_test_0-${DIV}_${i}.csv \
                --primers ${WD}/simulations/data/Primers.csv \
                --output ${WD}/simulations/results/SimilarConsensus/amplipiper_test_0-${DIV}_${i}_SC${SC}/analysis \
                --quality 10 \
                --nreads 10000 \
                --similar_consensus ${SC} \
                --threads 1 \
                --kthreshold 0.05 \
                --minreads 50 \
                --sizerange 100 \
                --force &

            conda deactivate

        done
    done
done

## summarize
printf "Divergence SimilarConsensus Replicate Groups Haplotypes\n" >${WD}/simulations/results/summary/HapCount_SC.txt
for i in {1..20}; do
    for DIV in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1; do

        for SC in 96 96.5 97 97.5 97.5 97.5.5 99; do
            C=0
            for j in ${WD}/simulations/results/SimilarConsensus/amplipiper_test_0-${DIV}_${i}_SC${SC}/analysis/results/consensus_seqs/amplipiper_test_0-${DIV}/COX1/COX1_?_consensussequences.fasta; do
                HAP=$(grep '^>' ${j} | wc -l)
                let C++
                echo $DIV $SC $i $C $HAP >>${WD}/simulations/results/summary/HapCount_SC.txt
            done
        done
    done
done

## plot
Rscript -e """
library(tidyr)
library(tidyverse)
DATA<-read.table('${WD}/simulations/results/summary/HapCount_SC.txt',header=T)

DATA.mean <- DATA %>%
    group_by(Divergence,SimilarConsensus) %>%
    summarise(Mean=mean(Haplotypes),SD=sd(Haplotypes))

# Define your three colors
three_colors <- c('#21614A', '#FF8E3E','#7AC1DC')

# Create a discrete color palette with n levels (same as the number of categories)
n_levels <- length(levels(as.factor(DATA.mean$SimilarConsensus)))
my_palette <- colorRampPalette(three_colors)(n_levels)

PLOT<-ggplot(DATA.mean,aes(x=as.factor(Divergence),y=Mean,fill=as.factor(SimilarConsensus)))+
    geom_bar( stat='identity', position='dodge',alpha=0.7) +
    geom_errorbar( stat='identity', position='dodge', aes(
            ymin=Mean-SD, 
            ymax=Mean+SD),  
        colour='black', 
        alpha=0.9)+
    xlab('Divergence')+
    ylab('Number of Haplotypes')+
    scale_fill_manual(name='SC',
    values=my_palette)+
    theme_bw() 
PLOT

ggsave('${WD}/simulations/results/summary/HapCount_SC.png',
    PLOT,
    width=10,
    height=2.5)

ggsave('${WD}/simulations/results/summary/HapCount_SC.pdf',
    PLOT,
    width=10,
    height=2.5)
"""

printf "Divergence,Replicate,SimilarConsensus,Frequency,NumberHaplotypes,NumberMM,PropMM\n" >${WD}/simulations/results/summary/SeqSimilarity.csv

## make 10 replicates
for i in {1..20}; do

    ## loop through all divergences
    for DIV in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1; do

        ## test range of --similar_consensus parameters
        for SC in 96 96.5 97 97.5 97.5 97.5.5 99; do

            python ${WD}/simulations/scripts/compare2Sanger4sim.py \
                --infile ${WD}/simulations/results/SimilarConsensus/amplipiper_test_0-${DIV}_${i}_SC${SC}/analysis/results/consensus_seqs/amplipiper_test_0-${DIV}/COX1/COX1_consensussequences.fasta \
                --sanger ${WD}/simulations/data/Sy_tor_31_COI_short_${DIV}.fa \
                --primer GGTCAACAAATCATAAAGATATTGG,CCAAAAAATCARAATARRTGYTG \
                --frequency 0.5 \
                --divergence ${DIV} \
                --SC ${SC} \
                --replicate ${i} >>${WD}/simulations/results/summary/SeqSimilarity.csv

        done
    done
done

## plot
Rscript -e """
library(tidyr)
library(tidyverse)
DATA<-read.csv('${WD}/simulations/results/summary/SeqSimilarity.csv',header=T)

DATA.mean <- DATA %>%
    group_by(Divergence,SimilarConsensus) %>%
    summarise(Mean=mean(PropMM),SD=sd(PropMM))

# Define your three colors
three_colors <- c('#21614A', '#FF8E3E','#7AC1DC')

# Create a discrete color palette with n levels (same as the number of categories)
n_levels <- length(levels(as.factor(DATA.mean$SimilarConsensus)))
my_palette <- colorRampPalette(three_colors)(n_levels)

PLOT<-ggplot(DATA.mean,aes(x=as.factor(Divergence),y=Mean,fill=as.factor(SimilarConsensus)))+
    geom_bar( stat='identity', position='dodge',alpha=0.7) +
    geom_errorbar( stat='identity', position='dodge', aes(
            ymin=Mean-SD, 
            ymax=Mean+SD),  
        colour='black', 
        alpha=0.9)+
    xlab('Divergence')+
    ylab('Proportion Mismatches')+
    scale_fill_manual(name='SC',
    values=my_palette)+
    theme_bw() 
PLOT

ggsave('${WD}/simulations/results/summary/SeqSimilarity.png',
    PLOT,
    width=10,
    height=2.5)

ggsave('${WD}/simulations/results/summary/SeqSimilarity.pdf',
    PLOT,
    width=10,
    height=2.5) 
"""

### now mutate frequencies

## do 10 replicates
for i in {11..20}; do

    for Freq in 1 2 3 4 5; do

        FIRST=$((40000 * ${Freq} / 10))
        SECOND=$((40000 - FIRST))

        ## Pick 40k random aligned reads WITHOUT mutations (seed set to ${i} for reproducablity)
        (
            gunzip -c ${WD}/simulations/simulation_trans_0/COI_simulated_0_aligned_reads.fastq.gz |
                paste - - - - |
                shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:${i} -nosalt </dev/zero 2>/dev/null) -n ${FIRST} |
                tr '\t' '\n'
        ) >${WD}/simulations/simulation_trans_0/COI_simulated_0_aligned_reads_40k_${i}_${Freq}.fastq

        ## loop over all divergences
        for DIV in 0 0.01 0.02 0.03 0.04 0.05; do

            mkdir -p ${WD}/simulations/results/Frequency/amplipiper_test_0-${DIV}_${i}/data

            ## combined
            (
                cat ${WD}/simulations/simulation_trans_0/COI_simulated_0_aligned_reads_40k_${i}_${Freq}.fastq
                ## randomly subsample 40k reads from the input file and use $i as seed
                gunzip -c ${WD}/simulations/simulation_trans_${DIV}/COI_simulated_${DIV}_aligned_reads.fastq.gz |
                    paste - - - - |
                    shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:${i} -nosalt </dev/zero 2>/dev/null) -n ${SECOND} |
                    tr '\t' '\n'
            ) | gzip >${WD}/simulations/results/Frequency/amplipiper_test_0-${DIV}_${i}/data/amplipiper_test_0-${DIV}_${i}_${Freq}.fq.gz

            ## create samples file
            echo "amplipiper_test_0-${DIV},"${WD}/simulations/results/Frequency/amplipiper_test_0-${DIV}_${i}/data/amplipiper_test_0-${DIV}_${i}_${Freq}.fq.gz \
                >${WD}/simulations/results/Frequency/amplipiper_test_0-${DIV}_${i}/data/amplipiper_test_0-${DIV}_${i}_${Freq}.csv

            conda deactivate

            conda activate /media/inter/pipelines/AmpliPiper/envs/python_dependencies

            ## run AmpliPiper
            bash /media/inter/pipelines/AmpliPiper/shell/AmpliPiper.sh \
                --samples ${WD}/simulations/results/Frequency/amplipiper_test_0-${DIV}_${i}/data/amplipiper_test_0-${DIV}_${i}_${Freq}.csv \
                --primers ${WD}/simulations/data/Primers.csv \
                --output ${WD}/simulations/results/Frequency/amplipiper_test_0-${DIV}_${i}_${Freq}/analysis_97.5 \
                --quality 10 \
                --nreads 10000 \
                --threads 1 \
                --kthreshold 0.05 \
                --similar_consensus 97.5 \
                --minreads 50 \
                --sizerange 100 \
                --force &

            conda deactivate

        done
    done
done

## summarize
printf "Divergence Frequency Replicate Groups Haplotypes\n" >${WD}/simulations/results/summary/HapCount_Frequency_97.5.txt
for i in {1..20}; do
    for Freq in 1 2 3 4 5; do
        for DIV in 0 0.01 0.02 0.03 0.04 0.05; do
            C=0
            for j in ${WD}/simulations/results/Frequency/amplipiper_test_0-${DIV}_${i}_${Freq}/analysis_97.5/results/consensus_seqs/amplipiper_test_0-${DIV}/COX1/COX1_?_consensussequences.fasta; do
                HAP=$(grep '^>' ${j} | wc -l)
                let C++
                echo $DIV $Freq $i $C $HAP >>${WD}/simulations/results/summary/HapCount_Frequency_97.5.txt
            done
        done
    done
done

## plot
echo """
library(tidyr)
library(tidyverse)
DATA<-read.table('${WD}/simulations/results/summary/HapCount_Frequency_97.5.txt',header=T)
DATA$Frequency<-DATA$Frequency/10

DATA.mean <- DATA %>%
    group_by(Divergence,Frequency) %>%
    summarise(Mean=mean(Haplotypes),SD=sd(Haplotypes))

# Define your three colors
three_colors <- c('grey90', '#FF8E3E')

# Create a discrete color palette with n levels (same as the number of categories)
n_levels <- length(levels(as.factor(DATA.mean$Frequency)))
my_palette <- colorRampPalette(three_colors)(n_levels)

PLOT<-ggplot(DATA.mean,aes(x=as.factor(Divergence),y=Mean,fill=as.factor(Frequency)))+
    geom_bar( stat="identity", position="dodge",alpha=0.7) +
    geom_errorbar( stat="identity", position="dodge", aes(
            ymin=Mean-SD, 
            ymax=Mean+SD),  
        colour="black", 
        alpha=0.9)+
    xlab("Divergence")+
    ylab("Number of Haplotypes")+
    scale_fill_manual(name="Frequency",
    values=my_palette)+
    theme_bw() 
PLOT
ggsave('${WD}/simulations/results/summary/HapCount_Frequency_97.5.png',
    PLOT,
    width=10,
    height=2.5)
ggsave('${WD}/simulations/results/summary/HapCount_Frequency_97.5.pdf',
    PLOT,
    width=10,
    height=2.5)
"""

printf "Divergence,Replicate,SimilarConsensus,Frequency,NumberHaplotypes,NumberMM,PropMM\n" >${WD}/simulations/results/summary/SeqSimilarity_Freq.csv

## make 10 replicates
for i in {1..20}; do
    for Freq in 1 2 3 4 5; do
        for DIV in 0 0.01 0.02 0.03 0.04 0.05; do
            python ${WD}/simulations/scripts/compare2Sanger4sim.py \
                --infile ${WD}/simulations/results/Frequency/amplipiper_test_0-${DIV}_${i}_${Freq}/analysis_97.5/results/consensus_seqs/amplipiper_test_0-${DIV}/COX1/COX1_consensussequences.fasta \
                --sanger ${WD}/simulations/data/Sy_tor_31_COI_short_${DIV}.fa \
                --primer GGTCAACAAATCATAAAGATATTGG,CCAAAAAATCARAATARRTGYTG \
                --frequency ${Freq} \
                --divergence ${DIV} \
                --SC 97.5 \
                --replicate ${i} >>${WD}/simulations/results/summary/SeqSimilarity_Freq.csv

        done
    done
done

## plot
Rscript -e """library(tidyr)
library(tidyverse)
DATA<-read.csv('${WD}/simulations/results/summary/SeqSimilarity_Freq.csv',header=T)
DATA$Frequency<-DATA$Frequency/10

DATA.mean <- DATA %>%
    group_by(Divergence,Frequency) %>%
    summarise(Mean=mean(PropMM),SD=sd(PropMM))

# Define your three colors
three_colors <- c('grey90', '#FF8E3E')

# Create a discrete color palette with n levels (same as the number of categories)
n_levels <- length(levels(as.factor(DATA.mean$Frequency)))
my_palette <- colorRampPalette(three_colors)(n_levels)

PLOT<-ggplot(DATA.mean,aes(x=as.factor(Divergence),y=Mean,fill=as.factor(Frequency)))+
    geom_bar( stat='identity', position='dodge',alpha=0.7) +
    geom_errorbar( stat='identity', position='dodge', aes(
            ymin=Mean-SD, 
            ymax=Mean+SD),  
        colour='black', 
        alpha=0.9)+
    xlab('Divergence')+
    ylab('Proportion Mismatches')+
    scale_fill_manual(name='Frequency',
    values=my_palette)+
    theme_bw() 
PLOT

ggsave('${WD}/simulations/results/summary/SeqSimilarity_Frequency.png',
    PLOT,
    width=10,
    height=2.5)

ggsave('${WD}/simulations/results/summary/SeqSimilarity_Frequency.pdf',
    PLOT,
    width=10,
    height=2.5)
"""

### finally, test influence of read dept

## do 10 replicates
for i in {1..20}; do

    for Freq in 10 20 50 100 500 1000 5000; do

        ## loop over all divergences
        for DIV in 0 0.01 0.02 0.03 0.04 0.05; do

            conda activate /media/inter/pipelines/AmpliPiper/envs/python_dependencies

            ## run AmpliPiper
            bash /media/inter/pipelines/AmpliPiper/shell/AmpliPiper.sh \
                --samples ${WD}/simulations/results/Divergence/amplipiper_test_0-${DIV}_${i}/data/amplipiper_test_0-${DIV}_${i}.csv \
                --primers ${WD}/simulations/data/Primers.csv \
                --output ${WD}/simulations/results/ReadDepth/amplipiper_test_0-${DIV}_${i}_${Freq}/analysis_97.5 \
                --quality 10 \
                --nreads ${Freq} \
                --threads 1 \
                --kthreshold 0.05 \
                --similar_consensus 97.5 \
                --minreads 50 \
                --sizerange 100 \
                --force &

            conda deactivate

        done
    done
done

## summarize
printf "Divergence ReadDepth Replicate Groups Haplotypes\n" >${WD}/simulations/results/summary/HapCount_ReadDepth.txt
for i in {1..20}; do
    for Freq in 10 20 50; do
        for DIV in 0 0.01 0.02 0.03 0.04 0.05; do
            C=0
            for j in ${WD}/simulations/results/ReadDepth/amplipiper_test_0-${DIV}_${i}_${Freq}/analysis_97.5/results/consensus_seqs/amplipiper_test_0-${DIV}/NCOX1/NCOX1_consensussequences.fasta; do
                HAP=$(grep '^>' ${j} | wc -l)
                let C++
                echo $DIV $Freq $i $C $HAP >>${WD}/simulations/results/summary/HapCount_ReadDepth.txt
            done
        done
    done
done

for i in {1..20}; do
    for Freq in 100 500 1000 5000; do
        for DIV in 0 0.01 0.02 0.03 0.04 0.05; do
            C=0
            for j in ${WD}/simulations/results/ReadDepth/amplipiper_test_0-${DIV}_${i}_${Freq}/analysis_97.5/results/consensus_seqs/amplipiper_test_0-${DIV}/*COX1/*COX1_?_consensussequences.fasta; do
                HAP=$(grep '^>' ${j} | wc -l)
                let C++
                echo $DIV $Freq $i $C $HAP >>${WD}/simulations/results/summary/HapCount_ReadDepth.txt
            done
        done
    done
done

## plot
Rscript -e """
library(tidyr)
library(tidyverse)
DATA<-read.table('${WD}/simulations/results/summary/HapCount_ReadDepth.txt',header=T)

DATA.mean <- DATA %>%
    group_by(Divergence,ReadDepth) %>%
    summarise(Mean=mean(Haplotypes),SD=sd(Haplotypes))

# Define your three colors
three_colors <- c('#FF8E3E','black')

# Create a discrete color palette with n levels (same as the number of categories)
n_levels <- length(levels(as.factor(DATA.mean$ReadDepth)))
my_palette <- colorRampPalette(three_colors)(n_levels)

PLOT<-ggplot(DATA.mean,aes(x=as.factor(Divergence),y=Mean,fill=as.factor(ReadDepth)))+
    geom_bar( stat='identity', position='dodge',alpha=0.7) +
    geom_errorbar( stat='identity', position='dodge', aes(
            ymin=Mean-SD, 
            ymax=Mean+SD),  
        colour='black', 
        alpha=0.9)+
    xlab('Divergence')+
    ylab('Number of Haplotypes')+
    scale_fill_manual(name='Read Depth',
    values=my_palette)+
    theme_bw() 
PLOT

ggsave('${WD}/simulations/results/summary/HapCount_ReadDepth.png',
    PLOT,
    width=10,
    height=2.5)

ggsave('${WD}/simulations/results/summary/HapCount_ReadDepth.pdf',
    PLOT,
    width=10,
    height=2.5)
"""
conda activate /media/inter/pipelines/AmpliPiper/envs/python_dependencies

printf "Divergence,Replicate,SimilarConsensus,ReadDepth,NumberHaplotypes,NumberMM,PropMM\n" >${WD}/simulations/results/summary/SeqSimilarity_ReadDepth.csv

## make 10 replicates
for i in {1..20}; do
    for Freq in 10 20 50; do
        for DIV in 0 0.01 0.02 0.03 0.04 0.05; do
            if [[ ! -f "${WD}/simulations/results/ReadDepth/amplipiper_test_0-${DIV}_${i}_${Freq}/analysis_97.5/results/consensus_seqs/amplipiper_test_0-${DIV}/NCOX1/NCOX1_consensussequences.fasta" ]]; then
                printf "${DIV},${i},97.5,${Freq},0,0,0\n" >>${WD}/simulations/results/summary/SeqSimilarity_ReadDepth.csv
            else
                python ${WD}/simulations/scripts/compare2Sanger4sim.py \
                    --infile ${WD}/simulations/results/ReadDepth/amplipiper_test_0-${DIV}_${i}_${Freq}/analysis_97.5/results/consensus_seqs/amplipiper_test_0-${DIV}/NCOX1/NCOX1_consensussequences.fasta \
                    --sanger ${WD}/simulations/data/Sy_tor_31_COI_short_${DIV}.fa \
                    --primer GGTCAACAAATCATAAAGATATTGG,CCAAAAAATCARAATARRTGYTG \
                    --frequency ${Freq} \
                    --divergence ${DIV} \
                    --SC 97.5 \
                    --replicate ${i} >>${WD}/simulations/results/summary/SeqSimilarity_ReadDepth.csv
            fi
        done
    done
done

for i in {1..20}; do
    for Freq in 100 500 1000 5000; do
        for DIV in 0 0.01 0.02 0.03 0.04 0.05; do
            python ${WD}/simulations/scripts/compare2Sanger4sim.py \
                --infile ${WD}/simulations/results/ReadDepth/amplipiper_test_0-${DIV}_${i}_${Freq}/analysis_97.5/results/consensus_seqs/amplipiper_test_0-${DIV}/COX1/COX1_consensussequences.fasta \
                --sanger ${WD}/simulations/data/Sy_tor_31_COI_short_${DIV}.fa \
                --primer GGTCAACAAATCATAAAGATATTGG,CCAAAAAATCARAATARRTGYTG \
                --frequency ${Freq} \
                --divergence ${DIV} \
                --SC 97.5 \
                --replicate ${i} >>${WD}/simulations/results/summary/SeqSimilarity_ReadDepth.csv

        done
    done
done

## plot
RScript -e """
library(tidyr)
library(tidyverse)
DATA<-read.csv('${WD}/simulations/results/summary/SeqSimilarity_ReadDepth.csv',header=T)


DATA.mean <- DATA %>%
    group_by(Divergence,ReadDepth) %>%
    summarise(Mean=mean(PropMM),SD=sd(PropMM))

# Define your three colors
three_colors <- c( '#FF8E3E','black')

# Create a discrete color palette with n levels (same as the number of categories)
n_levels <- length(levels(as.factor(DATA.mean$ReadDepth)))
my_palette <- colorRampPalette(three_colors)(n_levels)

PLOT<-ggplot(DATA.mean,aes(x=as.factor(Divergence),y=Mean,fill=as.factor(ReadDepth)))+
    geom_bar( stat='identity', position='dodge',alpha=0.7) +
    geom_errorbar( stat='identity', position='dodge', aes(
            ymin=Mean-SD, 
            ymax=Mean+SD),  
        colour='black', 
        alpha=0.9)+
    xlab('Divergence')+
    ylab('Proportion Mismatches')+
    scale_fill_manual(name='Read Depth',
    values=my_palette)+
    theme_bw() 
PLOT

ggsave('${WD}/simulations/results/summary/SeqSimilarity_ReadDepth.png',
    PLOT,
    width=10,
    height=2.5)

ggsave('${WD}/simulations/results/summary/SeqSimilarity_ReadDepth.pdf',
    PLOT,
    width=10,
    height=2.5)
"""

### Finally, compare COI consensus in real dataset to

for k in 0.00 0.03 0.05 0.08 0.1 0.2; do

    for mr in 100 200 500 1000; do

        bash /media/inter/pipelines/AmpliPiper/shell/AmpliPiper.sh \
            --samples /media/inter/pipelines/AmpliPiper/testdata/data/samples.csv \
            --primers /media/inter/pipelines/AmpliPiper/testdata/data/primers.csv \
            --output ${WD}/simulations/results/SyrBench/k${k}_nr${mr} \
            --quality 10 \
            --nreads ${mr} \
            --similar_consensus 97.5 \
            --threads 1 \
            --kthreshold ${k} \
            --minreads 50 \
            --sizerange 100 \
            --force &

    done

done

mkdir ${WD}/simulations/shell/qsub

printf "k\tm\tType\tID\tAbsDiff\tRelDiff\n" >${WD}/simulations/results/summary/Accuracy.txt

conda activate /media/inter/pipelines/AmpliPiper/envs/python_dependencies

for k in 0.00 0.03 0.05 0.08 0.1 0.2; do

    for mr in 100 200 500 1000; do

        python /media/inter/mkapun/projects/HAPLOTYPES-Benchmark/scripts/compare2Sanger.py \
            -i ${WD}/simulations/results/SyrBench/k${k}_nr${mr}/results/haplotypes/NCOX1/NCOX1.fasta \
            -s /media/inter/mkapun/projects/HAPLOTYPES-Benchmark/Syrphid/data/Syrphidae_COI_Sanger_red.fasta |
            awk -v k=${k} -v T="FlowCell" -v m=${mr} '{print k"\t"m"\t"T"\t"$0}' \
                >>${WD}/simulations/results/summary/Accuracy.txt

    done

done

conda deactivate

## plot
Rscript -e """
library(tidyr)
library(tidyverse)
DATA<-read.table('${WD}/simulations/results/summary/Accuracy.txt',header=T)

# Define your three colors
three_colors <- c('#21614A', '#FF8E3E','#7AC1DC')

# Create a discrete color palette with n levels (same as the number of categories)
n_levels <- length(levels(as.factor(DATA$k)))
my_palette <- colorRampPalette(three_colors)(n_levels)

PLOT<-ggplot(DATA,aes(x=ID,y=AbsDiff,fill=as.factor(k)))+
    geom_bar( stat='identity', position='dodge',alpha=0.7) +
    xlab('SampleID')+
    ylim(0,693)+
    ylab('No. Mismatches')+
    scale_fill_manual(name='k-threshold',
    values=my_palette)+
    theme_bw()+
    facet_grid(m~.,scales = 'free_y') + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

PLOT

ggsave('${WD}/simulations/results/summary/Accuracy.png',
    PLOT,
    width=10,
    height=6)

ggsave('${WD}/simulations/results/summary/Accuracy.pdf',
    PLOT,
    width=10,
    height=6)
"""
