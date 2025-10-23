#!/bin/bash

#===============================================================================
# AmpliPiper - Amplicon Sequencing Analysis Pipeline v2.0
#===============================================================================
# 
# DESCRIPTION:
#   AmpliPiper is a comprehensive pipeline for processing amplicon sequencing data
#   from raw FASTQ files to phylogenetic trees and species identification.
#   Version 2.0 includes outlier taxa filtering for improved phylogenetic analysis.
#
# WORKFLOW:
#   1. Quality filtering and demultiplexing by locus
#   2. Consensus sequence reconstruction with AmpliconSorter
#   3. Multiple sequence alignment with MAFFT
#   4. Genetic distance calculation and outlier detection
#   5. Phylogenetic reconstruction with IQ-TREE (outlier-filtered)
#   6. Species identification via BOLD/BLAST APIs
#   7. Species delimitation with ASAP (outlier-filtered)
#   8. HTML report generation
#
# DEPENDENCIES:
#   - Conda environments (automatically installed if missing)
#   - Python dependencies, R packages, bioinformatics tools
#
# AUTHORS: Astra Bertelli, Sonja Steindl & Martin Kapun
# VERSION: 2.0
# LICENSE: GPL3
#===============================================================================

## Usage function
usage() {
    echo "Usage: AmpliPiper -s, --samples SAMPLES_CSV -p, --primers PRIMERS_CSV -o,--output OUTPUT_FOLDER  [-q, --quality QUALITY] 
    [-n, --nreads NUMBER_OF_READS] [-t, --threads NUMBER_OF_THREADS] [-f,--force] [-b,--blast]

    AmpliPiper is a comprehensive pipeline for amplicon sequencing analysis, from raw 
    FASTQ files to phylogenetic trees and species identification.

    REQUIRED ARGUMENTS:
    -o | --output        : Path to the output folder where all results will be stored
    -p | --primers       : CSV file containing primer information (ID, forward_seq, reverse_seq, expected_size)
                          Example format: ID,Forward,Reverse,Size
                                         COX1,TTGATTTTTTGGTCATCCAGA,TAAACTTCAGGGTGACCAAAA,658
    -s | --samples       : CSV file containing sample information (ID, fastq_file_path)
                          Example format: ID,Path
                                         Sample1,/path/to/sample1.fastq.gz
                                         Sample2,/path/to/sample2.fastq
    
    OPTIONAL ARGUMENTS:
    -b | --blast         : Email address for BLAST API access (enables BLAST instead of BOLD for species ID)
                          Example: --blast your.email@institution.edu (default: uses BOLD API)
    -c | --similar_consensus : Similarity threshold (%) for AmpliconSorter consensus clustering
                          Lower values = more clusters, higher values = fewer clusters (default: 96)
    -d | --sd-distance    : Set the standard deviation genetic distance for outlier detection
                          Default is 3.0, adjust based on expected genetic variation
    -e | --exclude       : Text file listing sample-locus combinations to exclude from analysis
                          Format: one line per exclusion as SampleID,LocusID
    -f | --force         : Force overwrite of existing output directory (default: abort if exists)
    -g | --outgroup      : Sample ID(s) to use as outgroup for phylogenetic trees
                          Multiple outgroups: --outgroup Sample1,Sample2,Sample3
    -i | --partition     : Use partitioned substitution models for concatenated IQ-TREE analysis
                          Improves accuracy when combining multiple loci (default: single model)
    -k | --kthreshold    : Maximum mismatch rate for primer alignment during demultiplexing
                          Range 0.0-1.0, higher values = more permissive matching (default: 0.05)
    -m | --minreads      : Minimum number of reads required for consensus sequence reconstruction
                          Loci with fewer reads are excluded from analysis (default: 100)
    -n | --nreads        : Maximum number of top-quality reads used for consensus generation
                          Can be absolute number (e.g., 500) or percentage (e.g., 80%) (default: 500)
    -q | --quality       : Minimum Phred base quality score for read filtering
                          Reads with lower quality are discarded (default: 7)
    -r | --sizerange     : Allowed size deviation (bp) from expected amplicon length
                          Amplicons outside expected_size ± sizerange are discarded (default: 100)
    -t | --threads       : Number of CPU threads for parallel processing (default: 10)
    -w | --nowatermark   : Remove AmpliPiper watermark from phylogenetic tree figures
    -y | --freqthreshold : Minimum frequency threshold for consensus sequence retention
                          Range 0.0-1.0, sequences below this frequency are filtered out (default: 0.1)

    EXAMPLES:
    # Basic run with required parameters only:
    AmpliPiper -s samples.csv -p primers.csv -o results_folder

    # Advanced run with custom parameters:
    AmpliPiper -s samples.csv -p primers.csv -o results_folder \\
              -q 15 -n 1000 -t 20 -c 95 -f --blast your.email@uni.edu

    # Run with outgroup and partitioned analysis:
    AmpliPiper -s samples.csv -p primers.csv -o results_folder \\
              -g Reference_sample -i -y 0.05

  
  Input AmpliPiper -h,--help to show this help message"
    exit 1
}

#===============================================================================
# PARAMETER INITIALIZATION AND VALIDATION
#===============================================================================

## Initialize variables with default values
quality=7           # -q: Minimum base quality for filtering
similarconsensus=96 # -c: Similarity threshold for consensus clustering (%)
sdthreshold=3.0     # -d: Standard deviation threshold for outlier detection
nreads=500          # -n: Number of top quality reads for consensus
sizerange=100       # -r: Allowed size buffer around expected fragment length
minreads=100        # -m: Minimum reads required for consensus reconstruction
threads=10          # -t: Number of parallel threads
kthres=0.05         # -k: Maximum mismatches for primer alignment
force="no"          # -f: Force overwrite existing output
blast="no"          # -b: Use BLAST instead of BOLD for species ID
partition="no"      # -i: Use partition model for concatenated phylogeny
outgroup="no"       # -g: Outgroup sample(s) for phylogenetic analysis
freqthreshold=0.1   # -y: Minimum frequency threshold for consensus sequences

os="$(uname -s)"

#===============================================================================
# DIRECTORY SETUP AND DEPENDENCY CHECKS
#===============================================================================

## Find basedir - determine the base directory of the pipeline
tmp=$(dirname $0)
wd=${tmp%/*} ######## <- wd is now the base directory below the shell/ folder

## Clean up previous log directory and create fresh one

if
    [[ -d $wd/logs ]] \
        ;
then
    rm -rf $wd/logs
    mkdir $wd/logs
fi

## Test if software already installed, if not start installation
if
    [[ ! -d $wd/envs ]] \
        ;
then
    bash $wd/shell/setup.sh
fi

## Test if dependencies correctly installed - check for installation errors
if
    [[ -f $wd/logs/dep.err ]] \
        ;
then
    printf "\n######################\nError during installation:\n\n"

    # print errors, if there are any and quit
    while read -r line; do
        echo $line
    done <$wd/logs/dep.err
    printf "######################\n\n\n"
    usage
    exit 1
fi

## Test if software correctly installed - check for software package errors
if
    [[ -f $wd/envs/logs/setup.err ]] \
        ;
then
    printf "\n######################\nError during installation of the following software packages:\n\n"

    # print errors, if there are any and quit
    while read -r line; do
        echo $line
    done <$wd/envs/logs/setup.err
    printf "######################\n\n\n"
    usage
    exit 1
fi

## Initialize conda environment for loading software packages
eval "$(conda shell.bash hook)"

#===============================================================================
# COMMAND LINE ARGUMENT PARSING
#===============================================================================
# Parse all command line options and store in variables for later use
while
    [[ $# -gt 0 ]] \
        ;
do
    case "$1" in
    -h | --help)
        usage
        ;;
    -b | --blast)
        blast="$2"  # Email address for BLAST API access
        shift 2
        ;;
    -c | --similar_consensus)
        similarconsensus="$2"  # Similarity threshold for consensus clustering
        shift 2
        ;;
    -d | --sd-distance)
        sdthreshold="$2"  # Standard deviation threshold for outlier detection
        shift 2
        ;;
    -e | --exclude)
        exclude="$2"  # File with sample-locus combinations to exclude
        shift 2
        ;;
    -f | --force)
        force="yes"  # Force overwrite existing output directory
        shift 1
        ;;
    -g | --outgroup)
        outgroup="$2"  # Outgroup sample(s) for phylogenetic rooting
        shift 2
        ;;
    -i | --partition)
        partition="yes"  # Use partitioned model for concatenated analysis
        shift 1
        ;;
    -k | --kthreshold)
        kthres="$2"  # Primer alignment mismatch threshold
        shift 2
        ;;
    -m | --minreads)
        minreads="$2"  # Minimum reads required for consensus
        shift 2
        ;;
    -n | --nreads)
        nreads="$2"  # Number of reads for consensus generation
        shift 2
        ;;
    -o | --output)
        output="$2"  # Output directory path
        shift 2
        ;;
    -p | --primers)
        primers="$2"  # CSV file with primer information
        shift 2
        ;;
    -q | --quality)
        quality="$2"  # Minimum base quality threshold
        shift 2
        ;;
    -r | --sizerange)
        sizerange="$2"  # Allowed size deviation from expected fragment length
        shift 2
        ;;
    -s | --samples)
        samples="$2"  # CSV file with sample information
        shift 2
        ;;
    -t | --threads)
        threads="$2"  # Number of parallel processing threads
        shift 2
        ;;
    -y | --freqthreshold)
        freqthreshold="$2"  # Frequency threshold for consensus sequences
        shift 2
        ;;
    -w | --nowatermark)
        nowatermark="yes"  # Remove watermark from phylogenetic trees
        shift 2
        ;;
    *)
        echo "Unknown option: $1"
        usage
        ;;
    esac
done

#===============================================================================
# INPUT VALIDATION AND PREPROCESSING
#===============================================================================

## Validate that required sample input file is provided
if
    [[ -z "${samples}" ]] \
        ;
then
    printf "\n######################\nMissing required argument: INPUT FILE with names and paths to raw data in csv format\n######################\n\n\n"
    usage
fi

## Validate that required primer input file is provided
if
    [[ -z "${primers}" ]] \
        ;
then
    printf "\n######################\nMissing required argument: PRIMERS_TABLE\n######################\n\n\n"
    usage
fi

## Validate that required output directory path is defined
if
    [[ -z "${output}" ]] \
        ;
then
    printf "\n######################\nMissing required argument: path to OUTPUT folder\n######################\n\n\n"
    usage
fi

## Handle force overwrite option for pre-existing output folder
if
    [[ -d "${output}" ]] &&
        [[ ${force} == "yes" ]] \
        ;
then
    rm -rf ${output}  # Remove existing directory and its contents
    rm -rf ${output}  # Double removal to ensure cleanup

fi

## Check if output directory exists and force flag not set
if
    [[ -d "${output}" ]] \
        ;
then
    printf "\n######################\nOutput directory already exists, use -f to override\n######################\n\n\n"
    usage
fi

## Configure watermark setting for phylogenetic tree plots
if
    [[ ${nowatermark} = "yes" ]] \
        ;
then
    WM="NO"  # Disable watermark
else
    WM="YES"  # Enable watermark (default)
fi

## Store parameter settings as comma-separated string for later use
parametersettings="${quality},${similarconsensus},${nreads},${sizerange},${minreads},${threads},${kthres},${force},${blast},${partition},${outgroup}"

## Preprocess input files to ensure correct formatting
# Convert Windows line endings (\r\n) to Unix format (\n)
${wd}/envs/python_dependencies/bin/sed -i 's/\r$//' ${samples}
${wd}/envs/python_dependencies/bin/sed -i 's/\r$//' ${primers}

# Remove any spaces that might cause parsing issues
${wd}/envs/python_dependencies/bin/sed -i 's/ //g' ${samples}
${wd}/envs/python_dependencies/bin/sed -i 's/ //g' ${primers}

# Remove empty lines that might interfere with processing
${wd}/envs/python_dependencies/bin/sed -i '/^$/d' ${samples}
${wd}/envs/python_dependencies/bin/sed -i '/^$/d' ${primers}

# Validate that all FASTQ input files specified in samples file actually exist
while IFS=$"," read -r samplename file; do

    ## Skip header row
    if
        [[ ${samplename} == "ID" ]] \
            ;
    then
        continue
    fi

    ## Check if the specified FASTQ file exists
    if
        [[ ! -f ${file} ]] \
            ;
    then

        printf "\n######################\n${file} does not exist, quitting\n######################\n\n\n"
        usage
    fi

done <${samples}

## Parse exclusion file and fill array for samples/loci to exclude from analysis
declare -a EXCLUDE=()
if
    [[ ! -z "${exclude}" ]] \
        ;
then
    # Read exclude file and create array of sample-locus combinations
    mapfile -t EXCLUDE < <(awk -F',' '{print $1$2}' "${exclude}")
fi

## Set Search Engine variable for species identification
if [[ ${blast} != "no" ]]; then
    SE="BLAST"  # Use BLAST API for species identification
else
    SE="BOLD"   # Use BOLD API for species identification (default)
fi

#===============================================================================
# HELPER FUNCTION: FILTER OUTLIER TAXA
#===============================================================================
# Filter out outlier taxa from FASTA files based on genetic distance analysis
# Creates filtered versions of alignment files for downstream phylogenetic analysis

filter_outlier_taxa() {
    local primername=$1
    local input_fasta=$2
    local output_fasta=$3
    local outlier_file="${output}/results/haplotypes/outlier_taxa.csv"
    
    # Check if outlier file exists
    if [[ ! -f "${outlier_file}" ]]; then
        echo "Warning: Outlier file ${outlier_file} not found. Using original FASTA."
        cp "${input_fasta}" "${output_fasta}"
        return 0
    fi
    
    # Get outlier taxa for this locus (skip header, extract taxa for this primername)
    outlier_taxa=$(awk -F',' -v locus="${primername}" '
        NR > 1 && $1 == locus { print $2 }
    ' "${outlier_file}")
    
    if [[ -z "${outlier_taxa}" ]]; then
        echo "No outlier taxa found for locus ${primername}. Using original FASTA."
        cp "${input_fasta}" "${output_fasta}"
        return 0
    fi
    
    echo "Filtering outlier taxa for locus ${primername}: $(echo ${outlier_taxa} | tr '\n' ' ')"
    
    # Create a temporary file with outlier taxa (one per line)
    temp_outliers=$(mktemp)
    echo "${outlier_taxa}" | tr ' ' '\n' > "${temp_outliers}"
    
    # Filter the FASTA file using awk
    awk -v outliers="${temp_outliers}" '
    BEGIN {
        # Read outlier taxa into array
        while ((getline outlier < outliers) > 0) {
            outlier_list[outlier] = 1
        }
        close(outliers)
        keep_seq = 1
    }
    /^>/ {
        # Extract sequence name (remove > and everything after first space)
        seq_name = substr($1, 2)
        gsub(/[ \t].*/, "", seq_name)
        
        if (seq_name in outlier_list) {
            keep_seq = 0
            print "Excluding outlier taxon:", seq_name > "/dev/stderr"
        } else {
            keep_seq = 1
            print $0
        }
        next
    }
    keep_seq { print $0 }
    ' "${input_fasta}" > "${output_fasta}"
    
    # Clean up temporary file
    rm -f "${temp_outliers}"
    
    # Count sequences before and after filtering
    orig_count=$(grep -c "^>" "${input_fasta}")
    filt_count=$(grep -c "^>" "${output_fasta}")
    
    echo "Filtered ${primername}: ${orig_count} -> ${filt_count} sequences (removed $((orig_count - filt_count)) outliers)"
}

#===============================================================================
# PIPELINE EXECUTION BEGINS
#===============================================================================
now=$(date +%Y-%m-%d' '%H:%M:%S)
echo "+++++++++ Program started at ${now} +++++++++"

#===============================================================================
# STEP 1: PRIMER SIMILARITY ANALYSIS
#===============================================================================
# Analyze primer sequences for potential cross-reactivity and similarity

## Test primer sequence similarity
echo "***** Test primer sequence similarity *****"

## Create output folder structure for the entire pipeline
mkdir -p ${output}/results/summary/primers
mkdir -p ${output}/data/raw
mkdir -p ${output}/log/demulti
mkdir ${output}/data/demultiplexed
mkdir -p ${output}/shell/demult1
mkdir -p ${output}/log/summary
mkdir ${output}/data/filtered

# Activate Python environment and run primer comparison analysis
conda activate ${wd}/envs/python_dependencies

${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/CompPrimers.py \
    -p ${primers} \
    -o ${output}/results/summary/primers/primers_dist.csv \
    >${output}/results/summary/min_edit.dist 2>${output}/log/summary/primerdist.log

conda deactivate

echo "finished"

## Count total number of loci specified in primers file for downstream processing
LOCI=$(awk '!/^ID,/' ${primers} | wc -l)

#===============================================================================
# STEP 2: QUALITY FILTERING AND DEMULTIPLEXING
#===============================================================================
# Process raw FASTQ files through three main stages:
# 1. Copy raw FASTQ files to output directory (compress if needed)
# 2. Filter reads by minimum base quality using Chopper
# 3. Demultiplex reads by locus using primer sequences with fuzzy matching
#
# This step generates individual shell scripts for each sample to enable 
# parallel processing across multiple samples simultaneously.
# The advanced demultiplexing script allows for:
# - Identify potential chimeras consisting in multiple amplicons in one read 
# - First-amplicon-only retention to avoid chimeras
# - Size range filtering based on expected amplicon length 

## Process all input files: copy, filter, and demultiplex by locus
echo "***** Copying files, starting filtering and demultiplexing by locus *****"

## Generate individual shell scripts for each sample to enable parallel processing
# Each script handles: file copying/compression, quality filtering, and demultiplexing

while IFS=$"," read -r samplename file; do

    ## skip header
    if
        [[ ${samplename} == "ID" ]] \
            ;
    then
        continue
    fi

    echo """
    
    eval \"\$(conda shell.bash hook)\"

    ## Copy raw data to output folder, compressing if necessary

    if [[ ${file} == *.gz ]]; then
        # File already compressed, copy directly
        cp -n ${file} ${output}/data/raw/${samplename}.fastq.gz
    else
        # Compress uncompressed FASTQ file during copy
        conda activate ${wd}/envs/chopper
        pigz -c ${file} > ${output}/data/raw/${samplename}.fastq.gz
        conda deactivate
    fi

    ## filter raw sequences with chopper
    conda activate ${wd}/envs/chopper
    pigz -dc ${output}/data/raw/${samplename}.fastq.gz |
        chopper -q ${quality} 2>> ${output}/log/demulti/${samplename}_demulti.log | pigz \
        >${output}/data/filtered/${samplename}-filt.fastq.gz
        
    conda deactivate

    ## demultiplex per locus and sample
    mkdir ${output}/data/demultiplexed/${samplename}

    conda activate ${wd}/envs/python_dependencies
    ${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/DemultFastqAdvanced.py \
        -i ${output}/data/filtered/${samplename}-filt.fastq.gz \
        -p $primers \
        -o ${output}/data/demultiplexed/${samplename} \
        -th ${kthres} \
        -sr ${sizerange} \
        --first-amplicon-only \
        -mr ${minreads} \
        -rp ${nreads} >> ${output}/log/demulti/${samplename}_demulti.log 2>&1
    conda deactivate

    #echo '${samplename} finished'

    """ >${output}/shell/demult1/${samplename}.sh

done <${samples}

## Execute shell scripts in parallel
conda activate ${wd}/envs/parallel
echo "will cite" | parallel --citation >/dev/null 2>&1
parallel --bar -j${threads} bash ::: ${output}/shell/demult1/*.sh
conda deactivate

## Consensus sequence reconstruction
echo "***** Files demultiplexed, starting consensus haplotype reconstruction *****"

#===============================================================================
# STEP 3: CONSENSUS SEQUENCE RECONSTRUCTION WITH AMPLICONSORTER
#===============================================================================
# Generate consensus sequences from demultiplexed reads using AmpliconSorter:
# - Groups similar reads into clusters based on similarity threshold
# - Reconstructs consensus sequences for each cluster
# - Filters clusters by minimum read count and frequency thresholds
# - Outputs FASTA files with consensus sequences for downstream analysis
#
# Each sample-locus combination is processed independently in parallel.

mkdir -p ${output}/results/consensus_seqs
mkdir -p ${output}/log/ampliconsorter
mkdir -p ${output}/shell/demult2

while IFS=$"," read -r samplename file; do

    ## skip header
    if [[ ${samplename} == "ID" ]]; then
        continue
    fi

    while IFS=$"," read -r primername fwd rev size; do

        ## skip samples/loci in exclude list
        if [[ " ${EXCLUDE[*]} " == *" $samplename$primername "* ]]; then
            echo "skipping "$samplename": "$primername
            continue
        fi

        ## skip header
        if
            [[ ${primername} == "ID" ]] \
                ;
        then
            continue
        fi

        echo """

        ## source conda
        eval \"\$(conda shell.bash hook)\"

        ## make shortcut for path and create directories
        SamPrim=${output}/results/consensus_seqs/${samplename}/${primername}
        mkdir -p ${output}/results/consensus_seqs/${samplename}
        cd ${output}/results/consensus_seqs/${samplename}

        ## AmpliconSorter consensus sequence reconstruction
        conda activate ${wd}/envs/python_dependencies
        ## Adjust AmpliconSorter parameters for macOS compatibility
        if
            [[ ${os} == "Darwin" ]] \
                ;
        then
            # macOS version with -mac flag for compatibility
            amplicon_sorter.py \
                -i ${output}/data/demultiplexed/${samplename}/${primername}.fastq \
                -np 1 \
                -mac \
                --similar_consensus ${similarconsensus} \
                -maxr ${nreads} \
                -o \${SamPrim} \
            >> ${output}/log/ampliconsorter/${samplename}_${primername}_AS.log 2>&1
        else
            # Linux/Unix version without -mac flag
            amplicon_sorter.py \
                -i ${output}/data/demultiplexed/${samplename}/${primername}.fastq \
                -np 1 \
                --similar_consensus ${similarconsensus} \
                -maxr ${nreads} \
                -o \${SamPrim} \
            >> ${output}/log/ampliconsorter/${samplename}_${primername}_AS.log 2>&1
        fi
        """ >${output}/shell/demult2/${samplename}_${primername}.sh

    done <${primers}

done <${samples}

## Execute shell scripts in parallel
conda activate ${wd}/envs/parallel
echo "will cite" | parallel --citation >/dev/null 2>&1
parallel --bar -j${threads} bash ::: ${output}/shell/demult2/*.sh
conda deactivate

#===============================================================================
# STEP 4: SUMMARY STATISTICS AND HAPLOTYPE SELECTION
#===============================================================================
# Process AmpliconSorter outputs to generate summary statistics and select
# representative haplotypes for downstream analysis:
# 1. Parse AmpliconSorter output files to extract read counts and frequencies
# 2. Generate summary CSV with statistics for all sample-locus combinations
# 3. Select consensus sequences meeting frequency thresholds
# 4. Create visualization of missing data patterns across samples and loci

## make CSV summary
echo "***** Summarize Ampliconsorter output *****"

${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/ParseSummary.py \
    --path ${output} \
    --primer ${primers} \
    --samples ${samples} \
    >${output}/results/summary/summary.csv

conda activate ${wd}/envs/python_dependencies

## Select consensus sequences meeting frequency and read count thresholds
${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/ChooseCons.py \
    --input ${output}/results/summary/summary.csv \
    --path ${output}/results/consensus_seqs \
    --output ${output}/results/haplotypes \
    --FreqTH ${freqthreshold} \
    >>${output}/log/ampliconsorter/Summary.log 2>&1

conda deactivate

## Update summary file with ploidy information 
mv ${output}/results/summary/summary.csv.ploidy ${output}/results/summary/summary.csv

## Generate heatmap visualization of missing data patterns across samples and loci
conda activate ${wd}/envs/R
${wd}/envs/R/bin/Rscript ${wd}/scripts/MissingDataHeatmap.r \
    ${output}/results/summary/summary.csv \
    >>${output}/log/ampliconsorter/Summary_AS.log 2>&1
conda deactivate

echo " finished"

#===============================================================================
# STEP 5: MULTIPLE SEQUENCE ALIGNMENT
#===============================================================================
# Align consensus sequences for each locus using MAFFT:
# - Uses accurate alignment mode with direction adjustment
# - Processes each locus independently 
# - Requires minimum of 4 sequences per locus for meaningful alignment
# - Outputs aligned FASTA files for phylogenetic analysis

## align haplotypes
echo "***** align haplotypes *****"

mkdir -p ${output}/log/tree

while IFS=$"," read -r primername fwd rev size; do

    ## skip header
    if
        [[ ${primername} == "ID" ]] \
            ;
    then
        continue
    fi

    ## Skip if FASTA file doesn't exist for this locus
    if
        [[ ! -f ${output}/results/haplotypes/${primername}/${primername}.fasta ]] \
            ;
    then
        continue
    fi

    ## Skip if insufficient sequences for meaningful phylogenetic analysis (need ≥4)
    if
        [[ $(grep "^>" ${output}/results/haplotypes/${primername}/${primername}.fasta | wc -l) -lt 4 ]] \
            ;
    then
        continue
    fi

    ## Perform multiple sequence alignment using MAFFT with direction adjustment
    conda activate ${wd}/envs/mafft

    mafft \
        --adjustdirectionaccurately \
        ${output}/results/haplotypes/${primername}/${primername}.fasta 2>>${output}/log/tree/${primername}_TREE.log |
        ${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/fixIDAfterMafft.py \
            --input ${output}/results/haplotypes/${primername}/${primername}.fasta \
            --Alignment - \
            >${output}/results/haplotypes/${primername}/${primername}_aln.fasta

    conda deactivate

    echo "locus ${primername} finished"

done <${primers}

#===============================================================================
# STEP 6: GENETIC DISTANCE CALCULATIONS
#===============================================================================
# Calculate pairwise genetic distances between haplotypes using R:
# - Computes sequence divergence metrics for each aligned locus
# - Generates distance matrices for phylogenetic analysis
# - Identifies outlier taxa for exclusion from downstream analyses
# - Outputs summary statistics of genetic diversity

echo "***** Calculate Genetic Distances *****"

conda activate ${wd}/envs/R

${wd}/envs/R/bin/Rscript ${wd}/scripts/GeneticDist_v2.r \
    ${output}/results/haplotypes ${outgroup} ${sdthreshold} \
    >>${output}/log/ampliconsorter/GeneticDistance_AS.log 2>&1

echo "finished"

conda deactivate

#===============================================================================
# STEP 7: SPECIES IDENTIFICATION
#===============================================================================
# Identify species for standard barcoding loci (COX1, ITS, MATK_RBCL) using:
# - BOLD API (default): Queries Barcode of Life database
# - BLAST API (optional): Queries NCBI GenBank database
# 
# Only processes recognized barcoding markers that are present in the dataset.
# Results include taxonomic assignments with confidence scores and metadata.

## Species identification with BOLD/BLAST
mkdir ${output}/log/SpecID
PRINT=0

while IFS=$"," read -r primername fwd rev size; do

    if
        [[ ${primername} == "COX1" ]] ||
            [[ ${primername} == "ITS" ]] ||
            [[ ${primername} == "MATK_RBCL" ]] \
            ;
    then

        ## Print Header
        if
            [[ ${PRINT} == 0 ]] \
                ;
        then
            echo "***** Species ID from ${SE} *****"
            PRINT=1
        fi

        ## use BOLD/BLAST API for Species identification
        conda activate ${wd}/envs/python_dependencies

        if [[ ${blast} != "no" ]]; then
            mkdir -p ${output}/results/SpeciesID/${SE}/${primername}/summarized_outputs
            ${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/BLASTapi.py \
                -i ${output}/results/haplotypes/${primername}/${primername}_aln.fasta \
                -e ${blast} \
                -o ${output}/results/SpeciesID/${SE}/${primername}/summarized_outputs \
                >>${output}/log/SpecID/${primername}_SI.log 2>&1
        else
            mkdir -p ${output}/results/SpeciesID/${SE}/${primername}
            ${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/bold_api/BOLDapi.py \
                -i ${output}/results/haplotypes/${primername}/${primername}_aln.fasta \
                -p ${primername} \
                -c ${wd}/scripts/style.css \
                -n 10 \
                -o ${output}/results/SpeciesID/${SE}/${primername} \
                >>${output}/log/SpecID/${primername}_SI.log 2>&1
        fi
        conda deactivate

        echo ${primername} "finished"

    fi

done <${primers}

## Determine the best locus for species identification (preferentially COX1)
# Search for species identification results in order of preference: COX1 > ITS > MATK_RBCL
for locus in COX1 ITS MATK_RBCL; do
    if
        [[ -f ${output}/results/SpeciesID/${SE}/${locus}/summarized_outputs/final.csv ]] \
            ;
    then
        ID=${locus}  # Store the locus with available species IDs for tree labeling
        break
    fi
done

#===============================================================================
# STEP 8: MULTI-LOCUS CONCATENATION
#===============================================================================
# Concatenate alignments from multiple loci for combined phylogenetic analysis:
# - Only performed when multiple loci are present and frequency threshold > 0
# - Creates concatenated alignment with partition information
# - Enables analysis of phylogenetic signal across multiple markers

if [[ ${LOCI} -gt 1 && $(ls -l ${output}/results/haplotypes/*/*_aln.fasta | wc -l) -gt 1 && ${freqthreshold} != "0" && ${freqthreshold} != "0.0" && ${freqthreshold} != "0.00" ]]; then

    ## Species delineation with ASAP for all concatenated Haplotypes
    echo "***** concatenate all loci *****"

    conda activate ${wd}/envs/asap
    mkdir -p ${output}/results/haplotypes/Concatenated_loci

    ${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/MergeAln.py \
        --input "${output}/results/haplotypes/*/*_aln.fasta" \
        --output ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci

    conda deactivate
    echo "finished"
fi

#===============================================================================
# STEP 9: PHYLOGENETIC RECONSTRUCTION WITH OUTLIER FILTERING
#===============================================================================
# Reconstruct maximum likelihood phylogenetic trees using IQ-TREE:
# - Individual trees for each locus (requires ≥4 sequences)
# - Optional concatenated tree across all loci
# - Bootstrap support values (1000 replicates)
# - Tree visualization with R/ggplot2
# - Integration of species names from BOLD/BLAST if available
# - Outlier taxa filtering based on genetic distance analysis

## reconstruct trees if > 3 haplotypes in input
echo "***** reconstruct ML trees per locus (with outlier filtering) *****"

mkdir -p ${output}/log/tree

## Adjust tree plot dimensions based on whether frequency information is displayed
# When freqthreshold=0, frequencies are shown in sequence names, requiring wider plots
if [[ ${freqthreshold} == 0 ]]; then
    WIDTH=10  # Wider plot to accommodate frequency information
else
    WIDTH=8   # Standard width for sequence names only
fi

while IFS=$"," read -r primername fwd rev size; do

    ## skip header
    if
        [[ ${primername} == "ID" ]] \
            ;
    then
        continue
    fi

    if
        [[ ! -f ${output}/results/haplotypes/${primername}/${primername}.fasta ]] \
            ;
    then
        continue
    fi

    ## skip if less than four aligned sequences
    if
        [[ $(grep "^>" ${output}/results/haplotypes/${primername}/${primername}.fasta | wc -l) -lt 4 ]] \
            ;
    then
        continue
    fi

    ## Filter outlier taxa from alignment before tree reconstruction
    mkdir -p ${output}/results/tree/${primername}/
    
    # Create filtered alignment excluding outlier taxa
    filter_outlier_taxa "${primername}" \
        "${output}/results/haplotypes/${primername}/${primername}_aln.fasta" \
        "${output}/results/tree/${primername}/${primername}_filtered.fasta"
    
    ## Check if we still have enough sequences after filtering
    if
        [[ $(grep "^>" ${output}/results/tree/${primername}/${primername}_filtered.fasta | wc -l) -lt 4 ]] \
            ;
    then
        echo "Warning: Less than 4 sequences remaining for ${primername} after outlier filtering. Skipping tree reconstruction."
        continue
    fi

    ## Maximum likelihood phylogenetic reconstruction using IQ-TREE with bootstrap support
    cp ${output}/results/tree/${primername}/${primername}_filtered.fasta \
        ${output}/results/tree/${primername}/${primername}

    conda activate ${wd}/envs/iqtree

    # Run IQ-TREE with automatic model selection and 1000 bootstrap replicates
    iqtree \
        -s ${output}/results/tree/${primername}/${primername} \
        -ntmax ${threads} \
        -B 1000 \
        >>${output}/log/tree/${primername}_TREE.log 2>&1

    conda deactivate

    conda activate ${wd}/envs/R

    ## Dynamically adjust tree plot height based on number of sequences
    if
        [[ $(($(grep "^>" ${output}/results/tree/${primername}/${primername} | wc -l) / 3)) -gt 8 ]] \
            ;
    then
        HEIGHT=$(($(grep "^>" ${output}/results/tree/${primername}/${primername} | wc -l) / 3))
    else
        HEIGHT=8  # Minimum height for readability
    fi

    ## Integrate species names from BOLD/BLAST results if available
    OFFSET=0.3  # Default x-axis offset for tip labels
    outgroupNew="no"
    if
        [[ ! -z ${ID} ]] \
            ;
    then
        OFFSET=0.7  # Increased offset to accommodate longer species names
        # Replace sequence IDs with species names in tree file
        outgroupNew=$(${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/RenameTreeLeaves.py \
            --primername ${primername} \
            --input ${output}/results/tree/${primername}/${primername}.treefile \
            --name ${output}/results/SpeciesID/${SE}/${ID}/summarized_outputs/final.csv \
            --outgroup ${outgroup})
    fi

    ## Generate publication-ready tree plots using R/ggplot2
    ${wd}/envs/R/bin/Rscript ${wd}/scripts/PlotTree.r \
        ${output}/results/tree/${primername}/${primername}.treefile \
        ${output}/results/tree/${primername}/${primername} \
        ${primername} \
        ${OFFSET} \
        ${WIDTH} \
        ${HEIGHT} \
        ${outgroupNew} \
        $WM \
        >>${output}/log/tree/${primername}_TREE.log 2>&1

    conda deactivate

    echo "ML tree for locus ${primername} finished"

done <${primers}

## Process concatenated alignment for phylogenetic reconstruction
# Concatenated analysis provides additional phylogenetic resolution by combining
# signal from multiple loci, with optional partitioned substitution models

if
    [[ -s "${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta" ]] \
        ;
then

    if
        [[ ${partition} == "yes" ]] \
            ;
    then
        ## Phylogeny using IQtree with 100 bootsrapping rounds
        mkdir -p ${output}/results/tree/Concatenated_loci/
        cp ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta \
            ${output}/results/tree/Concatenated_loci/Concatenated_loci.fasta
        cp ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.part \
            ${output}/results/tree/Concatenated_loci/Concatenated_loci

        conda activate ${wd}/envs/iqtree

        cd ${output}/results/haplotypes/Concatenated_loci

        iqtree \
            -s ${output}/results/tree/Concatenated_loci/Concatenated_loci.fasta \
            -ntmax ${threads} \
            -B 1000 \
            -p ${output}/results/tree/Concatenated_loci/Concatenated_loci \
            >>${output}/log/tree/Concatenated_loci_TREE.log 2>&1

        conda deactivate

    else

        ## Phylogeny using IQtree with 100 bootsrapping rounds
        mkdir -p ${output}/results/tree/Concatenated_loci/
        cp ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta \
            ${output}/results/tree/Concatenated_loci/Concatenated_loci

        conda activate ${wd}/envs/iqtree

        cd ${output}/results/haplotypes/Concatenated_loci

        iqtree \
            -s ${output}/results/tree/Concatenated_loci/Concatenated_loci \
            -ntmax ${threads} \
            -B 1000 \
            >>${output}/log/tree/Concatenated_loci_TREE.log 2>&1

        conda deactivate

    fi

    conda activate ${wd}/envs/R

    ## adjust tree height based on samples in dataset
    if
        [[ $(($(grep "^>" ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta | wc -l) / 3)) -gt 8 ]] \
            ;
    then
        HEIGHT=$(($(grep "^>" ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta | wc -l) / 3))
    else
        HEIGHT=8
    fi

    ## append BOLD names if available and adjust x-axis offset to account for longer names
    OFFSET=0.3
    outgroupNew="no"
    if
        [[ ! -z ${ID} ]] \
            ;
    then
        OFFSET=0.7
        outgroupNew=$(${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/RenameTreeLeaves.py \
            --primername combined \
            --input ${output}/results/tree/Concatenated_loci/Concatenated_loci.treefile \
            --name ${output}/results/SpeciesID/${SE}/${ID}/summarized_outputs/final.csv \
            --outgroup ${outgroup})

    fi

    echo ${outgroup} ${outgroupNew}
    ## plot trees with ggplot
    ${wd}/envs/R/bin/Rscript ${wd}/scripts/PlotTree.r \
        ${output}/results/tree/Concatenated_loci/Concatenated_loci.treefile \
        ${output}/results/tree/Concatenated_loci/Concatenated_loci \
        Concatenated_loci \
        ${OFFSET} \
        ${WIDTH} \
        ${HEIGHT} \
        ${outgroupNew} \
        ${WM} \
        >>${output}/log/tree/Concatenated_loci_TREE.log 2>&1

    conda deactivate

    echo "ML tree for Concatenated_loci dataset finished"

fi

#===============================================================================
# STEP 10: ASTRAL SPECIES TREE RECONSTRUCTION  
#===============================================================================
# Reconstruct species tree using ASTRAL (coalescent-based method):
# - Combines individual gene trees into a species tree
# - Accounts for incomplete lineage sorting
# - Only performed when multiple loci are available
# - Requires frequency threshold > 0 to ensure adequate data

if [[ ${LOCI} -gt 1 && ${freqthreshold} != "0" && ${freqthreshold} != "0.0" && ${freqthreshold} != "0.00" ]]; then
    ## Create ASTRAL consensus tree from individual gene trees
    mkdir -p ${output}/results/tree/ASTRAL

    ## Collect all individual gene trees for species tree reconstruction
    cat ${output}/results/tree/*/*.treefile >${output}/results/astral_input.tree
    mv ${output}/results/astral_input.tree ${output}/results/tree/ASTRAL

    ## Verify that gene trees are available for analysis
    if [[ ! -z $(grep '[^[:space:]]' ${output}/results/tree/ASTRAL/astral_input.tree) ]]; then

        echo "***** reconstruct ASTRAL tree across all loci *****"

        ## Calculate appropriate plot dimensions based on number of samples
        if
            [[ $(($(cat ${samples} | wc -l) / 3)) -gt 8 ]] \
                ;
        then
            HEIGHT= $(($(cat ${samples} | wc -l) / 3))
        else
            HEIGHT=8  # Minimum height for readability
        fi

        ## Set plot offset based on whether species names are available
        if
            [[ ! -z ${ID} ]] \
                ;
        then
            OFFSET=0.5  # Space for species names
        else
            OFFSET=0.1  # Minimal space for sequence IDs only
        fi

        ## Reconstruct coalescent-based species tree using WASTRAL
        conda activate ${wd}/envs/aster

        wastral \
            -i ${output}/results/tree/ASTRAL/astral_input.tree \
            -o ${output}/results/tree/ASTRAL/ASTRAL.tree \
            >>${output}/log/tree/Astral.log 2>&1

        conda deactivate

        ## Generate ASTRAL tree visualization
        conda activate ${wd}/envs/R

        ${wd}/envs/R/bin/Rscript ${wd}/scripts/PlotTree_astral.r \
            ${output}/results/tree/ASTRAL/ASTRAL.tree \
            ${output}/results/tree/ASTRAL/ASTRAL \
            ASTRAL \
            ${OFFSET} \
            ${WIDTH} \
            ${HEIGHT} \
            ${outgroupNew} \
            ${WM} \
            >>${output}/log/tree/Astral.log 2>&1

        conda deactivate

        echo "finished"

    fi
fi

#===============================================================================
# STEP 11: SPECIES DELIMITATION WITH ASAP (OUTLIER FILTERED)
#===============================================================================  
# Perform automatic species delimitation using ASAP (Assemble Species by 
# Automatic Partitioning):
# - Identifies potential species boundaries within datasets
# - Uses multiple genetic distance thresholds
# - Provides statistical support for delimitation hypotheses
# - Processes both individual loci and concatenated datasets
# - Excludes outlier taxa identified by genetic distance analysis

## Species delineation with ASAP
mkdir ${output}/log/SpecDelim

echo "***** Species delimitation with ASAP for each locus (with outlier filtering) *****"

conda activate ${wd}/envs/asap

while IFS=$"," read -r primername fwd rev size; do

    ## skip header
    if
        [[ ${primername} == "ID" ]] \
            ;
    then
        continue
    fi

    if
        [[ ! -f ${output}/results/haplotypes/${primername}/${primername}_aln.fasta ]] \
            ;
    then
        continue
    fi

    ## Create filtered alignment for ASAP analysis
    mkdir -p ${output}/results/SpeciesDelim/${primername}
    
    # Filter outlier taxa before ASAP analysis
    filter_outlier_taxa "${primername}" \
        "${output}/results/haplotypes/${primername}/${primername}_aln.fasta" \
        "${output}/results/SpeciesDelim/${primername}/${primername}_filtered.fasta"
    
    ## Check if we have enough sequences after filtering
    if
        [[ $(grep "^>" ${output}/results/SpeciesDelim/${primername}/${primername}_filtered.fasta | wc -l) -lt 3 ]] \
            ;
    then
        echo "Warning: Less than 3 sequences remaining for ${primername} after outlier filtering. Skipping ASAP analysis."
        continue
    fi
    
    cd ${output}/results/SpeciesDelim/${primername}

    asap \
        -a ${output}/results/SpeciesDelim/${primername}/${primername}_filtered.fasta \
        -o ${output}/results/SpeciesDelim/${primername} \
        >>${output}/log/SpecDelim/${primername}_SD.log 2>&1

    echo "locus ${primername} finished"

done <${primers}

## run ASAP on concatenated FASTA if file not empty
if
    [[ -s "${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta" ]] \
        ;
then

    mkdir ${output}/results/SpeciesDelim/Concatenated_loci

    asap \
        -a ${output}/results/haplotypes/Concatenated_loci/Concatenated_loci.fasta \
        -o ${output}/results/SpeciesDelim/Concatenated_loci \
        >>${output}/log/SpecDelim/Concatenated_loci_SD.log 2>&1
    
    conda deactivate

fi
conda deactivate

#===============================================================================
# STEP 12: HTML REPORT GENERATION
#===============================================================================
# Generate comprehensive HTML report containing all pipeline results:
# - Interactive multiple sequence alignments
# - Phylogenetic trees with species annotations
# - Species identification results
# - Species delimitation outcomes  
# - Summary statistics and data quality metrics
# - Pipeline parameter settings and metadata

## make HTML summary

mkdir -p ${output}/Output

echo "***** Summarize in HTML file *****"

conda activate ${wd}/envs/python_dependencies

mkdir ${output}/results/html
mkdir ${output}/log/html

## visualize MSA
unset SESSION_MANAGER

${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/msa_to_html.py \
    -hap ${output}/results/haplotypes/ \
    -py "${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/msaviz.py" \
    -hf ${output}/results/html/MSA_analysis.html \
    -sf ${wd}/scripts/scripts/style.css \
    >>${output}/log/html/MSA.log 2>&1

unset SESSION_MANAGER

samplestba=$(awk -F',' '{print$1}' "${samples}")
primerstba=$(awk -F',' 'NR>1 {print $1}' "${primers}")

## make complete HTML output

mkdir ${output}/.logos
cp ${wd}/imgs/tettris.png ${output}/.logos/tettris.png
cp ${wd}/imgs/nhm.svg.png ${output}/.logos/nhm.svg.png

${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/altersvg.py \
    -f ${output}/results/SpeciesDelim/ \
    >>${output}/log/html/HTML.log 2>&1

${wd}/envs/python_dependencies/bin/python3 ${wd}/scripts/DisplayOutput.py \
    -p ${parametersettings} \
    -r ${output}/results \
    -out ${output}/Output \
    -loci ${primerstba} \
    -samples ${samplestba} \
    >>${output}/log/html/HTML.log 2>&1

conda deactivate

## Copy species identification results to final output directory
if [[ -d ${output}/results/SpeciesID ]]; then
    for PA in ${output}/results/SpeciesID/${SE}/*; do
        IDlocus=${PA##*/}
        cp ${PA}/summarized_outputs/final.csv ${output}/Output/summary/SpeciesID_${IDlocus}.csv
    done
fi

## Copy final haplotype alignments to output directory for user access
while IFS=$"," read -r primername fwd rev size; do

    ## skip header
    if
        [[ ${primername} == "ID" ]] \
            ;
    then
        continue
    fi

    if
        [[ ! -f ${output}/results/haplotypes/${primername}/${primername}_aln.fasta ]] \
            ;
    then
        continue
    fi

    cp ${output}/results/haplotypes/${primername}/${primername}_aln.fasta \
        ${output}/Output/haplotypes/${primername}

done <${primers}

## Clean up temporary R plot files
rm -f ${output}/results/haplotypes/Rplots.pdf
rm -f ${output}/Output/haplotypes/Rplots.pdf

#===============================================================================
# PIPELINE COMPLETION
#===============================================================================
# The AmpliPiper pipeline has completed successfully!
# 
# MAIN OUTPUTS:
# - ${output}/Output/: Final results directory containing:
#   * index.html: Interactive HTML report with all results
#   * haplotypes/: Aligned consensus sequences for each locus  
#   * summary/: Summary statistics and species identification results
#   * trees/: Phylogenetic trees in multiple formats (PDF, SVG, Newick)
#   * species_delimitation/: ASAP species delimitation results
# 
# - ${output}/results/: Detailed intermediate results from each analysis step
# - ${output}/log/: Log files for troubleshooting and quality control
#
# NEW IN VERSION 2.0:
# - Outlier taxa filtering based on genetic distance analysis
# - Improved phylogenetic reconstruction excluding problematic sequences
# - Enhanced species delimitation with outlier-filtered datasets
#
# For questions or issues, please refer to the pipeline documentation or
# contact the development team.

## Finished
now=$(date +%Y-%m-%d' '%H:%M:%S)
echo " +++++++++ Program finished at ${now} +++++++++ "
