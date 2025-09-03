## Run testdataset

# (1) Define path to AmpliPiper folder
# replace '<path to AmpliPiper folder>' with real path
WD='<path to AmpliPiper folder>'
## e.g., WD=/media/inter/pipelines/AmpliPiper

# (2) generate samples.csv file

## print header for samples.csv
printf "ID,PATH\n" >${WD}/testdata/data/samples.csv

## loop through input FASTQ files
for Filepath in ${WD}/testdata/reads/*fastq.gz; do

    ## get Filename
    Filename=${Filepath##*/}

    ## get File ID
    ID=${Filename%.fastq.gz*}

    ## print to samples.csv
    echo ${ID},${Filepath} >>${WD}/testdata/data/samples.csv
done

# (3a) run AmpliPiper with default settings
bash ${WD}/shell/AmpliPiper.sh \
    --samples ${WD}/testdata/data/samples.csv \
    --primers ${WD}/testdata/data/primers.csv \
    --output ${WD}/testdata/results/demo_250819 \
    --blast your@email.com \
    --threads 150 \
    --outgroup He_mor_41 \
    --force

# (3b) run AmpliPiper with custom settings
bash ${WD}/shell/AmpliPiper_v2.sh \
    --samples ${WD}/testdata/data/samples.csv \
    --primers ${WD}/testdata/data/primers.csv \
    --output ${WD}/testdata/results/demo_250820_v2 \
    --quality 10 \
    --nreads 1000 \
    --blast your@email.com \
    --similar_consensus 97 \
    --threads 100 \
    --kthreshold 0.05 \
    --minreads 50 \
    --sizerange 100 \
    --outgroup He_mor_41 \
    --sd-distance 1.0 \
    --force

    MicrMu_IC76

bash ${WD}/shell/AmpliPiper_v2.sh \
    --samples ${WD}/testdata/data/samples.csv \
    --primers ${WD}/testdata/data/primers.csv \
    --output ${WD}/testdata/results/demo_250821_v2 \
    --quality 10 \
    --nreads 1000 \
    --blast your@email.com \
    --similar_consensus 97 \
    --threads 100 \
    --kthreshold 0.05 \
    --minreads 50 \
    --sizerange 100 \
    --outgroup He_mor_41,MicrMu_IC76,MicrDe_IC99 \
    --sd-distance 1.0 \
    --force

# (3b) run AmpliPiper and keep all consensus sequences

bash ${WD}/shell/AmpliPiper.sh \
    --samples ${WD}/testdata/data/samples.csv \
    --primers ${WD}/testdata/data/primers.csv \
    --output ${WD}/testdata/results/demo_allCons_250818 \
    --quality 10 \
    --nreads 1000 \
    --blast your@email.com \
    --similar_consensus 97 \
    --threads 100 \
    --kthreshold 0.05 \
    --minreads 50 \
    --sizerange 100 \
    --outgroup He_mor_41 \
    --freqthreshold 0 \
    --force
