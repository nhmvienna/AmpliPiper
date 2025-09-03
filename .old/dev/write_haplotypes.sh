#!/bin/bash

# Usage function
usage() {
    echo "Usage: bash ./write_haplotypes.sh -r, --ref REFERENCE_FASTA_FILE -o, --output OUTPUT_FILE [ -H, --haplotype <int: 1 or 2> -v, --vcf GZIPPED_VCF_FILE ]

  REQUIRED ARGUMENTS:
  -REFERENCE_FASTA_FILE: Provide the path to the reference fasta file, that should be already faidx-indexed with -H 2 (it is is not, write_haplotypes.sh takes charge of the indexing)
  -OUTPUT_FILE: Provide the path to the output file
  OPTIONAL ARGUMENTS:
  -GZIPPED_VCF_FILE: Provide the path to the vcf file, that should be compressed (.gz format, alternatively, write_haplotypes.sh takes charge of the compression) and already indexed (alternatively, write_haplotypes.sh takes charge of the indexing); if -H 2 is used, this argument is required!
  -HAPLOTYPE: Provide the haplotype as 1 (haploid) or 2 (diploid): default is 2
  
  Input bash ./write_haplotypes.sh -h or ./write_haplotypes.sh --help to show this message again"
    exit 1
}

# Initialize conda
Conda=$CONDA_PREFIX
source ${Conda}/etc/profile.d/conda.sh

# Initialize variables with default values
ref=""
haplotype=2
vcf=""
output=""

# Loop through the arguments
while [ $# -gt 0 ]; do
    case "$1" in
    -h | --help)
        usage
        ;;
    -v | --vcf)
        vcf="$2"
        shift 2
        ;;
    -r | --ref)
        ref="$2"
        shift 2
        ;;
    -H | --haplotype)
        haplotype="$2"
        shift 2
        ;;
    -o | --output) ######## <- new output folder
        output="$2"
        shift 2
        ;;
    *)
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
done

# Check if provided arguments are valid

## Return help message
if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    usage
fi

## Set haplotypes to 2 if not provided
if [ -z "${haplotype}" ] && [ ! "$1" = "-h" ] && [ ! "$1" = "--help" ]; then
    ${haplotype}=2
fi

## Check if ref is provided
if [ -z "${ref}" ] && [ ! "$1" = "-h" ] && [ ! "$1" = "--help" ]; then
    echo "Missing required argument: REFERENCE_FASTA_FILE"
    usage
fi

## Check if output is provided
if [ -z "${output}" ] && [ ! "$1" = "-h" ] && [ ! "$1" = "--help" ]; then
    echo "Missing required argument: OUTPUT_FILE"
    usage
fi

## Check if vcf is provided
if [ $haplotype -eq 2 ] && [ -z "${vcf}" ] && [ ! "$1" = "-h" ] && [ ! "$1" = "--help" ]; then
    echo "Missing required argument: GZIPPED_VCF_FILE"
    usage
fi

## If vcf and ref are provided but one of them or both do not exist, break and return usage
if [ $haplotype -eq 2 ] && [ ! -f "${vcf}" ] && [ ! -z "${vcf}" ] || [ $haplotype -eq 2 ] && [ ! -f "${ref}" ] && [ ! -z "${ref}" ]; then
    echo "Make sure that provided files exist and/or are actually files"
    usage
fi

## Index ref file if not indexed
if [ $haplotype -eq 2 ] && [ ! -z "${ref}" ] && [ -f "${ref}" ] && [ ! -f "${ref}.fai" ]; then
    echo "$ref should be provided in a .gz format, proceeding with compression..."
    conda activate samtools
    samtools faidx $ref
    conda deactivate
fi

## Compress vcf if not compressed
if [ $haplotype -eq 2 ] && [ ! -z "${vcf}" ] && [ -f "${vcf}" ] && [[ ! "${vcf}" == *.gz ]]; then
    echo "$vcf should be provided in a .gz format, proceeding with compression..."
    bgzip $vcf
    vcf=${vcf}.gz
fi

## Index vcf if not indexed
if [ $haplotype -eq 2 ] && [ ! -z "${vcf}" ] && [ -f "${vcf}" ] && [[ "${vcf}" == *.gz ]] && [ ! -f "${vcf}.tbi" ]; then
    echo "$vcf should be tabix-indexed, proceeding with indexing..."
    conda activate bcftools ### NEEDS TO BE INCLUDED IN setup.sh !!!
    bcftools tabix -p vcf $vcf
    conda deactivate
fi

# Main program
if [ ! -z "${vcf}" ] && [ -f "${vcf}" ] && [[ "${vcf}" == *.gz ]] && [ -f "${vcf}.tbi" ] && [ $haplotype -eq 2 ]; then
    dir=$(dirname $vcf)
    mkdir -p ${dir}/haplotypetmp ## create temporary directory
    conda activate bcftools
    for i in 1 2; do
        echo ">haplotype_${i}" >${dir}/haplotypetmp/header.txt                                                   ## Generate new header
        bcftools consensus -H ${i}pIu -f $ref $vcf >${dir}/haplotypetmp/tmpseq.fasta                             ## Reconstruct consensus from vcf
        tail +2 ${dir}/haplotypetmp/tmpseq.fasta >${dir}/haplotypetmp/seq.txt                                    ## Isolate sequence
        cat ${dir}/haplotypetmp/header.txt ${dir}/haplotypetmp/seq.txt >${dir}/haplotypetmp/haplotype_${i}.fasta ## Put header and sequence together
    done
    conda deactivate
    cat ${dir}/haplotypetmp/haplotype_1.fasta ${dir}/haplotypetmp/haplotype_2.fasta >$output ## Put haplotypes together
    rm ${dir}/haplotypetmp/*.*                                                               ## Remove everything from temporary directory
    rmdir ${dir}/haplotypetmp                                                                ## Remove temporary directory
elif [ $haplotype -eq 1 ]; then
    cat $ref >$output
else
    echo "HAPLOTYPE should be 1 or 2, not $haplotype"
fi
