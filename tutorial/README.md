<h1 align="center">AmpliPiper - Usage Tutorial</h1>
<br>
<div align="center">
    <img src="../imgs/AmpliPiper_logo.png">
</div>
<br>
<div style="display: flex;">
    <img src="../imgs/tettris.png" alt="TETTRIS project logo" style="width: 60%; height: auto; padding-right: 20px;">
    <img src="../imgs/nhm.svg.png" alt="NHM logo" style="width: 40%; height: auto; padding-left: 20px;">
</div>
<br>

## Table of contents

- [What is AmpliPiper - a short introduction](#what-is-amplipiper---a-short-introduction)
- [Installation - Let's get started](#installation---lets-get-started)
    + [Pre-requirements](#pre-requirements)
    + [Let's install](#lets-install-)
- [Usage - Prepare the data](#usage---prepare-the-data)
    + [Basecalling](#basecalling)
    + [Samples file and primers table](#samples-file-and-primers-table)
- [Usage - The pipeline at work](#usage---the-pipeline-at-work️)
- [Results - What to expect](#results---what-to-expect)
    + [The `data` subfolder](#data-subfolder)
    + [The `shell` subfolder](#shell-subfolder)
    + [The `results` subfolder](#results-subfolder)
    + [The `log` subfolder](#log-subfolder)
    + [Visualization](#visualization)
- [Results - How to correctly interpret and use them](#results---how-to-correctly-interpret-and-use-them)
    + [Phylogenetic trees](#phylogenetic-trees)
    + [Missing data](#missing-data)
    + [Species Delimitation](#species-delimitation)
    + [MSA and Genetic distances](#msa-and-genetic-distances)
    + []
- [Other useful resources](#resources)
    + [Install Linux on Windows as WSL](./install_wsl.md)
    + [Install Conda and Mamba](./install_conda_mamba.md)
    + [Install git](./install_git.md)
    + [Supplementary material](./suppl.md)

## What is AmpliPiper - a short introduction🤗

**AmpliPiper** is a comprehensive and modular pipeline that is able to perform a wide variety of downstream tasks on raw, basecalled, Oxford Nanopore long reads.

It is a pipeline developed by the Naturhistorisches Museum Wien to tackle most of the challenges posed by amplicon sequencing, and its ultimate goal is to provide the user with a fully proficient tool that is able to automatically perform the following tasks in a single run:

- Quality filtering and size selection
- Demultiplexing
- Consensus sequence generation
- Variant calling
- Haplotype calling
- Phylogenetic reconstruction
- Species delimitation
- Species identication

Browse [the documentation](https://astrabert.github.io/AmpliPiper-docs/) to learn more about AmpliPiper!

## Installation - Let's get started!🎬

### Pre-requirements

The pipeline was developed for POSIX-like operative systems (Linux distributions, e.g.), and tested on CentOS (AlmaLinux) and on Ubuntu. If you don't meet this basic requirement, you should consider installing a Linux subsystem on your machine, following the instructions [in this page](./install_wsl.md).

To install the pipeline, we will need [Mamba](https://mamba.readthedocs.io/en/latest/) and [Conda](https://docs.conda.io/en/latest/) up and running on our machine: these are two package installation and local environment managers, that will play a fundamental role in the modularization of the pipeline steps. If you don't have them, consider following the instructions to download the two softwares [in this page](./install_conda_mamba.md).

Once you are compliant with the aforementioned requirements, ensure that your shell has `git` installed, running:

```bash
git --help
```

If you don't have `git`, install it following [these instructions](./install_git.md).

### Let's install :)

The installation step are really simple and straightforward:

1. You first clone the GitHub repository:

```bash
git clone https://github.com/nhmvienna/AmpliPiper.git
```

2. Secondly, you go to the newly created `AmpliPiper` directory:

```bash
cd AmpliPiper
```

3. Lastly, you run the installation script!

```bash
bash shell/setup.sh
```

> _⚠️BEWARE!⚠️: The installation can take **really** long if you have an under-performing hardware (<16GB RAM, <4 core CPU), or if you have a slow connection_

Once the installation is completed, you will find an `envs` directory containing all the conda environments with the dependencies you need to run the pipeline: to go deeper, refer to [_Table 1_](./suppl.md#table-1).


## Usage - Prepare the data🍕

### Basecalling

When using AmpliPiper, you should employ raw basecalled Oxford Nanopore long reads in FASTQ format.

AmpliPiper **does not** come with a base-caller: we suggest, if you have an hardware good enough to support it (with GPU power), GUPPY (v.6.2.1 GPU) in super-high accuracy mode.

The resulting files should be demultiplexed with EXP-PBC096 kit, and compressed in gzip format.

### Samples file and primers table

With all your compressed files, you are now ready to start the real data preparation. 

> _For this example, we will be referring to small, artifically-generated FASTQ files contained in [tutorial_resources](./tutorial_resources/)_

AmpliPiper enables you to run analyses on multiple FASTQ files at a time, which should be listed in a file named `samples.csv` (naming is not strict, but it is used to make things clearer). In our case, it would be:

```
IND1,/path/to/AmpliPiper/tutorial/tutorial_resources/data/individual1.fastq.gz
IND2,/path/to/AmpliPiper/tutorial/tutorial_resources/data/individual2.fastq.gz
IND3,/path/to/AmpliPiper/tutorial/tutorial_resources/data/individual3.fastq.gz
IND4,/path/to/AmpliPiper/tutorial/tutorial_resources/data/individual4.fastq.gz
```
As you can see, the CSV file is **headerless** and contains two fields: the name of the sample and the path to the FASTQ file. It is not mandatory, but it is strongly advised to make use of absolute rather than relative path. 

AmpliPiper allows multi-locus analyses in one go, simply starting from your raw, non-per-locus-demultiplexed FASTQ files, as it performs demultiplexing: in order to do so, you should provide a primers table contained in a file named `primers.csv`. This table should not only encompass the IDs, forward and reverse sequence of the primers, but they should also contain the expected ploidy and length of the locus they refer to. This CSV file comes with a header, as you can see here:

```
ID,FWD,REV,PLOIDY,SIZE
COX1,AACTTTATACTTTCTCTTCGGAGCC,GAATAAATGTTGATAAAGAATTGGA,1,660
S7,ATGACCTCCCAGAATCACGGCTTTG,GCTTGCACTGCCAGCGTTCTCAGT,2,860
28S,CGGGCTTGGAAAAATTAGCGGGGAA,TGATAGGAAGAGCCGACATCGAAG,2,650
```

The naming of the header is **strict**.

And that's all the data preparation we need to start the analysis: let's now dive deep into the usage of the pipeline!

> _**PRO TIP!**: If you are editing the files from a non-Posix OS such as Windows, you may want to use the `dos2unix` command to make them Linux-accessible_

## Usage - The pipeline at work⚙️

The pipeline can be launched via command-line from your bash terminal, and it comes with many options that allow you to tweak and twist the results: you can find a detailed explanation of them in the [README of this repository](../README.md#run-the-pipeline) or in the [complete documentation](https://astrabert.github.io/AmpliPiper-docs/docs/index-test/#usage).

We will, for this tutorial, focus on the most important options and try to get a hang of them:

* `-k` or `--kthreshold`: Define the threshold _k_ for the maximum allowed proportion of mismatches for primer alignment during demultiplexing (default: 0.05). It should be set as quite strict (< 0.3) to ensure improved adherence of the demultiplexed sequences to the desired ones.
* `-r` or `--sizerange`: Define the allowed size buffer in basepairs around the expected locus length (default: 100). It should be kept fairly liberal ( 50 < r < 300): you should nevertheless be cautious when defining it, as you want to avoid including very short PCR byproducts or very long chimeras into the demultiplexed files
* `-n` or `--nreads`: Provide the absolute number or percentage of top-quality reads to consider for consensus sequence generation and variant calling (default: 500). 500 to 1000 reads is generally a sweet spot for the underlying consensus-generating software. However, consider exploring other possibilities if this does not produce the desired output.
* `-m` or `--minreads`: Set the minimum number of reads required for consensus sequence reconstruction (default: 100). A minimum threshold of less than 100 reads can lead to a lack of output from the underlying consensus generation software, while being too requiring (> 500 reads) may lead to lack of demultiplexed files.

Let's now build our command based on these parameters, and on the mandatory ones:

```bash
bash <path_to>/AmpliPiper/shell/AmpliPiper.sh \
    --samples <path_to>/AmpliPiper/tutorial/tutorial_resources/data/samples.csv \
    --primers <path_to>/AmpliPiper/tutorial/tutorial_resources/data/primers.csv \
    --output <path_to>/AmpliPiper/tutorial/tutorial_resources/results \
    --quality 2 \
    --nreads 500 \
    --sizerange 250 \
    --minreads 100 \
    --kthres 0.15 \
    --threads 10 \
    --force
```

As you can see in the command, in addition to the options we discussed beforehand you may also want to apply a quality filter (`-q` or `--quality`), set the number of threads on which to spam the computation workload of the pipeline (`-t` or `--threads`) and allow to forcibly rewrite the results folder if it already exists (`-f` or `--force`).

When you enter the command, you should start seeing some minimal logging, informative of the progressions of the pipeline.

If you want to monitor real-time and more in depth the various analysis step, you can always refer to the `log` folder which is created under the results directory shortly after the start of the pipeline execution, and that collects all the outputs from all the softwares executed throughout the pipeline.

## Results - What to expect🐛

After you run the pipeline, the output folder you chose at the beginning should look like this:

```
.
├── data
│   ├── demultiplexed
│   ├── filtered
│   └── raw
├── log
│   ├── SpecDelim
│   ├── SpecID
│   ├── ampliconsorter
│   ├── demulti
│   ├── html
│   ├── summary
│   ├── tree
│   └── variantcalling
├── results
│   ├── SpeciesDelim
│   ├── SpeciesID
│   ├── astraltree
│   ├── consensus_seqs
│   ├── haplotypes
│   ├── html
│   ├── summary
│   ├── tree
|   └── results.html
└── shell
    ├── demult1
    └── demult2
```

### `data` subfolder
This subfolder contains:

* **raw**: raw fastq files, which are a copy of the ones contained in the samples csv file.
* **filtered**: quality-filtered fastq files.
* **demultiplex**: demultiplexed files (with a subfolder for each sample, containing a file for each demultiplexed locus).

### `shell` subfolder
This subfolder contains:

* **demult1**: shell script files for quality filtering and demultiplexing (executed in parallel when allowed).
* **demult2**: shell script files for consensus generation and summarization (executed in parallel when allowed).

### `results` subfolder
This subfolder contains:

* **consensus_seqs**: in this sub-subfolder, we can find a directory for each sample, with a subdirectory for each locus: the Amplicon Sorter-generated consensus sequences are contained in this last level.
* **haplotypes**: in this sub-subfolder, we can find a directory for each locus. In this, we will have a fasta file with the concatenated consensus sequences across all samples for the specific locus, along with their alignments and the MSA visualization image.
* **tree**: in this sub-subfolder, we can find an ML tree for each locus, contained in a dedicated folder. You can also find the tree comparison heatmap and csv file.
* **astral**: in this sub-subfolder, we can find the Astral-generated tree for the concatenated loci, both in pdf and png format.
* **SpeciesDelim**: in this sub-subfolder, we can find the species delimitation outputs by ASAP for each locus and for Astral-concatenated loci. Here, there are partition-based representations, distance-based histograms and cumulative distributions, and text-based reports.
* **SpeciesID**: in this sub-subfolder, we can find the species identification outputs divided according to the database from which they were obtained and subdivided among the loci. For what concerns BOLD outputs, we'll have to pay attention to the _final.csv_ output, which summarizes the overall species identification results.
* **html**: in this sub-subfolder, we can find the html file to ease access to MSA visualization images.
* **results.html**: more about this file in the [Visualization](#visualization) section.

### `log` subfolder

This subfolder contains all the logs for the processes executed throughout the pipeline. It is advised that, if there is any problem with the pipeline, you trace it back here before reporting it, if possible.

### Visualization

The pipeline produces a `results.html` file to facilitate the visualization of the large amount of data output by the pipeline. Note that there are some compatibility issues with Chrome, so it is recommended to access it through other browsers, such as Firefox:

```bash
firefox results.html
```

## Results - How to correctly interpret and use them🔎

### Phylogenetic trees
The pipeline produces several maximum-likelihood phylogenetic trees, thanks to IQtree (single loci) and Astral (combined dataset). Nevertheless, all the resulting plots should be handled with care, as many of them may be the result of missing data, imperfect consensus reconstruction, chimeric pollution...

Anyone who would like to make use of them should first evaluate thoroughly their validity and test their solidity by comparing them with other approaches. 

These are the reasons why a watermark is placed on every tree image.

### Missing data
In the context of haplotypes reconstruction, the pipeline produces two outputs:

- A CSV summary containing information about size, ploidy and read depth for each reconstructed haplotype
- A heatmap indicating missing data across all the samples

This second output is really important, as it gives the user a complete landscape of available loci per sample for the analysis. It is higly reccomended to take a look to this plot, in order to have a clearer picture also of the phylogenetic trees output by the pipeline.


### Species delimitation

### MSA and genetic distances

When examining the pipeline's output in the case of Multiple Sequences Alignments, we can have a easier glimpse of the situation by looking at the visual representation of MSAs: if, when looking at them, you encounter densly gapped regions, you should:

1. Consider looking at the distance among the primers, in order to make an informed choice on how to set kthresh (might be too strict, might be too liberal)
2. View the genetic distances heatmap, as it might be informative in terms of remarkable outliers
3. Check the size and the size range you gave for the primers. 

Remember that some loci may be more sensitive to size range and k threshold than others, and that also the presence of chimeric or low-quality reads may affect the MSA output. In general, make sure to be working with a clean dataset.

## Resources

- [Install Linux on Windows as WSL](./install_wsl.md)
- [Install Conda and Mamba](./install_conda_mamba.md)
- [Install git](./install_git.md)
- [Supplementary material](./suppl.md)