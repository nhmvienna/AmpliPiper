## Install all required software

## function to assess realpath to setup.sh
realpath_custom() {
    local path="$1"

    # Check if the path is absolute or relative
    if [[ "$path" != /* ]]; then
        path="$PWD/$path"
    fi

    # Normalize the path by resolving '..' and '.'
    while [[ "$path" =~ (/[^/]+/\.{2}|/\./|//) ]]; do
        path="${path//\/.\//\/}"             # Replace '/./' with '/'
        path="${path//\/\/+/\/}"             # Replace '//' with '/'
        path="${path/\/[^\/]+\/\.\.\/\//\/}" # Resolve '/../'
    done

    # Resolve symlinks
    while [ -L "$path" ]; do
        link=$(readlink "$path")
        if [[ "$link" == /* ]]; then
            path="$link"
        else
            path="${path%/*}/$link"
        fi
    done

    # Output the final resolved path
    echo "$path"
}

## change to home directory of the shell folder which contains the setup script
tmp=$(realpath_custom $(dirname $0))
BASEDIR=${tmp%/*}
echo "Installing necessary software"
echo "*********************"
echo ""
echo "testing if Conda and Mamba are installed..."

command -v conda >/dev/null 2>&1 || {
    echo >&2 "The installation pipeline requires Anaconda/Miniconda but it is not installed. Please check here: https://anaconda.org/ for more details. Aborting."
    echo "Conda not installed" >>${BASEDIR}/logs/dep.err
    exit 1
}

command -v mamba >/dev/null 2>&1 || {
    echo >&2 "The installation pipeline requires Mamba but it is not installed. Please check here: https://github.com/conda-forge/miniforge#mambaforge for more details. Aborting."
    echo "Mamba not installed" >>${BASEDIR}/logs/dep.err
    exit 2
}

printf "\n...Success!\n"
echo "*********************"
echo ""
echo "Have a cup of coffee, this may take a while... "
echo '''
   ( (
    ) )
  ........
  |      |]
  \      /
   `----´ '''

echo ""
sleep 2
os="$(uname -s)"

if [ -d ${BASEDIR}/envs ]; then
    rm -rf ${BASEDIR}/envs
fi

if [ -d ${BASEDIR}/envs/logs ]; then
    rm -rf ${BASEDIR}/envs/logs
fi

mkdir -p ${BASEDIR}/envs/logs
eval "$(conda shell.bash hook)"

## install python_dependencies
mamba create \
    -p ${BASEDIR}/envs/python_dependencies \
    -y \
    -c conda-forge \
    -c bioconda \
    python=3.10 \
    sed >>${BASEDIR}/envs/logs/setup.log 2>&1

# activate python_dependencies environments

conda activate \
    ${BASEDIR}/envs/python_dependencies

mamba install \
    -y \
    -c conda-forge \
    -c bioconda \
    matplotlib biopython edlib cairosvg pymsaviz ete3 beautifulsoup4 scipy scikit-learn pandas >>${BASEDIR}/envs/logs/setup.log 2>&1

${BASEDIR}/envs/python_dependencies/bin/pip3 install --no-cache-dir numpy pandas >>${BASEDIR}/envs/logs/setup.log 2>&1
echo "Python dependencies installation done"

## install AmpliconSorter
# clone python_dependencies GitHub repo
cd ${BASEDIR}/envs/python_dependencies
git clone https://github.com/avierstr/amplicon_sorter >>${BASEDIR}/envs/logs/setup.log 2>&1

# make exectuable and change the interpreter to the python installation in the environment
echo \#\!${BASEDIR}/envs/python_dependencies/bin/python3.10 \
    >bin/amplicon_sorter.py
awk 'NR>1' ${BASEDIR}/envs/python_dependencies/amplicon_sorter/amplicon_sorter.py >>bin/amplicon_sorter.py
chmod 755 bin/amplicon_sorter.py

if [ ! -d ${BASEDIR}/envs/python_dependencies/amplicon_sorter/ ]; then
    echo "Ampliconsorter installation failed" >>${BASEDIR}/envs/logs/setup.err
fi
echo "Ampliconsorter installation done"

conda deactivate

## install ASAP from local package
conda create \
    -p ${BASEDIR}/envs/asap \
    ${BASEDIR}/.resources/asap-linux-v0.1.2-0.tar.bz2 >>${BASEDIR}/envs/logs/setup.log 2>&1

if [ ! -d ${BASEDIR}/envs/asap ]; then
    echo "asap installation failed" >>${BASEDIR}/envs/logs/setup.err
fi
echo "asap installation done"

## install ASTER
mamba create \
    -p ${BASEDIR}/envs/aster \
    -y \
    -c bioconda \
    -c conda-forge \
    aster >>${BASEDIR}/envs/logs/setup.log 2>&1

if [ ! -d ${BASEDIR}/envs/aster ]; then
    echo "ASTER installation failed" >>${BASEDIR}/envs/logs/setup.err
fi
echo "ASTER installation done"

## install GNU parallel
mamba create \
    -p ${BASEDIR}/envs/parallel \
    -y \
    -c conda-forge \
    parallel >>${BASEDIR}/envs/logs/setup.log 2>&1

if [ ! -d ${BASEDIR}/envs/parallel ]; then
    echo "GNU parallel installation failed" >>${BASEDIR}/envs/logs/setup.err
fi
echo "GNU parallel installation done"

## install IQtree
mamba create \
    -p ${BASEDIR}/envs/iqtree \
    -y \
    -c bioconda \
    -c conda-forge \
    iqtree=2.3.6 >>${BASEDIR}/envs/logs/setup.log 2>&1

if [ ! -d ${BASEDIR}/envs/iqtree ]; then
    echo "IQtree installation failed" >>${BASEDIR}/envs/logs/setup.err
fi
echo "IQtree installation done"

## install MAFFT
mamba create \
    -p ${BASEDIR}/envs/mafft \
    -y \
    -c bioconda \
    -c conda-forge \
    mafft=7.525 >>${BASEDIR}/envs/logs/setup.log 2>&1

if [ ! -d ${BASEDIR}/envs/mafft ]; then
    echo "MAFFT installation failed" >>${BASEDIR}/envs/logs/setup.err
fi
echo "MAFFT installation done"

## install Chopper
mamba create \
    -p ${BASEDIR}/envs/chopper \
    -y \
    -c bioconda \
    -c conda-forge \
    chopper >>${BASEDIR}/envs/logs/setup.log 2>&1

if [ ! -d ${BASEDIR}/envs/chopper ]; then
    echo "chopper installation failed" >>${BASEDIR}/envs/logs/setup.err
fi

conda activate \
    ${BASEDIR}/envs/chopper

mamba install \
    -y \
    -c conda-forge \
    -c bioconda \
    pigz >>${BASEDIR}/envs/logs/setup.log 2>&1

conda deactivate
echo "chopper installation done"

## install Pigz
mamba create \
    -p ${BASEDIR}/envs/pigz \
    -y \
    -c bioconda \
    -c conda-forge \
    pigz >>${BASEDIR}/envs/logs/setup.log 2>&1

if [ ! -d ${BASEDIR}/envs/pigz ]; then
    echo "PIGZ installation failed" >>${BASEDIR}/envs/logs/setup.err
fi
echo "PIGZ installation done"

## install R
mamba create \
    -p ${BASEDIR}/envs/R \
    -y \
    -c bioconda \
    -c conda-forge \
    r-base=4.1.3 >>${BASEDIR}/envs/logs/setup.log 2>&1

if [ ! -d ${BASEDIR}/envs/R ]; then
    echo "R installation failed" >>${BASEDIR}/envs/logs/setup.err
fi

conda activate \
    ${BASEDIR}/envs/R

# install dependencies
mamba install \
    -y \
    -c conda-forge \
    -c bioconda \
    -c r \
    bioconductor-ggtree \
    r-tidyverse \
    bioconductor-treeio \
    r-phytools \
    r-ape \
    r-treedist \
    r-reshape2 \
    r-seqinr \
    r-gridextra \
    r-stringr \
    >>${BASEDIR}/envs/logs/setup.log 2>&1

conda deactivate

echo "R installation done"
