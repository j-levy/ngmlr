# About assessing ngmlr results

You can build, run and assess ngmlr results in a relatively painless way.

- install miniconda3
- clone this and create the conda environment from `conda.yml`
- build ngmlr
- find some DNA you want to test with
- run the thing

### 1. install miniconda3

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda3
rm miniconda.sh
```

You'd probably want to add the registering hook to your ~/.bashrc, ~/.zshrc or whatever

```bash
$HOME/miniconda3/bin/conda init $SHELL
```

### 2. clone this ans create conda env

This will let you download all the tools you need to simulate reads, simulate SVs, etc. Lots of stuff. Plus Python.

```bash
git clone git@github.com:j-levy/ngmlr
cd ngmlr
conda env create -f conda.yml
conda activate ngmlr
```

### 3. build ngmlr

You'll need usual build tools installed.

```bash
mkdir build
cd build
cmake .. && make -j4
```

### 4. Find some DNA you'd like to test with

For example c_elegans DNA (it's a worm)

```bash
NCBIURL=ftp://ftp.wormbase.org/pub/wormbase/releases/WS245/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS245.genomic.fa.gz
NCBIFILENAME=$(basename $NCBIURL | sed "s/.gz$//")
OUTDIR=c_elegans
wget --show-progress -qO- "$NCBIURL" | zcat >> "$OUTDIR/$NCBIFILENAME"
ln -s *.fna c_elegans.fasta
```

### Run the thing

edit default path in `scripts/subsegments_run.py` to put the one from your DNA.

You can look at the functions in `align.py` and build your own run script. `subsegments_run` is more like a first attempt. 

```bash
conda activate ngmlr
cd scripts
python -m subsegments_run
```

Sit back and wait for the script to poop out the result.
