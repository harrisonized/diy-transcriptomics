## DIY Transcriptomics

Course link: [https://diytranscriptomics.com/](https://diytranscriptomics.com/)

This is where I followed along this RNA-seq data analysis course.


## Installation

Install the following R packages:

```R
BiocManager::install("rhdf5")
BiocManager::install("tximport")
BiocManager::install("ensembldb")
BiocManager::install("beepr")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("BSgenome.Mfuro.UCSC.musFur1")  # Mustela putorius furo
# BiocManager::install("biomaRt")  # already installed
BiocManager::install("edgeR")
install.packages('hexbin')
install.packages('cowplot')
install.packages('gt')
install.packages("webshot2")
install.packages('DT')
install.packages("tidyverse")
install.packages("plotly")
```

Installing the RNAseq conda environment:
See: [https://protocols.hostmicrobe.org/conda](https://protocols.hostmicrobe.org/conda)

```bash
conda create --name rnaseq
conda activate rnaseq
conda config --add channels bioconda
conda install pip

# install brew and add to $PATH
git clone https://github.com/Homebrew/brew homebrew
eval "$(homebrew/bin/brew shellenv)"
brew update --force --quiet
export PATH="/Users/wanghc/homebrew/bin:$PATH"

# packages
pip install multiqc
brew install kallisto
conda install fastqc
pip install sourmash

# install Rosetta 2
/usr/sbin/softwareupdate --install-rosetta --agree-to-license

# Download Centrifuge from here: http://www.ccb.jhu.edu/software/centrifuge/index.shtml
# Releases > Mac OS X x86_64 binary
# Add this to .zshrc:
export PATH="$HOME/centrifuge-1.0.3-beta:$PATH


```

Note: Kallisto is also the name of a package from AstraZeneca used to calculate atomic features, which is what you will get if you run `pip install kallisto`, and it is the wrong package to install.


## Scripts

To be run in order:

| Step | Script | Source | Description |
| :--- | ------ | ------ | ----------- |
| 1 | bash/build\_kallisto\_index.sh | Lecture 2 | Use Kallisto to build an index on leishmania dataset |
| 2 | bash/map\_reads.sh | Lab 2 | Use Kallisto to map reads to the index file on leishmania dataset |
| 3 | bash/estimate\_sequence\_similarity.sh | Lab 5 | Use Sourmash to identify non-human reads. Use Centrifuge to search for gene signatures. |
| 4 | R/leishmania_annotate\_gene\_expression.R | Lab 4 | Use the ferret genome to annotate the leishmania dataset |
| 5 | schistosoma_eda.R | Lab 6 | Concatenates the abundance.tsv files in the schistosoma dataset and creates some violin plots || 6 | schistosoma_pca.R | Lab 6 | Performs PCA, then there are examples from gt, DT, and plotly. |
| 8 | lemis_eda.R | Lab 7 | EDA on the lemis dataset. |
| 9 | Rscript R/lab_9-malaria.R | Lab 9 ||

