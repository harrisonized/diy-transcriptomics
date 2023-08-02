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
install.packages('DT')
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

| Step | Script | Description | Source |
| :--- | ------ | ----------- | ------ |
| 1 | bash/build\_kallisto\_index.sh | Use Kallisto to build an index | Lecture 2 |
| 2 | bash/map\_reads.sh | Use Kallisto to map reads to the index file | Lab 2 |
| 3 | bash/estimate\_sequence\_similarity.sh | Use Sourmash to identify non-human reads. Use Centrifuge to search for gene signatures. | Lab 5|
| 4 | `Rscript R/lab_4-annotate_gene_expression.R` || Lab 4|
| 5 | `Rscript R/Step1_TxImport.R` || Lab 6 || 6 | `Rscript R/Step2_dataWrangling.R` || Lab 6 || 7 | `Rscript R/Step3_multivariate.R` || Lab 6 |
| 8 | `Rscript R/lab_7.R  # standalone` || Lab 7 |
| 9 | `Rscript R/lab_9-malaria.R` || Lab 9 |

