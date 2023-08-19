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
BiocManager::install("DropletUtils")

install.packages('hexbin')
install.packages('cowplot')
install.packages('gt')
install.packages("webshot2")
install.packages('DT')
install.packages("tidyverse")
install.packages("plotly")
install.packages('textshape')
install.packages('R2HTML')
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

Listed in chronological order:

| Item | Source | Script Name | Description |
| :--- | ------ | ------ | ----------- |
| 1 | Lecture 2 | bash/build\_kallisto\_index.sh  | Use Kallisto to build an index on leishmania dataset |
| 2 | Lab 2 | bash/map\_reads.sh | Use Kallisto to map reads to the index file on leishmania dataset |
| 3 | Lab 5 | bash/estimate\_sequence\_similarity.sh | Use Sourmash to identify non-human reads. Use Centrifuge to search for gene signatures. |
| 4 | Lab 4 | R/leishmania_annotate\_gene\_expression.R | Use the ferret genome to annotate the leishmania dataset |
| 5 | Lab 6 | R/schistosoma_eda.R | Concatenates the abundance.tsv files in the schistosoma dataset and creates some violin plots || 6 | Lab 6 | R/schistosoma_pca.R | Performs PCA, then there are examples from gt, DT, and plotly. |
| 8 | Lab 7 | R/lemis_eda.R | EDA on the lemis dataset. |
| 9 | Lab 9 | R/malaria_eda.R | EDA on the malaria dataset. This is almost the same as the schistosoma_eda.R. |
| 10 | Lab 9 | R/malaria_pca.R | EDA on the malaria dataset. Some of this overlaps with schistosoma_pca.R, but this also includes a heatmap with dendrogram. |
| 11 | Lab 10 | R/covid19_eda.R | EDA on the covid19 dataset. Placeholder for now. |
| 12 | Lab 13 | R/covid19_scrnaseq.R | scRNAseq analysis on covid19_scrnaseq dataset. Placeholder for now. |
