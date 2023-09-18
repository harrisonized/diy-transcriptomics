## DIY Transcriptomics

Course link: [https://diytranscriptomics.com/](https://diytranscriptomics.com/)

This is where I followed along Daniel Beiting's RNA-seq data analysis course. In my opinion, this was a good introductory course for getting started with using some of the basic libraries and functions to analyze some simple bioinformatics datasets, especially if you have no previous coding experience. I liked that all the scripts were designed to be able to run locally, which makes it easy to really examine the objects in RStudio and figure out how the code works.

However, during the course, I personally found it difficult to follow along, because the provided scripts contained a lot of antipatterns, which obscured the logic of scripts that should have otherwise been straightforward. While this may be acceptable for beginner users learning how to code for the first time, it makes the code far less modifiable or reusable. I found myself having to really work at getting each line of code to work. Having to focus most of my time and effort on fixing broken lines of code and fixing library installations took time away from the most important part of the course: learning what the code actually does, learning which algorithms to use, and learning how to improve the visualizations.

After the class ended, I still wanted to make sure I made it through every aspect of the class. Therefore, in this repo, I meticulously rewrote each script in such a way that it runs in the terminal directly out of the box. Major sections are demarcated by log messages so that there is feedback as the code runs, and these major sections are standardized so you know exactly what to expect:

1. imports at the top
2. command line options after imports
3. import data
4. data wrangling
5. output figures

That way, if you have a custom analysis, you should be able to reuse the same script, swap out the data, and make minimal adjustments to do roughly the same analysis. Furthermore, you will be able to copy code chunks to quickly put together custom scripts, much like you can assemble a castle out of legos. 

One caveat is that if you are currently taking DIYTranscriptomics and the assignments are still the same, keep in mind that this repo is NOT meant to be a cheat sheet, even if it generates the figures that are supposed to be "the solution." This is only meant to make your life a little bit easier, not solve all your problems for you. You will still have to understand the data you're working with. You will still have to manually figure out which subset of the data you're interested in. You will still have to do work to make your figures look pretty. I left this repo at the bare minimum of the solution, because if I were to invest time making things pretty, I'd rather do it for my own datasets and projects (sorry!).

If you try this code and you run into problems, please reach out to me through my Penn email, and I will be happy to answer your questions.

## Library Installations

This is by far the most frustrating part of the course, so I made sure to document every step in installing it in my 2021 Macbook Pro (M1 chip). Sorry to the Windows folks. Hope you guys get a good TA that uses Windows!

1. Install Miniconda3 from [here](https://docs.conda.io/en/latest/miniconda.html). I downloaded the `Miniconda3-latest-MacOSX-arm64.sh` and used `chmod +x` to make it exectuable. See: [https://protocols.hostmicrobe.org/conda](https://protocols.hostmicrobe.org/conda).

	In general, because trying to resolve issues with the base environment are challenging and can potentially destroy your entire conda installation, it is best practice to never use the base environment for anything and to only ever work in isolated environments. Therefore, create the `rnaseq` conda environment.
	
	```bash
	conda create --name rnaseq
	conda activate rnaseq
	conda config --add channels bioconda
	conda install pip
	
	# packages
	pip install multiqc
	conda install fastqc
	pip install sourmash
	```
	
2. Install Kallisto. On the Apple M1 Pro, you will need homebrew.

	```
	# install brew and add to $PATH
	git clone https://github.com/Homebrew/brew homebrew
	eval "$(homebrew/bin/brew shellenv)"
	brew update --force --quiet
	export PATH="~/homebrew/bin:$PATH"  # you can replace the ~ with your path to home if this doesn't work
	
	brew install kallisto
	```
	
	Note that this installation is completely separate from conda. If you use use `pip install kallisto`, you will instead get a package from AstraZeneca with the same name, which is used to calculate atomic features. This is not what we want.
	
3. Install Centrifuge. On the Apple M1 Pro, you will need Rosetta 2 to make this work.

	```
	# install Rosetta 2
	/usr/sbin/softwareupdate --install-rosetta --agree-to-license
	
	# Download Centrifuge from here: http://www.ccb.jhu.edu/software/centrifuge/index.shtml
	# Releases > Mac OS X x86_64 binary
	# Add this to .zshrc:
	export PATH="$HOME/centrifuge-1.0.3-beta:$PATH
	```
	
4. Install the following packages in R:

	```R
	# regular packages
	install.packages('hexbin')
	install.packages('cowplot')
	install.packages('gt')
	install.packages("webshot2")
	install.packages('DT')
	install.packages("tidyverse")
	install.packages("plotly")
	install.packages('textshape')
	install.packages('R2HTML')
	
	# biology-focused packages
	BiocManager::install("rhdf5")
	BiocManager::install("tximport")
	BiocManager::install("ensembldb")
	BiocManager::install("beepr")
	BiocManager::install("EnsDb.Hsapiens.v86")
	BiocManager::install("BSgenome.Mfuro.UCSC.musFur1")  # Mustela putorius furo
	# BiocManager::install("biomaRt")  # already installed
	BiocManager::install("edgeR")
	BiocManager::install("DropletUtils")
	BiocManager::install("scran")
	BiocManager::install("SingleR")
	BiocManager::install('celldex')
	
	# difficult-to-install packages
	# try this later
	# BiocManager::install("densvis")  # ERROR: compilation failed for package ‘densvis’
	# BiocManager::install("scater")  # ERROR: dependency ‘densvis’ is not available for package ‘scater’
	```

5. Install Tensorflow and CellAssign.

	Important! Read through this entire section before attempting the installation. This was by far one of the most challenging installations I've ever encountered.
	
	In R, the `reticulate` package is used to allow R to access python through conda. By default, it will look for python at `/usr/bin/python3` or `~/miniconda3/bin/python` (your base environment), and then it will be extremely confusing when you run `library('cellassign')` and you get the following error, even though you swear you installed it:
	
	```
	Error: package or namespace load failed for ‘cellassign’:
	  .onLoad failed in loadNamespace() for 'cellassign', details:
	  call: fun(libname, pkgname)
	  error: Tensorflow installation not detected. Please run 'tensorflow::install_tensorflow()' to continue...
	```
	
	Because it is best practice not to work in the base environment, I recommend creating the `r-reticulate` environment. The name comes from the default environment R will create for you if cellassign can't find Tensorflow. Run the following command in your terminal:
	
	```bash
	conda create --name r-reticulate python=3.9
	```
	
	Next, install tensorflow using RStudio. This is because if you instead opt to install tensorflow directly using conda in the terminal, it is possible that at the end of it, R may be unable to detect your installation.
	
	```R
	install.packages("tensorflow")
	reticulate::use_condaenv('r-reticulate')  # IMPORTANT!
	tensorflow::install_tensorflow(extra_packages='tensorflow-probability')
	```
	Also note that unfortunately, you cannot specify which version of tensorflow you'll get. If you try, you will get the following warning:
	
	```
	# Only `version = 'default'` supported on Arm Macs at this time.
	Request for Tensorflow version '2.1.0' ignored.
	```
	When installing tensorflow through RStudio, you will be left with some inconsistencies, which will need to be cleaned up. In your terminal, run the following:
	
	```bash
	conda activate r-reticulate
	
	# this causes an annoying import error
	pip uninstall tensorflow-metal
	
	# fix the inconsistencies
	conda install -c anaconda tensorflow-base
	conda install -c conda-forge tensorflow-probability=0.18.0
	
	# reinstall keras
	pip uninstall keras
	pip install keras==2.13.1
	```
	
	On my system, this gave me tensorflow                2.13.0, tensorflow-base=2.10.0, tensorflow-deps=2.9.0,tensorflow-probability=0.18.0, and keras=2.13.1. Do NOT install tensorflow using pip via `pip install tensorflow`. This will result in an import error.
	
	Check that it installed properly. In your R console:
	
	```
	reticulate::use_condaenv('r-reticulate')
	reticulate::import("tensorflow")
	reticulate::import("tensorflow_probability")
	reticulate::py_discover_config("keras")
	```
	
	When the imports work, you'll get a nice return message, eg. `Module(tensorflow)`. If Keras isn't installed properly, when you run `reticulate::py_discover_config("keras")`, you'll see `Keras: [NOT FOUND]`. If this happens to you, check that there isn't an extra pointer in your conda environment. When you run `conda list` and `pip list`, Keras should only show up on `pip list`. If it also shows up on `conda list`, especially if there's a different version, navigate to `~/miniconda3/envs/r-reticulate/lib/python3.9/site-packages` and remove the extra one. For example, for me, it was `rm -rf keras-2.9.0.dist-info`. This happens because you use both conda and pip to install it, and pip may not remove the conda pointer when overwriting the folder.
	
	Finally, after all that, you are ready to install `cellassign`. But not so fast! You cannot install it directly from `Irrationone/cellassign`, because doing so will cause you to run into the following error when you go to run cellassign:
	
	```
	ValueError: Tried to convert 'shape' to a tensor and failed.
	Error: Cannot convert a partially known TensorShape (1, ?) to a Tensor.
	```
	If you see this error, you are on the right track, because the solution is documented in [this github issue](https://github.com/Irrationone/cellassign/issues/92#issuecomment-1154355997). Following the solution, what you should do instead is fork the repo, then edit [inference-tensorflow.R line 65](https://github.com/Irrationone/cellassign/blob/master/R/inference-tensorflow.R#L165) to be this:
	
	```R
	p_y_on_c_norm <- tf$reshape(tf$reduce_logsumexp(p_y_on_c_unorm, 0L), as_tensor(shape(1,NULL)))
	```
	
	Finally, install tensorflow from your forked repo:
	
	```R
	devtools::install_github("your_fork/cellassign")
	```
	
	Alternatively, you could save yourself some hassle by installing from my fork:
	
	```
	devtools::install_github("harrisonized/cellassign")`
	```
	However, if you do this, just keep in mind that while it's generally okay to install from well-vetted git repos such as cellassign, installing repos from unknown developers such as myself is a security risk, because you never know if someone committed  malicious code that will execute when you try to run the program. Obviously, I didn't do that, but you get the point.
	
6. The only library I was unable to install is scater. As such, I hashed out all the lines requiring `scater::plotUMAP`. This is a work-in-progress, I'll get around to it at some point.
	
7. Tensorflow installation errors I encountered:

	If your tensorflow-probability version is too low (eg. 0.14.0), you will encounter the following error message:
	
	```
	"ImportError: cannot import name 'deserialize_keras_object' from partially initialized module 'keras."
	```
	
	If your tensorflow-probability version is too high (eg. 0.20.0), you will encounter the following error message:
	
	```
	Error in py_module_import(module, convert = convert) : 
	AttributeError: module 'tensorflow.python.framework.type_spec' has no attribute '_NAME_TO_TYPE_SPEC'
	```
	
	If keras is not installed or detected, you will see the following error:
	
	```
	Error in py_module_import(module, convert = convert) : 
	AttributeError: module 'keras.api._v2.keras' has no attribute 'layers'
	```

    If your keras version is inconsistent with one of your other packages, then you will see the following error:
	
	```
	Error in py_module_import(module, convert = convert):
	AttributeError: module 'tensorflow.python.data.ops.from_tensor_slices_op' has no attribute '_TensorSliceDataset'
	```
	
	If you have the right version of keras (>=2.10.0), then you will see the following error:
	
	```
	ValueError: Tried to convert 'shape' to a tensor and failed.
	Error: Cannot convert a partially known TensorShape (1, ?) to a Tensor.
	```
	
	If you have a low version of keras (<=2.4.3), you will see the following error:
	
	```
	ImportError: Keras requires TensorFlow 2.2 or higher. Install TensorFlow via `pip install tensorflow`
	```
	
	If you start to install tensorflow without first running `reticulate::use_condaenv('r-reticulate')`, tensorflow will be installed in your base environment, which is NOT what you want. The reason is that you can potentially run into the following issue, and then you'll have to create a new environment anyway.
	
	```
	Solving environment: failed with initial frozen solve. Retrying with flexible solve.
	    Solving environment: / 
	    Found conflicts!
	```
	
	If you are an advanced user, you can always try to reset your base environment by running the following command: `conda install --revision=0`, but I can't guarantee this will resolve the issue.
	
	Suppose you started to install tensorflow and then you realized that it was installing into the wrong environment. Unfortunately, you cannot switch halfway in between. R will give you the following error:
	
	```
	ERROR: The requested version of Python ('~/miniconda3/envs/r-reticulate/bin/python') cannot be used,
	as another version of Python ('~/miniconda/bin/python') has already been initialized.
	```	
	
	If you get this, simply restart R, and make sure you point reticulate to the correct environment upon startup.
	
	When I was trying to resolve conflicts while installing Tensorflow on my Apple M1 Pro, I messed up in a dramatic way and accidentally destroyed my base environment! Specifically, I ran `pip uninstall -y -r <(pip freeze)` after reading [this StackOverflow post](https://stackoverflow.com/questions/41914139/how-to-reset-anaconda-root-environment), which left my conda base environment in an inconsistent, irrcoverable state. Afterward, I ended up having to do a clean install. When you do a clean install, it is possible to save your existing environments by moving the `miniconda3/envs` folder somewhere else. Then, after reinstalling conda, move that folder back. (Unfortunately, I forgot to do this too.) Please make sure you do not repeat these mistakes!


## Scripts

Listed in chronological order:

| Script Number | Source | Script Name | Description |
| :--- | ------ | ------ | ----------- |
| 1 | Lecture 2 | bash/build\_kallisto\_index.sh  | Use Kallisto to build an index on leishmania dataset |
| 2 | Lab 2 | bash/map\_reads.sh | Use Kallisto to map reads to the index file on leishmania dataset |
| 3 | Lab 5 | bash/estimate\_sequence\_similarity.sh | Use Sourmash to identify non-human reads. Use Centrifuge to search for gene signatures. |
| 4 | Lab 4 | R/query\_biomart.R | This standalone script provides an example for how to query bioMart's database. |
| 5 | Lab 4 | R/leishmania_eda.R | Starter script for EDA on the covid19 dataset. This script imports the kalliso outputs from the leishmania dataset, using 'EnsDb.Hsapiens.v86' to annotate each row of the resulting dataframe. |
| 6 | Lab 6 | R/schistosoma_eda.R | Concatenates the abundance.tsv files in the schistosoma dataset and creates some violin plots || 7 | Lab 6 | R/schistosoma_pca.R | Performs PCA, then there are examples from gt, DT, and plotly. |
| 8 | Lab 7 | R/lemis_eda.R | EDA on the lemis dataset. |
| 9 | Lab 9 | R/malaria_eda.R | EDA on the malaria dataset. This is almost the same as the schistosoma_eda.R. |
| 10 | Lab 9 | R/malaria_pca.R | Performs PCA on the malaria dataset. Some of this overlaps with schistosoma_pca.R, but this also includes a heatmap with dendrogram. |
| 11 | Lab 9 | R/malaria\_deg\_analysis.R | Performs differential gene expression analysis on the malaria dataset. |
| 12 | Lab 10 | R/covid19\_eda.R | Starter script for EDA on the covid19 dataset. |
| 13 | Lecture 14 | R/covid19\_scrnaseq\_qc.R | scRNAseq analysis on covid19_scrnaseq dataset. Filters empty drops and creates some QC plots. |
| 14 | Lecture 14 | R/covid19\_scrnaseq\_seurat.R | scRNAseq analysis on covid19_scrnaseq dataset. Creates a Seurat object from the empty-drop-filtered cellranger output. Plots UMAPs and heatmaps of DEGs. |
| 15 | Lecture 14 | R/covid19\_scrnaseq\_assign\_clusters.R | scRNAseq analysis on covid19_scrnaseq dataset after filtering. Autoassigns cluster identities based on publically available datasets. |
| 16 | Lecture 14 | R/toxoplasma\_data\_integration.R | Compares two seurat objects from the toxoplasma dataset. |
| 17 | Lecture 10 | R/malaria\_functional\_enrichment.R | This is the only script that has not been addressed, because it was skipped during the main course. |

In addition, there are two utility files (`scrnaseq_qc_plots.R` and `utils.R`) used for reusable functions used in some of the scripts.