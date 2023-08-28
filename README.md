## DIY Transcriptomics

Course link: [https://diytranscriptomics.com/](https://diytranscriptomics.com/)

This is where I followed along this RNA-seq data analysis course.


## Installation

1. Install Miniconda3 from  [here](https://docs.conda.io/en/latest/miniconda.html). I downloaded the `Miniconda3-latest-MacOSX-arm64.sh` and used `chmod +x` to make it exectuable. See: [https://protocols.hostmicrobe.org/conda](https://protocols.hostmicrobe.org/conda).

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
	export PATH="/Users/wanghc/homebrew/bin:$PATH"
	
	brew install kallisto  # Careful! See below.
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
	
	Check that it installed properly:
	
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
	If you see this error, you are on the right track, because the solution is documented in [this github issue](https://github.com/Irrationone/cellassign/issues/92#issuecomment-1154355997). Following the solution, fork the repo, then edit [inference-tensorflow.R line 65](https://github.com/Irrationone/cellassign/blob/master/R/inference-tensorflow.R#L165) to be this:
	
	```R
	p_y_on_c_norm <- tf$reshape(tf$reduce_logsumexp(p_y_on_c_unorm, 0L), as_tensor(shape(1,NULL)))
	```
	
	Finally, install tensorflow from your forked repo:
	
	```R
	devtools::install_github("your_fork/cellassign")
	```
	
	Alternatively, you could save yourself some hassle by installing from my fork:
	
	```
	devtools::install_github("harrison/cellassign")`
	```
	However, if you do this, just keep in mind that while it's generally okay to install from well-vetted git repos such as cellassign, installing repos from unknown developers such as myself is a security risk, because you never know if someone committed  malicious code that will execute when you try to run the program. Obviously, I didn't do that, but you get the point.
	
6. Installation errors I encountered:

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
	
7. Advanced troubleshooting:
	
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
| 4 | Lab 4 | R/leishmania_annotate\_gene\_expression.R | Use the ferret genome to annotate the leishmania dataset |
| 5 | Lab 6 | R/schistosoma_eda.R | Concatenates the abundance.tsv files in the schistosoma dataset and creates some violin plots || 6 | Lab 6 | R/schistosoma_pca.R | Performs PCA, then there are examples from gt, DT, and plotly. |
| 8 | Lab 7 | R/lemis_eda.R | EDA on the lemis dataset. |
| 9 | Lab 9 | R/malaria_eda.R | EDA on the malaria dataset. This is almost the same as the schistosoma_eda.R. |
| 10 | Lab 9 | R/malaria_pca.R | EDA on the malaria dataset. Some of this overlaps with schistosoma_pca.R, but this also includes a heatmap with dendrogram. |
| 11 | Lab 10 | R/covid19\_eda.R | EDA on the covid19 dataset. This is supposed to proceed similarly as the schistosoma and malaria datasets, but I'm going to keep this just as a placeholder. |
| 12 | Lecture 14 | R/covid19\_scrnaseq\_clustering | scRNAseq analysis on covid19_scrnaseq dataset. Does the clustering and heatmap generation. |
| 13 | Lecture 14 | R/covid19\_scrnaseq\_clustering | scRNAseq analysis on covid19_scrnaseq dataset. Does automatic cluster assignment.
| 14 | Lecture 14 | R/toxoplasma\_data\_integration | Compares two seurat objects from the toxoplasma dataset. |

The only script not explictly addressed is `Step7_functionalEnrichment.R` from Lecture 10. Also need to figure out how to install scater. All of lines requiring `scater::plotUMAP` have been hashed out.
