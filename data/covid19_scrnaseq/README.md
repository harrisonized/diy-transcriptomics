Links to the data are available here:

1. [https://diytranscriptomics.com/project/lecture-13](https://diytranscriptomics.com/project/lecture-13)
2. [https://diytranscriptomics.com/project/lecture-14](https://diytranscriptomics.com/project/lecture-14)

Download the following files:

1. [pbmc\_1k\_v3\_scRNAseq_processed](https://drive.google.com/drive/folders/1RO45z5DEVpuaq5qwlF5QNdhc0tbGVK7l). This is a directory containing the following items:

	```
	pbmc_1k_v3_scRNAseq_processed/
	├─ counts_unfiltered/
	│   ├─ cellranger/
	│   │   ├─ barcodes.tsv
	│   │   ├─ genes.tsv
	│   │   └─ matrix.mtx
	│   ├─ cells_x_genes.barcodes.txt
	│   ├─ cells_x_genes.genes.txt
	│   └─ cells_x_genes.mtx
	├─ inspect.json
	└─ run_info.json
	```
2. [pbmc\_1k\_raw](https://drive.google.com/drive/folders/1DbLRO4kv-y3W06adFR26RdSaDPmfB4UA). This is a directory containing the following items:

	```
	pbmc_1k_raw/
	├─ pbmc_1k_v3_S1_mergedLanes_R1.fastq.gz
	└─ pbmc_1k_v3_S1_mergedLanes_R2.fastq.gz
	```
3. [t2g.txt](http://diytranscriptomics.github.io/Code/files/t2g.txt). I put this in the pbmc\_1k\_raw folder.

Your final directory should look like this:

```
covid19_scrnaseq/
├─ pbmc_1k_v3_scRNAseq_processed/
│   ├─ counts_unfiltered/
│   │   ├─ cellranger/
│   │   │   ├─ barcodes.tsv
│   │   │   ├─ genes.tsv
│   │   │   └─ matrix.mtx
│   │   ├─ cells_x_genes.barcodes.txt
│   │   ├─ cells_x_genes.genes.txt
│   │   └─ cells_x_genes.mtx
│   ├─ inspect.json
│   └─ run_info.json
├─ pbmc_1k_raw/
│   ├─ pbmc_1k_v3_S1_mergedLanes_R1.fastq.gz
│   ├─ pbmc_1k_v3_S1_mergedLanes_R2.fastq.gz
│   └─ t2g.txt
└─ pbmc_marker_list.json
```