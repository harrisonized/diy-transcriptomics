#!/usr/bin/env bash
set -e  # exit on error


# params
export index_filepath='data/lab_2/ref/fasta/Homo_sapiens.GRCh38.cdna.all.index'
export input_dir="data/lab_2/fastq"
export output_dir='data/lab_2/mapped_reads'
export log_dir='logs'
export troubleshooting=true

# set working directory
export wd="$HOME/github/diy-transcriptomics"
cd $wd
pwd
echo ""

# set conda
export CONDA=$HOME/anaconda3/etc/profile.d/conda.sh
export CONDA_ENV=rnaseq
source $CONDA
conda activate $CONDA_ENV


# make output_dir if not exists
if ! test -d "$output_dir" ; then
    mkdir $output_dir
    echo "TRUE"
fi


# for loop
echo "start time: $(date +'%Y%m%d_%H%M%S')"
input_files_glob="$input_dir/*"
declare -a input_files=($input_files_glob)
for input_file in ${input_files[@]}
do
    # get input_filename
    declare -a input_path_levels=($(echo $input_file | tr "/" " "))  # input_file.split('/')
    input_filename=${input_path_levels[${#input_path_levels[@]}-1]}  # get last element
    
    filename_no_ext=${input_filename/.fastq.gz/""}
    output_subdir="$output_dir/$filename_no_ext"
    log_filepath="$filename_no_ext.log"
    current_time=$(date +'%Y%m%d_%H%M%S')
    
    # troubleshooting
    if $troubleshooting; then
        echo "current_time: $current_time"
        echo "index_filepath: $index_filepath"
        echo "output_subdir: $output_subdir"
        echo "input_file: $input_file"
        echo "log_filepath: logs/$log_filepath"
        echo ""
    else
        echo "Processing $input_file ... $current_time"
        kallisto quant --index=$index_filepath --output-dir=$output_subdir \
        	-t 8 --single -l 250 -s 30 $input_file &> "$log_dir/$log_filepath"
    fi

done

current_time=$(date +'%Y%m%d_%H%M%S')
echo "END $current_time"
