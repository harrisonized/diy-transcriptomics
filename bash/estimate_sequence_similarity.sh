#!/usr/bin/env bash
set -e  # exit on error

# set working directory
export wd="$HOME/github/diy-transcriptomics"
cd $wd
pwd

# conda
export CONDA=$HOME/miniconda3/etc/profile.d/conda.sh
source $CONDA
conda activate 'rnaseq'

# this follows lab 5: https://diytranscriptomics.com/lab/lab-05

# ----------------------------------------------------------------------
# Task 1: Use fastqc to get a summary

export DATA_DIR="data/leishmania_revisited"
export FASTQ_FILENAME="SRR8668774_dehosted.fastq.gz"
export FASTQ_FILE="$DATA_DIR/$FASTQ_FILENAME"
export GENOME_FILE="$DATA_DIR/reference_sequences/genbank-k31.lca.json.gz"

fastqc data/leishmania_revisited/SRR8668774_dehosted.fastq.gz -t 8

# ----------------------------------------------------------------------
# Task 2

# time = ~2min to sketch ~9M reads
sourmash sketch dna -p scaled=10000,k=31,abund $FASTQ_FILE --name-from-first
mv "$FASTQ_FILENAME.sig" "data/leishmania_revisited/$FASTQ_FILENAME.sig"

# time = ~2min
# no output?
sourmash gather -k 31 "data/leishmania_revisited/$FASTQ_FILENAME.sig" $GENOME_FILE
# == This is sourmash version 4.6.1. ==
# == Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

# selecting specified query k=31
# loaded query: SRR8668774.33 33/1... (k=31, DNA)
# --
# loaded 93249 total signatures from 1 locations.
# after selecting signatures compatible with search, 93249 remain.

# Starting prefetch sweep across databases.
# Found 17 signatures via prefetch; now doing gather.

# overlap     p_query p_match avg_abund
# ---------   ------- ------- ---------
# 80.0 kbp       0.0%    0.4%       1.1    FR798975.1 Leishmania braziliensis MHOM/BR/75/M2904 comp...
# 70.0 kbp       0.3%    1.2%      13.9    FLZI01001256.1 Escherichia coli isolate CG05C.C2 genome ...
# found less than 50.0 kbp in common. => exiting

# found 2 matches total;
# the recovered matches hit 0.3% of the abundance-weighted query.
# the recovered matches hit 0.5% of the query k-mers (unweighted).



# once this is done, try rerunning with an additional argument 
# to relax the threshold used for classification: '--threshold-bp 100'
# no output?
sourmash gather -k 31 "data/leishmania_revisited/$FASTQ_FILENAME.sig" $GENOME_FILE \
--threshold-bp 100
# == This is sourmash version 4.6.1. ==
# == Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

# selecting specified query k=31
# loaded query: SRR8668774.33 33/1... (k=31, DNA)
# --
# loaded 93249 total signatures from 1 locations.
# after selecting signatures compatible with search, 93249 remain.

# Starting prefetch sweep across databases.
# Found 6180 signatures via prefetch; now doing gather.

# overlap     p_query p_match avg_abund
# ---------   ------- ------- ---------
# 80.0 kbp       0.0%    0.4%       1.1    FR798975.1 Leishmania braziliensis MHOM/BR/75/M2904 complete genome, chromosome 1
# 70.0 kbp       0.3%    1.2%      13.9    FLZI01001256.1 Escherichia coli isolate CG05C.C2 genome assembly, contig: 20151013_...
# 30.0 kbp       1.8%    0.6%     213.3    JTEB01000001.1 Klebsiella pneumoniae strain XH212 contig-100_0, whole genome shotgu...
# 80.0 kbp       0.2%    0.1%      30.0    KB453278.1 Leishmania panamensis MHOM/COL/81/L13 unplaced genomic scaffold Scaffold...
# 60.0 kbp       0.1%    0.4%      16.0    FLZB01000512.1 Escherichia coli isolate CT67C.C1 genome assembly, contig: 20151360_...
# 50.0 kbp       0.0%    0.7%       2.5    FREM01000944.1 Enterococcus faecium isolate Hp_6-9 genome assembly, contig: NODE_94...
# 20.0 kbp       0.0%    0.1%       1.5    LQMU01000001.1 Euglena gracilis var. bacillaris strain SAG 1224-5/15 As3_contig0001...
# 20.0 kbp       0.9%    0.5%     164.0    LYHM01000001.1 Acinetobacter baumannii strain XH740 XH740_contig_1, whole genome sh...
# 20.0 kbp       0.0%    0.0%       1.5    LLKM01000001.1 Hyaloperonospora arabidopsidis isolate Noks1 scf_59133_1.contig_1, w...
# 20.0 kbp       4.8%    0.5%     879.0    JQCU01000001.1 Acinetobacter idrijaensis strain MII contig_1, whole genome shotgun ...
# 60.0 kbp       0.0%    0.2%       6.0    FLZM01001124.1 Escherichia coli isolate CG05C.C1 genome assembly, contig: 20151012_...
# 60.0 kbp       0.0%    0.2%      14.0    FLZL01001110.1 Escherichia coli isolate CG10C.C1 genome assembly, contig: 20151316_...
# 50.0 kbp       0.0%    0.1%       1.0    LN609229.1 Leishmania peruviana genome assembly Leishmania peruviana LEM-1537_V1, c...
# 10.0 kbp       0.0%    0.6%       3.0    LKBH01000001.1 Acidiplasma cupricumulans strain BH2 Contig1, whole genome shotgun s...
# 10.0 kbp       0.0%    0.2%       1.0    MEEO01000001.1 Nitrosomonadales bacterium SCN 54-20 ABS69_C0001, whole genome shotg...
# 10.0 kbp       0.0%    0.3%       1.0    JPSE01000001.1 Plasmodium falciparum 365.1 ctg7180000077325, whole genome shotgun s...
# 10.0 kbp       0.0%    0.2%       1.0    JYOF01000001.1 Pseudomonas sp. 2(2015) contig1, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.1%       1.0    ABGO01000922.1 Gemmata obscuriglobus UQM 2246 gcontig_1106221710807, whole genome s...
# 10.0 kbp       0.0%    0.4%       3.0    CP003283.1 Ornithobacterium rhinotracheale DSM 15997, complete genome
# 10.0 kbp       1.7%    0.2%     604.0    LYHI01000001.1 Acinetobacter baumannii strain XH744 XH744_contig_1, whole genome sh...
# found less than 100 bp  in common. => exiting

# found 20 matches total;
# the recovered matches hit 9.8% of the abundance-weighted query.
# the recovered matches hit 1.3% of the query k-mers (unweighted).



# ----------------------------------------------------------------------
# Task 3

export CENTRIFUGE_OUTPUT='CL13_dehosted'
conda activate 'rnaseq'

# time = ~10min
# If you haven't done so already, unzip the file you downloaded above to get a folder. 
# Move the three files present in this folder to your working directory where your dehosted fastq file is located.
# see: https://www.biostars.org/p/399916/

# make output_dir if not exists
if ! test -d "$DATA_DIR/centrifuge" ; then
    mkdir $DATA_DIR/centrifuge
    echo "TRUE"
fi
centrifuge -x "$DATA_DIR/reference_sequences/p_compressed+h+v/p_compressed+h+v" \
    -U $FASTQ_FILE \
    --report-file "$DATA_DIR/centrifuge/${CENTRIFUGE_OUTPUT}_report.txt" \
    -S "$DATA_DIR/centrifuge/${CENTRIFUGE_OUTPUT}_results.txt"
# once done running, open the 'CL13_dehosted_results.txt'
# ouput file in Excel and sort based on 'numUniqueReads' column

# report file CL13_dehosted_report.txt
# Number of iterations in EM algorithm: 53
# Probability diff. (P - P_prev) in the last iteration: 7.22319e-11
# Calculating abundance: 00:00:00



# ----------------------------------------------------------------------
# Task 4

# Run sourmash and centrifuge on a new fastQ file: SRR11207265_1.fastq.gz
conda activate "rnaseq"
export FASTQ_FILENAME="SRR11207265_1.fastq.gz"
export FASTQ_FILE="$DATA_DIR/fastq/$FASTQ_FILENAME"

sourmash sketch dna -p scaled=10000,k=31,abund $FASTQ_FILE --name-from-first
mv "$FASTQ_FILENAME.sig" "$DATA_DIR/$FASTQ_FILENAME.sig"

sourmash gather -k 31 "$DATA_DIR/$FASTQ_FILENAME.sig" $GENOME_FILE --threshold-bp 100
# == This is sourmash version 4.6.1. ==
# == Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

# selecting specified query k=31
# loaded query: SRR11207265.1 1/1... (k=31, DNA)
# --
# loaded 93249 total signatures from 1 locations.
# after selecting signatures compatible with search, 93249 remain.

# Starting prefetch sweep across databases.
# Found 37955 signatures via prefetch; now doing gather.

# overlap     p_query p_match avg_abund
# ---------   ------- ------- ---------
# 7.4 Mbp       24.3%   66.3%      38.0    JSPL01000060.1 Escherichia coli strain blood-2011-0167 blood-2011-0167_ctg_1010, whole genome shotgun sequence
# 4.8 Mbp        0.7%   51.8%       1.7    CM004719.1 Saccharomyces cerevisiae strain wild006 chromosome I, whole genome shotgun sequence
# 4.6 Mbp       18.3%   97.8%      47.1    MXVD01000010.1 Salmonella enterica strain BCW_3049 NODE_10_length_171219_cov_3.05473, whole genome shotgun sequence
# 4.2 Mbp        7.0%  100.0%      19.2    KN049967.1 Bacillus subtilis strain BST genomic scaffold scaffold1, whole genome shotgun sequence
# 5.4 Mbp       12.0%   44.2%      43.4    FRYH01002035.1 Escherichia coli isolate Hp_23-13 genome assembly, contig: NODE_2035_length_128_cov_49_ID_40819, whole genome shotgun sequence
# 2.9 Mbp        0.4%   18.7%       1.4    AE017341.1 Cryptococcus neoformans var. neoformans JEC21 chromosome 1, complete sequence
# 2.8 Mbp        4.5%   97.3%      18.5    MKPI01000001.1 Listeria monocytogenes strain NRRL B-33794 PROKKA_contig000001, whole genome shotgun sequence
# 2.5 Mbp        6.0%   93.4%      27.5    JJDH01000001.1 Staphylococcus aureus 11S01420 adZAT-supercont1.1, whole genome shotgun sequence
# 2.2 Mbp        0.2%   12.2%       1.4    CP003820.1 Cryptococcus neoformans var. grubii H99 chromosome 1, complete sequence
# 6.7 Mbp        5.0%   26.6%      32.6    AYUC01000001.1 Pseudomonas aeruginosa ATCC 15442 contig000001, whole genome shotgun sequence
# 1.6 Mbp        2.9%  100.0%      20.5    GG669973.1 Lactobacillus fermentum ATCC 14931 genomic scaffold SCAFFOLD74, whole genome shotgun sequence
# 2.5 Mbp        2.8%   31.2%      37.2    GG692897.1 Enterococcus faecalis DS5 plasmid scaffold supercont1.43, whole genome shotgun sequence
# 4.5 Mbp        2.5%   10.9%      45.4    JOST01000001.1 Escherichia coli 5-366-08_S4_C1 e536608S4C1.contig.0_1, whole genome shotgun sequence
# 4.6 Mbp        0.1%    3.6%       2.3    CM005535.1 Saccharomyces cerevisiae strain bioethanol001 chromosome I, whole genome shotgun sequence
# 4.2 Mbp        0.0%    1.5%       1.8    AGVY01000003.1 Saccharomyces cerevisiae x Saccharomyces kudriavzevii VIN7 chromosome ScI VIN7_c12_3, whole genome shotgun sequence
# 4.4 Mbp        0.6%    2.8%      43.5    CYUC01000001.1 Achromobacter sp. ATCC35328 genome assembly, contig: ERS372666SCcontig000001, whole genome shotgun sequence
# 4.4 Mbp        0.0%    0.5%       1.4    LOQK01000001.1 Saccharomyces pastorianus strain Hybrid yeast 7 sequence_1, whole genome shotgun sequence
# 2.9 Mbp        0.0%    0.6%       1.4    CM000040.1 Cryptococcus neoformans var. neoformans B-3501A chromosome 1, whole genome shotgun sequence
# 2.7 Mbp        0.3%    1.6%      32.4    MOKD01000001.1 Escherichia coli strain 679 BN4_679_1_(paired)_contig_1, whole genome shotgun sequence
# 1.6 Mbp        0.0%    1.3%       1.4    MXZI01000100.1 Salmonella enterica subsp. salamae serovar 55:k:z39 strain BCW_2788 NODE_100_length_16738_cov_1.30046, whole genome shotgun ...
# 1.5 Mbp        0.2%    3.2%      46.8    LJFJ01000001.1 Lactobacillus fermentum strain HFB3 contig01, whole genome shotgun sequence
# 1.1 Mbp        0.0%    1.2%       1.0    MXPA01000010.1 Salmonella enterica strain BCW_1529 NODE_10_length_150447_cov_1.41913, whole genome shotgun sequence
# 4.6 Mbp        0.0%    0.5%       1.0    CM005711.1 Saccharomyces cerevisiae strain beer092 chromosome I, whole genome shotgun sequence
# 4.4 Mbp        0.1%    1.0%      25.0    JMGV01000001.1 Escherichia coli 2-011-08_S3_C3 e201108S3C3.contig.0_1, whole genome shotgun sequence
# 4.4 Mbp        0.0%    0.6%       1.8    CP011810.1 Saccharomyces cerevisiae strain NCIM3186 chromosome I sequence
# 3.5 Mbp        0.0%    1.0%       7.2    KB732580.1 Escherichia coli KTE233 genomic scaffold acEoc-supercont1.1, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.9%       1.2    MXXT01001000.1 Salmonella enterica subsp. enterica serovar Senftenberg strain BCW_2841 NODE_1000_length_248_cov_1.74054, whole genome shotg...
# 2.1 Mbp        0.0%    1.0%       1.0    JSJW01000002.1 Escherichia coli strain upec-244 upec-244_ctg_1216, whole genome shotgun sequence
# 1.9 Mbp        0.0%    0.9%       1.0    AEJX01000517.1 Escherichia sp. TW15838 contig_517, whole genome shotgun sequence
# 1.8 Mbp        0.1%    1.7%      24.2    KB944701.1 Enterococcus faecalis RMC1 genomic scaffold aczIg-supercont1.1, whole genome shotgun sequence
# 1.5 Mbp        0.0%    1.0%       1.2    MXON01000010.1 Salmonella enterica subsp. houtenae serovar 45:g,z51:- strain BCW_1545 NODE_10_length_188647_cov_0.269975, whole genome shot...
# 5.3 Mbp        0.0%    0.6%       1.0    FRVO01000145.1 Pseudomonas aeruginosa isolate 437.1 genome assembly, contig: 437_1.NODE_317_length_5119_cov_222.715179, whole genome shotgu...
# 4.5 Mbp        0.0%    0.6%       1.2    LLOE01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07380 WH-SGI-V-07380_contig1, whole genome shotgun sequence
# 4.0 Mbp        0.1%    0.8%      36.2    AMSK01000001.1 Escherichia coli AD30 1.ECAD30.1_1, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.9%       1.5    AUYQ02000001.1 Salmonella enterica subsp. enterica serovar Abortusovis str. SS44 Contig_001, whole genome shotgun sequence
# 1.5 Mbp        0.0%    0.6%       1.2    JTNE01000001.1 Pseudomonas aeruginosa strain AZPAE15042 AZPAE15042_contig_1, whole genome shotgun sequence
# 1.4 Mbp        0.0%    0.8%       2.0    MLXV01000001.1 Salmonella enterica subsp. diarizonae serovar 60:r:e,n,x,z15 strain CFSAN044987 CFSAN044987_1, whole genome shotgun sequence
# 5.4 Mbp        0.0%    0.5%       1.7    KI519281.1 Pseudomonas aeruginosa CF18 genomic scaffold adgfD-supercont1.1, whole genome shotgun sequence
# 4.6 Mbp        0.0%    0.3%       1.3    CP004450.2 Saccharomyces cerevisiae YJM1208 chromosome I genomic sequence
# 4.5 Mbp        0.0%    0.3%       2.0    LCTD01000001.1 Saccharomyces cerevisiae strain 4124-S60 PE150_denovo_1, whole genome shotgun sequence
# 4.5 Mbp        0.0%    0.5%       1.0    LLLV01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07166 WH-SGI-V-07166_contig1, whole genome shotgun sequence
# 4.3 Mbp        0.0%    0.4%       1.3    AEEZ01000003.1 Saccharomyces cerevisiae FostersO chromosome I FOSTERSO_c1_3, whole genome shotgun sequence
# 3.9 Mbp        0.0%    0.5%       1.3    LLPH01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07420 WH-SGI-V-07420_contig1, whole genome shotgun sequence
# 3.8 Mbp        0.1%    0.5%      25.0    MOHR01000001.1 Escherichia coli strain 508 BN4_508_1_(paired)_contig_1, whole genome shotgun sequence
# 3.7 Mbp        0.0%    0.6%       1.3    MXLL01000100.1 Salmonella enterica subsp. enterica serovar Typhimurium strain INSP 24 NODE_100_length_1640_cov_1.32287, whole genome shotgu...
# 3.6 Mbp        0.0%    0.6%       2.0    LLYH01000001.1 Escherichia coli strain AY65, whole genome shotgun sequence
# 3.6 Mbp        0.0%    0.6%       1.0    MYMI01000100.1 Salmonella enterica strain BCW_4985 NODE_100_length_11668_cov_11.4306, whole genome shotgun sequence
# 3.6 Mbp        0.0%    0.6%      19.0    AQCV01000001.1 Escherichia coli MP021017.11 gecMP02101711.contig.0, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.5%       1.0    JHFJ01000001.1 Escherichia coli O26:H11 str. 2010C-3871 contig1, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.6%       1.3    MXJB01000010.1 Salmonella enterica subsp. enterica serovar Worthington strain BCW_2037 NODE_10_length_238565_cov_4.0092, whole genome shotg...
# 3.3 Mbp        0.0%    0.6%       2.3    CP019183.1 Salmonella enterica subsp. enterica serovar Mbandaka str. ATCC 51958, complete genome
# 3.0 Mbp        0.0%    0.6%       1.0    MYYI01000010.1 Salmonella enterica subsp. enterica serovar 9,12:l,z28:- strain 01-0173 NODE_10_length_202304_cov_4.84396, whole genome shot...
# 2.9 Mbp        0.0%    0.6%       1.0    LFIH01000001.1 Salmonella enterica subsp. enterica serovar Napoli strain SN310 SNOUT310_1, whole genome shotgun sequence
# 2.9 Mbp        0.0%    0.6%       1.3    JH590510.1 Escherichia coli E101 genomic scaffold supercont1.1, whole genome shotgun sequence
# 2.7 Mbp        0.0%    0.7%       2.0    AKNA01000001.1 Shigella boydii 965-58 gss96558.contig.0, whole genome shotgun sequence
# 2.6 Mbp        0.0%    0.5%       1.7    JSRU01000006.1 Escherichia coli strain blood-09-0023rep2 blood-09-0023rep2_ctg_1047, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.6%       1.0    AICE01000001.1 Escherichia coli C9_92 C9_92_2, whole genome shotgun sequence
# 1.9 Mbp        0.0%    0.7%       1.0    MWZB01000100.1 Bacillus intestinalis strain 1731 NODE_100_length_250_cov_0.855491_ID_199, whole genome shotgun sequence
# 1.8 Mbp        0.0%    0.6%       2.3    MXPJ01000010.1 Salmonella enterica subsp. salamae serovar 48:d:z6 strain BCW_1519 NODE_10_length_145457_cov_2.08435, whole genome shotgun s...
# 1.7 Mbp        0.0%    1.0%       1.0    KV805526.1 Enterococcus sp. HMSC078F03 genomic scaffold Scaffold0, whole genome shotgun sequence
# 0.8 Mbp        0.0%    0.6%       1.0    LFHY01000001.1 Escherichia sp. B1147 EC2821_c1, whole genome shotgun sequence
# 210.0 kbp      0.0%    0.7%       1.3    FAMN01000001.1 Peptoclostridium difficile isolate VL_0413 genome assembly, contig: NODE_1_length_322931_cov_27.6498_ID_1, whole genome shot...
# 80.0 kbp       0.0%    0.6%       2.7    JSXG01000100.1 Stenotrophomonas maltophilia strain B418 Contig_100_, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.3%       1.0    KE137323.1 Pseudomonas aeruginosa PAK genomic scaffold acUKM-supercont-complete, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.3%       1.0    LJNV01000001.1 Pseudomonas aeruginosa strain ATCC 33349 scaffold_0, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.3%       1.0    JTMU01000001.1 Pseudomonas aeruginosa strain AZPAE15052 AZPAE15052_contig_1, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.3%       1.0    ASQY01000001.1 Pseudomonas aeruginosa str. PA 17 contig001, whole genome shotgun sequence
# 5.0 Mbp        0.0%    0.3%       1.5    BAMB01001377.1 Pseudomonas aeruginosa JCM 6119 DNA, contig: JCM6119.contig01377, whole genome shotgun sequence
# 4.6 Mbp        0.0%    0.3%       1.5    CUXK01000001.1 Pseudomonas aeruginosa genome assembly 244_15_311Ga_contgs_500, contig NODE_50_length_4975_cov_37.4803_ID_99, whole genome s...
# 4.3 Mbp        0.0%    0.4%      27.5    LZDY01000001.1 Escherichia coli strain EPEC 1316 EPEC-1316_contig_1, whole genome shotgun sequence
# 4.2 Mbp        0.0%    0.2%       1.0    CP020225.1 Saccharomyces cerevisiae strain UWOPS03-461.4 chromosome I, complete sequence
# 4.1 Mbp        0.0%    0.4%       1.0    JNRV01000001.1 Escherichia coli 3-073-06_S1_C2 e307306S1C2.contig.0_1, whole genome shotgun sequence
# 4.0 Mbp        0.0%    0.4%       2.0    AQBI01000001.1 Escherichia coli P0305260.11 gecP030526011.contig.0, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.5%       1.0    AJHO01000001.1 Bacillus subtilis subsp. spizizenii RFWG5B15 contig00001, whole genome shotgun sequence
# 3.6 Mbp        0.0%    0.4%       1.0    MLWM01000001.1 Salmonella enterica subsp. enterica serovar Saintpaul strain CFSAN045043 CFSAN045043_1, whole genome shotgun sequence
# 3.6 Mbp        0.0%    0.4%       1.0    JH953706.1 Escherichia coli 3006 genomic scaffold E3006.contig.0, whole genome shotgun sequence
# 3.5 Mbp        0.0%    0.4%       1.0    MXSZ01000100.1 Salmonella enterica subsp. enterica serovar Heidelberg strain NCTR-SF606 NODE_100_length_15548_cov_5.15738, whole genome sho...
# 3.5 Mbp        0.0%    0.4%       1.0    MYMD01000100.1 Salmonella enterica strain BCW_4992 NODE_100_length_2112_cov_1.97506, whole genome shotgun sequence
# 3.5 Mbp        0.0%    0.4%       1.5    MYTO01000100.1 Salmonella enterica subsp. enterica serovar Newport strain BCW_4343 NODE_100_length_8164_cov_2.10449, whole genome shotgun s...
# 3.5 Mbp        0.0%    0.4%       1.5    APXA01000001.1 Escherichia coli 178900 gec178900.contig.0, whole genome shotgun sequence
# 3.5 Mbp        0.0%    0.4%       1.0    KK082444.1 Salmonella enterica subsp. enterica serovar Namur str. 05-2929 genomic scaffold scaffold10_size1049, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.4%       1.0    MYMN01001000.1 Salmonella enterica strain BCW_4959 NODE_1000_length_869_cov_41.8564, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.3%       1.0    LLLT01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07050 WH-SGI-V-07050_contig1, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.4%       1.0    AM933173.1 Salmonella enterica subsp. enterica serovar Gallinarum str. 287/91 complete genome
# 3.3 Mbp        0.0%    0.4%       2.0    MIVM01000001.1 Escherichia coli strain VL2832 contig_0001, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.4%       1.0    MYYO01000010.1 Salmonella enterica subsp. enterica serovar Corvallis strain 95-0387 NODE_10_length_155334_cov_3.84565, whole genome shotgun...
# 2.9 Mbp        0.0%    0.4%       1.0    MZAP01000010.1 Salmonella enterica subsp. enterica serovar Edinburgh strain 2012K-0034 NODE_10_length_159416_cov_3.13705, whole genome shot...
# 2.8 Mbp        0.0%    0.4%       1.0    LM997148.1 Escherichia coli genome assembly FHI92, scaffold scaffold-1_contig-23.0_1_1096_[organism:Escherichia
# 2.7 Mbp        0.0%    0.3%       1.0    MOKL01000001.1 Escherichia coli strain 614 BN4_614_1_(paired)_contig_1, whole genome shotgun sequence
# 2.6 Mbp        0.0%    0.2%       2.0    AKNF01000001.1 Shigella flexneri 1235-66 strain 1236-66 gss123566.contig.0, whole genome shotgun sequence
# 2.5 Mbp        0.0%    0.4%       1.0    MIIK01000001.1 Escherichia coli strain B-8638 B8638_contig1, whole genome shotgun sequence
# 2.5 Mbp        0.0%    0.4%       1.5    LVOA01000001.1 Escherichia coli strain swine50 swine50_contig_1, whole genome shotgun sequence
# 2.4 Mbp        0.0%    0.3%       1.0    CP015229.1 Escherichia coli strain 06-00048, complete genome
# 2.4 Mbp        0.0%    0.4%       1.5    CP010157.1 Escherichia coli strain D10, complete genome
# 2.3 Mbp        0.0%    0.5%       1.0    ADUT01000067.1 Shigella dysenteriae 1617 gss1617.assembly.9, whole genome shotgun sequence
# 2.3 Mbp        0.0%    0.4%       1.0    ADUL01000083.1 Escherichia coli 2362-75 gec2362.assembly.9, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.4%       1.0    JSOT01000004.1 Escherichia coli strain upec-107 upec-107_ctg_2184, whole genome shotgun sequence
# 2.0 Mbp        0.0%    0.8%       1.5    MAQO01000001.1 Staphylococcus aureus strain Sa14-001 Sa14-001_SPAdes_contigs__contig_1, whole genome shotgun sequence
# 2.0 Mbp        0.0%    0.4%       1.0    AOQI01000001.1 Escherichia coli TOP293-3, whole genome shotgun sequence
# 1.9 Mbp        0.0%    0.8%       1.5    GG668818.1 Enterococcus faecalis ATCC 29200 genomic scaffold SCAFFOLD64, whole genome shotgun sequence
# 1.9 Mbp        0.0%    0.4%       1.0    MXPF01000010.1 Salmonella enterica subsp. salamae serovar 13,22:z29:enx strain BCW_1523 NODE_10_length_143223_cov_0.212613, whole genome sh...
# 1.7 Mbp        0.0%    0.4%       1.0    LHTH01000001.1 Salmonella enterica subsp. salamae serovar 56:z10:e,n,x str. 1369-73 ssp-II_O56_mirahybrid1_c1, whole genome shotgun sequence
# 1.7 Mbp        0.0%    0.4%       2.0    MYRL01000010.1 Salmonella enterica strain BCW_4659 NODE_10_length_134779_cov_2.70265, whole genome shotgun sequence
# 1.6 Mbp        0.0%    0.8%       1.5    FMPJ01000048.1 Staphylococcus aureus strain GKP136-11 genome assembly, contig: ERS189597SCcontig000048, whole genome shotgun sequence
# 1.5 Mbp        0.0%    0.7%       1.0    MKYZ01000010.1 Staphylococcus aureus strain SA_123 NODE_10_length_38017_cov_53.1023_ID_1585, whole genome shotgun sequence
# 1.2 Mbp        0.0%    0.4%       1.0    MXOI01000100.1 Salmonella enterica strain BCW_1551 NODE_100_length_446_cov_0.663912, whole genome shotgun sequence
# 1.1 Mbp        0.0%    0.6%       1.0    LUEQ01000001.1 Listeria monocytogenes strain CFSAN026587 CFSAN026587_contig0000, whole genome shotgun sequence
# 0.9 Mbp        0.0%    0.4%       1.5    AEJW01000495.1 Escherichia sp. TW09231 contig_495, whole genome shotgun sequence
# 0.6 Mbp        0.0%    0.4%       1.0    MXOH01000010.1 Salmonella bongori serovar 48:i:- strain BCW_1552 NODE_10_length_195952_cov_0.253221, whole genome shotgun sequence
# 220.0 kbp      0.0%    0.4%       1.5    CACD01000357.1 Citrobacter freundii str. ballerup 7851/39, whole genome shotgun sequence, contig00357
# 190.0 kbp      0.0%    0.4%       1.5    FABQ01000001.1 Peptoclostridium difficile isolate VL_0122 genome assembly, contig: NODE_1_length_100232_cov_31.9568_ID_1, whole genome shot...
# 130.0 kbp      0.0%    0.4%       1.0    LSME01000001.1 Kluyvera ascorbata strain WCH1410 NODE_1, whole genome shotgun sequence
# 120.0 kbp      0.0%    0.1%       9.0    LN625280.1 Cryptococcus gattii R265 genomic scaffold, contig_1, whole genome shotgun sequence
# 120.0 kbp      0.0%    0.1%       7.0    CP000286.1 Cryptococcus gattii WM276 chromosome A, complete sequence
# 120.0 kbp      0.0%    0.3%       1.5    BATO01000237.1 Pseudomonas alcaligenes MRY13-0052 DNA, contig: contig00237, whole genome shotgun sequence
# 100.0 kbp      0.0%    0.5%       2.0    FACE01000001.1 Peptoclostridium difficile isolate VL_0117 genome assembly, contig: NODE_1_length_15690_cov_11.3658_ID_1, whole genome shotg...
# 90.0 kbp       0.0%    0.5%       2.5    JMSQ01000001.1 Siccibacter colletis strain 1383 contig001, whole genome shotgun sequence
# 60.0 kbp       0.0%    0.3%       1.0    CP018743.1 Pseudomonas sp. DRA525 genome
# 40.0 kbp       0.0%    0.3%       1.0    CIDK01000001.1 Streptococcus pneumoniae genome assembly 9803_7#73, scaffold ERS232502SCcontig000001, whole genome shotgun sequence
# 5.5 Mbp        0.0%    0.1%       1.0    KI914476.1 Pseudomonas aeruginosa BWHPSA045 genomic scaffold adkfU-supercont1.1, whole genome shotgun sequence
# 5.5 Mbp        0.0%    0.1%       1.0    LLNJ01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07238 WH-SGI-V-07238_contig1, whole genome shotgun sequence
# 5.4 Mbp        0.0%    0.2%       1.0    JUAL01000001.1 Pseudomonas aeruginosa strain AZPAE12151 AZPAE12151_contig_1, whole genome shotgun sequence
# 5.4 Mbp        0.0%    0.1%       1.0    CTFM01000001.1 Pseudomonas aeruginosa genome assembly 294_Rw130_contgs_500, contig NODE_51_length_9078_cov_35.5116_ID_101, whole genome sho...
# 5.4 Mbp        0.0%    0.1%       1.0    LLLO01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07051 WH-SGI-V-07051_contig1, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.2%       1.0    KV829825.1 Pseudomonas sp. HMSC060G02 genomic scaffold Scaffold0, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.1%       1.0    LJZA01000001.1 Pseudomonas aeruginosa strain ATCC 33354 scaffold_0, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.1%       1.0    JUAP01000001.1 Pseudomonas aeruginosa strain AZPAE12147 AZPAE12147_contig_1, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.2%       4.0    LLUZ01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07310 WH_SGI_V_07310_1, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.2%       1.0    KI518732.1 Pseudomonas aeruginosa M8A.4 genomic scaffold adgfJ-supercont1.1, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.2%       1.0    AJXX01000001.1 Pseudomonas aeruginosa XMG contig000001, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.2%       1.0    LLPX01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07491 WH-SGI-V-07491_contig1, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.2%       1.0    LBIG01000001.1 Pseudomonas aeruginosa MRSN 321 contig1, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.2%       5.0    JTQU01000001.1 Pseudomonas aeruginosa strain AZPAE14947 AZPAE14947_contig_1, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.2%       2.0    JTUZ01000001.1 Pseudomonas aeruginosa strain AZPAE14835 AZPAE14835_contig_1, whole genome shotgun sequence
# 5.3 Mbp        0.0%    0.2%       1.0    LLMH01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07178 WH-SGI-V-07178_contig1, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.2%       2.0    LLQR01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07620 WH-SGI-V-07620_contig1, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.2%       1.0    LLMX01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07226 WH-SGI-V-07226_contig1, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.2%       1.0    JTZS01000001.1 Pseudomonas aeruginosa strain AZPAE12422 AZPAE12422_contig_1, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.2%       1.0    JTRH01000001.1 Pseudomonas aeruginosa strain AZPAE14934 AZPAE14934_contig_1, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.2%       1.0    AWYK01000001.1 Pseudomonas aeruginosa JD318, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.1%       1.0    JTWO01000001.1 Pseudomonas aeruginosa strain AZPAE14718 AZPAE14718_contig_1, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.1%       1.0    JTZW01000001.1 Pseudomonas aeruginosa strain AZPAE12418 AZPAE12418_contig_1, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.2%       1.0    LJZC01000001.1 Pseudomonas aeruginosa strain ATCC 33356 scaffold_0, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.2%       1.0    JTTD01000001.1 Pseudomonas aeruginosa strain AZPAE14885 AZPAE14885_contig_1, whole genome shotgun sequence
# 5.2 Mbp        0.0%    0.2%       1.0    FRRY01000054.1 Pseudomonas aeruginosa isolate 378 genome assembly, contig: 378.NODE_125_length_4338_cov_85.702629, whole genome shotgun seq...
# 5.2 Mbp        0.0%    0.2%       1.0    LLKX01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07053 WH-SGI-V-07053_contig1, whole genome shotgun sequence
# 5.1 Mbp        0.0%    0.2%       1.0    JTYM01000001.1 Pseudomonas aeruginosa strain AZPAE14410 AZPAE14410_contig_1, whole genome shotgun sequence
# 5.1 Mbp        0.0%    0.2%       2.0    MAZW01000023.1 Pseudomonas aeruginosa strain TNCF_106 NODE_10_length_205191_cov_15.4706_ID_19, whole genome shotgun sequence
# 5.1 Mbp        0.0%    0.2%       1.0    JTQB01000001.1 Pseudomonas aeruginosa strain AZPAE14967 AZPAE14967_contig_1, whole genome shotgun sequence
# 5.0 Mbp        0.0%    0.2%       1.0    JTRL01000001.1 Pseudomonas aeruginosa strain AZPAE14930 AZPAE14930_contig_1, whole genome shotgun sequence
# 4.9 Mbp        0.0%    0.2%       1.0    KK213160.1 Pseudomonas aeruginosa PS75 genomic scaffold adTyT-supercont1.1, whole genome shotgun sequence
# 4.9 Mbp        0.0%    0.2%       1.0    LCSQ01000001.1 Pseudomonas aeruginosa strain Pae_CF67.01g CF67.01g_contig_1, whole genome shotgun sequence
# 4.8 Mbp        0.0%    0.1%       1.0    MPCQ01000001.1 Pseudomonas aeruginosa strain PAC08 scaffold1-size527348, whole genome shotgun sequence
# 4.7 Mbp        0.0%    0.1%       1.0    JTYN01000001.1 Pseudomonas aeruginosa strain AZPAE14404 AZPAE14404_contig_1, whole genome shotgun sequence
# 4.7 Mbp        0.0%    0.2%       1.0    LRYG01000001.1 Pseudomonas aeruginosa strain AU10714 NODE_103_length_2700_cov_16.0163_ID_205, whole genome shotgun sequence
# 4.7 Mbp        0.0%    0.1%       1.0    FRKW01000131.1 Pseudomonas aeruginosa isolate 177 genome assembly, contig: 177.NODE_312_length_592_cov_125.832771, whole genome shotgun seq...
# 4.6 Mbp        0.0%    0.1%       1.0    LLOI01000001.1 Pseudomonas aeruginosa strain WH-SGI-V-07384 WH-SGI-V-07384_contig1, whole genome shotgun sequence
# 4.6 Mbp        0.0%    0.1%       2.0    CP008236.1 Saccharomyces cerevisiae strain HB_S_GIMBLETTROAD_16 chromosome I sequence
# 4.6 Mbp        0.0%    0.1%       2.0    CM005519.1 Saccharomyces cerevisiae strain bioethanol002 chromosome I, whole genome shotgun sequence
# 4.6 Mbp        0.0%    0.2%       3.0    KI914485.1 Pseudomonas aeruginosa BWHPSA043 genomic scaffold adkgk-supercont1.1, whole genome shotgun sequence
# 4.5 Mbp        0.0%    0.1%       1.0    DF396881.1 Saccharomyces cerevisiae NAM34-4C DNA, contig: scaffold001, whole genome shotgun sequence
# 4.5 Mbp        0.0%    0.2%       1.0    CM001062.1 Salmonella enterica subsp. enterica serovar Choleraesuis str. A50 chromosome, whole genome shotgun sequence
# 4.5 Mbp        0.0%    0.1%       1.0    KB822770.1 Saccharomyces cerevisiae M3838 unplaced genomic scaffold C398scaffold_1, whole genome shotgun sequence
# 4.5 Mbp        0.0%    0.1%       5.0    AHZG01000001.1 Saccharomyces cerevisiae Lalvin L2056 denovo_c1, whole genome shotgun sequence
# 4.5 Mbp        0.0%    0.1%       3.0    CM006511.1 Saccharomyces cerevisiae strain beer042 chromosome I, whole genome shotgun sequence
# 4.5 Mbp        0.0%    0.1%       2.0    CM004767.1 Saccharomyces cerevisiae strain wild001 chromosome I, whole genome shotgun sequence
# 4.5 Mbp        0.0%    0.1%       1.0    JSAC01000001.1 Saccharomyces cerevisiae strain 9464 utg7180000000000, whole genome shotgun sequence
# 4.5 Mbp        0.0%    0.1%       1.0    FUGS01001556.1 Saccharomyces cerevisiae isolate MTF2119 genome assembly, contig: SACE_Genowine_1014-F5.C15342, whole genome shotgun sequence
# 4.5 Mbp        0.0%    0.1%       1.0    CP008117.1 Saccharomyces cerevisiae strain MTKSKsf_E2 chromosome I sequence
# 4.4 Mbp        0.0%    0.2%       1.0    LGMR01000005.1 Escherichia coli strain 401900 NODE10, whole genome shotgun sequence
# 4.4 Mbp        0.0%    0.1%       1.0    JPXA01000001.1 Saccharomyces cerevisiae M5 contig_5, whole genome shotgun sequence
# 4.4 Mbp        0.0%    0.1%       2.0    FUGP01000721.1 Saccharomyces cerevisiae isolate MTF2114 genome assembly, contig: SACE_Genowine_FS2D.C101158, whole genome shotgun sequence
# 4.4 Mbp        0.0%    0.2%       1.0    AP012030.1 Escherichia coli DH1 (ME8569) DNA, complete genome
# 4.4 Mbp        0.0%    0.2%       1.0    AIFV01000061.1 Escherichia coli DEC6A gecDEC6A.contig.60_1, whole genome shotgun sequence
# 4.3 Mbp        0.0%    0.1%       5.0    CP004477.1 Saccharomyces cerevisiae YJM1401 chromosome I genomic sequence
# 4.3 Mbp        0.0%    0.2%       1.0    AIFZ01000053.1 Escherichia coli DEC6E gecDEC6E.contig.52_1, whole genome shotgun sequence
# 4.3 Mbp        0.0%    0.2%       2.0    JNRH01000001.1 Escherichia coli 2-474-04_S3_C1 e247404S3C1.contig.0_1, whole genome shotgun sequence
# 4.3 Mbp        0.0%    0.2%       2.0    CP010152.1 Escherichia coli strain D9, complete genome
# 4.3 Mbp        0.0%    0.2%       3.0    KE137010.1 Escherichia coli KTE61 genomic scaffold acHaS-supercont1.1, whole genome shotgun sequence
# 4.3 Mbp        0.0%    0.2%       2.0    MOYZ01000100.1 Escherichia coli strain 1-R1 NODE_100_length_10199_cov_173.972, whole genome shotgun sequence
# 4.3 Mbp        0.0%    0.2%       1.0    JSNY01000002.1 Escherichia coli strain upec-127 upec-127_ctg_1139, whole genome shotgun sequence
# 4.2 Mbp        0.0%    0.1%       1.0    LRVB01000001.1 Saccharomyces sp. 'boulardii' strain ATCC MYA-797 contig_1, whole genome shotgun sequence
# 4.2 Mbp        0.0%    0.2%       1.0    CYCG01000001.1 Escherichia coli genome assembly 7790_1#90, scaffold ERS085423SCcontig000001, whole genome shotgun sequence
# 4.2 Mbp        0.0%    0.2%       1.0    LOPR01000001.1 Escherichia coli strain G76 Escherichia_coli_G76.contig.001, whole genome shotgun sequence
# 4.2 Mbp        0.0%    0.2%       1.0    LVLP01000001.1 Escherichia coli strain cattle6 cattle6_contig_1, whole genome shotgun sequence
# 4.2 Mbp        0.0%    0.2%       1.0    LOOK01000001.1 Escherichia coli strain G220 Escherichia_coli_G220.contig.001, whole genome shotgun sequence
# 4.1 Mbp        0.0%    0.2%       2.0    LOOI01000001.1 Escherichia coli strain G214 Escherichia_coli_G214.contig.001, whole genome shotgun sequence
# 4.1 Mbp        0.0%    0.2%       1.0    LONZ01000001.1 Escherichia coli strain G190 Escherichia_coli_G190.contig.001, whole genome shotgun sequence
# 4.1 Mbp        0.0%    0.2%       2.0    KE702342.1 Escherichia coli UMEA 3318-1 genomic scaffold acYxB-supercont1.1, whole genome shotgun sequence
# 4.1 Mbp        0.0%    0.2%       1.0    KQ087687.1 Escherichia coli strain MGH108 genomic scaffold aetAM-supercont1.1, whole genome shotgun sequence
# 4.0 Mbp        0.0%    0.2%       4.0    AQBJ01000001.1 Escherichia coli P0305260.12 gecP030526012.contig.0, whole genome shotgun sequence
# 4.0 Mbp        0.0%    0.2%       1.0    MOZI01000100.1 Escherichia coli strain 1.2-R3 NODE_100_length_2630_cov_37.2526, whole genome shotgun sequence
# 3.9 Mbp        0.0%    0.2%       1.0    JPKK01000001.1 Escherichia coli strain G5 EL78_contig_1, whole genome shotgun sequence
# 3.9 Mbp        0.0%    0.2%       1.0    JHDG01000001.1 Escherichia coli 1-176-05_S3_C1 e117605S3C1.contig.0_1, whole genome shotgun sequence
# 3.9 Mbp        0.0%    0.2%       4.0    MOGN01000001.1 Escherichia coli strain 554 BN4_554_1_(paired)_contig_1, whole genome shotgun sequence
# 3.9 Mbp        0.0%    0.2%       1.0    JMJD01000001.1 Escherichia coli strain UCD_JA65_pb contig_1, whole genome shotgun sequence
# 3.9 Mbp        0.0%    0.2%       1.0    CP017631.1 Escherichia coli SLK172, complete genome
# 3.9 Mbp        0.0%    0.1%       2.0    GL875991.1 Saccharomyces cerevisiae FL100 unplaced genomic scaffold Scaffold0, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.2%       2.0    AQCZ01000001.1 Escherichia coli BCE019_MS-13 gecBCE019MS13.contig.0, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.2%       1.0    APZM01000001.1 Escherichia coli BCE006_MS-23 gecBCE006MS23.contig.0, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.2%       1.0    ASHB01000001.1 Escherichia coli O157 str. NCCP15738 Scaffold1_1, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.2%       1.0    MOKC01000001.1 Escherichia coli strain 682 BN4_682_1_(paired)_contig_1, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.2%       1.0    JTQK01000001.1 Pseudomonas aeruginosa strain AZPAE14957 AZPAE14957_contig_1, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.2%       2.0    KE702490.1 Escherichia coli UMEA 3682-1 genomic scaffold acYzk-supercont1.1, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.2%       1.0    AQDE01000001.1 Escherichia coli 2871950 gec2871950.contig.0, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.2%       1.0    JNRP01000001.1 Escherichia coli 3-020-07_S1_C2 e302007S1C2.contig.0_1, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.2%       1.0    KI522465.1 Escherichia coli 909945-2 genomic scaffold Scaffold0, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.2%       1.0    CXYA01000001.1 Escheichia coli genome assembly 1.EC2923.1, contig P816_1, whole genome shotgun sequence
# 3.8 Mbp        0.0%    0.2%       1.0    JHRM01000001.1 Escherichia coli 401140 401140_1, whole genome shotgun sequence
# 3.7 Mbp        0.0%    0.2%       5.0    MOJL01000001.1 Escherichia coli strain 499 BN4_499_1_(paired)_contig_1, whole genome shotgun sequence
# 3.7 Mbp        0.1%    0.2%      67.0    MOFU01000001.1 Escherichia coli strain 458 BN4_458_1_(paired)_contig_1, whole genome shotgun sequence
# 3.7 Mbp        0.0%    0.2%       1.0    LVNO01000001.1 Escherichia coli strain swine38 swine38_contig_1, whole genome shotgun sequence
# 3.7 Mbp        0.0%    0.2%       1.0    CXWS01000001.1 Escheichia coli genome assembly YE9, contig YE9_1, whole genome shotgun sequence
# 3.7 Mbp        0.0%    0.2%       1.0    MXXQ01000010.1 Salmonella enterica subsp. enterica serovar Bovismorbificans strain BCW_2845 NODE_10_length_155790_cov_2.37306, whole genome...
# 3.6 Mbp        0.0%    0.2%       2.0    LVGJ01000001.1 Salmonella enterica subsp. enterica serovar Typhimurium strain CFSAN033870 CFSAN033870_contig0000, whole genome shotgun sequ...
# 3.6 Mbp        0.0%    0.2%       1.0    LYRS01000412.1 Salmonella enterica subsp. enterica serovar Typhimurium strain ST452 plasmid unnamed ST452-Typhimurium_contig_49, whole geno...
# 3.6 Mbp        0.0%    0.2%       3.0    LVFS01000001.1 Salmonella enterica subsp. enterica serovar Typhimurium strain CFSAN033896 CFSAN033896_contig0000, whole genome shotgun sequ...
# 3.6 Mbp        0.0%    0.2%       1.0    MZFM01000010.1 Salmonella enterica subsp. enterica serovar Reading strain 97-0463 NODE_10_length_174941_cov_3.53836, whole genome shotgun s...
# 3.5 Mbp        0.0%    0.2%       1.0    MYCG01000010.1 Salmonella enterica subsp. enterica serovar Isangi strain BCW_2703 NODE_10_length_192175_cov_4.99584, whole genome shotgun s...
# 3.5 Mbp        0.0%    0.2%       1.0    CCBJ010000001.1 Salmonella enterica subsp. enterica serovar Manhattan WGS project CCBJ00000000 data, isolate 165051_5, contig 165051_5_1, w...
# 3.5 Mbp        0.0%    0.2%       1.0    MXZM01000010.1 Salmonella enterica subsp. enterica serovar Stanley strain BCW_2782 NODE_10_length_175894_cov_4.12484, whole genome shotgun ...
# 3.5 Mbp        0.0%    0.2%       1.0    AOQT01000001.1 Escherichia coli TOP2515, whole genome shotgun sequence
# 3.5 Mbp        0.0%    0.2%       1.0    LBBE01000001.1 Escherichia coli strain avian5 contig00001, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    CU928145.2 Escherichia coli 55989 chromosome, complete genome
# 3.4 Mbp        0.0%    0.2%       1.0    MZEU01000010.1 Salmonella enterica subsp. enterica serovar Uganda strain 97-0352 NODE_10_length_194042_cov_3.05646, whole genome shotgun se...
# 3.4 Mbp        0.0%    0.2%       2.0    LKIX01000001.1 Salmonella enterica subsp. enterica serovar Derby strain 2012CEB117SAL scaffold_0, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    MYWN01000010.1 Salmonella enterica strain BCW_6260 NODE_10_length_197223_cov_3.04389, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       2.0    MXYU01000100.1 Salmonella enterica subsp. enterica serovar Vinnorrhady strain BCW_2812 NODE_100_length_15785_cov_5.12718, whole genome shot...
# 3.4 Mbp        0.0%    0.2%       1.0    GG771868.1 Escherichia coli MS 115-1 genomic scaffold Scfld1168, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    JHMT01000001.1 Escherichia coli O111:NM str. 2010C-4746 contig1, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    CCPT01000001.1 Escherichia coli genome assembly FHI3, contig scaffold-1_contig-1.0_1_1279_[organism:Escherichia, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    JWYX01000001.1 Escherichia coli strain CVM N38428PS N38428PS_contig_1, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    MXVZ01000010.1 Salmonella enterica subsp. enterica serovar Amsterdam strain BCW_2891 NODE_10_length_163901_cov_3.8342, whole genome shotgun...
# 3.4 Mbp        0.0%    0.2%       1.0    JUCX01000001.1 Escherichia coli strain CVM N34351PS N34351PS_contig_1, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    LM995789.1 Escherichia coli genome assembly FHI20, scaffold scaffold-1_contig-25.0_1_1204_[organism:Escherichia
# 3.4 Mbp        0.0%    0.2%       1.0    LWQZ01000001.1 Escherichia coli strain BJ10 Contig_1, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    MYAC01000100.1 Salmonella enterica subsp. enterica serovar Potsdam strain BCW_2762 NODE_100_length_14541_cov_1.32904, whole genome shotgun ...
# 3.4 Mbp        0.0%    0.2%       1.0    CP012344.1 Salmonella enterica subsp. enterica serovar Choleraesuis str. ATCC 10708, complete genome
# 3.4 Mbp        0.0%    0.2%       1.0    MZBA01000100.1 Salmonella enterica subsp. enterica serovar Cubana strain 07-0717 NODE_100_length_274_cov_10.2461, whole genome shotgun sequ...
# 3.4 Mbp        0.0%    0.2%       1.0    CYFQ01000001.1 Escherichia coli genome assembly 8205_3#8, scaffold ERS139205SCcontig000001, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    LVPP01000001.1 Escherichia coli strain sheep20 sheep20_contig_1, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    LRKB01000001.1 Escherichia coli strain 202733 NODE_100_length_1047_cov_84.4638_ID_199, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    KE702261.1 Escherichia coli UMEA 3292-1 genomic scaffold acXut-supercont1.1, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    CCQO01000001.1 Escherichia coli genome assembly FHI67, contig scaffold-1_contig-1.0_1_2047_[organism:Escherichia, whole genome shotgun sequ...
# 3.4 Mbp        0.0%    0.2%       1.0    MYYE01000010.1 Salmonella enterica subsp. enterica serovar Meleagridis strain 03-0344 NODE_10_length_163263_cov_3.04593, whole genome shotg...
# 3.4 Mbp        0.0%    0.2%       1.0    LRKD01000001.1 Escherichia coli strain 700328 ctg7180000028696, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    LM996928.1 Escherichia coli genome assembly FHI85, scaffold scaffold-1_contig-32.0_1_1231_[organism:Escherichia
# 3.4 Mbp        0.0%    0.2%       2.0    LOBH01000111.1 Salmonella enterica subsp. enterica serovar Kentucky strain CFSAN011777 contig_1, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       2.0    AQBD01000001.1 Escherichia coli P0304816.6 gecP03048166.contig.0, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    LHCW01000001.1 Escherichia coli strain CFSAN026849 contig_1, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    CXXT01000001.1 Escheichia coli genome assembly 1.EC2997.1, contig E955_1, whole genome shotgun sequence
# 3.4 Mbp        0.0%    0.2%       1.0    LVPD01000002.1 Escherichia coli strain sheep8 sheep8_contig_1, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       1.0    LOIH01000001.1 Escherichia coli strain STEC 2573 STEC2573_contig_1, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       1.0    AQAJ01000001.1 Escherichia coli P0304777.13 gecP030477713.contig.0, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       1.0    AQAT01000001.1 Escherichia coli P0304816.10 gecP030481610.contig.0, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       1.0    MPUE01000018.1 Escherichia coli strain TS21/08 IMT32647_S8_L001_R1_001_contig_1, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       2.0    LM997282.1 Escherichia coli genome assembly FHI90, scaffold scaffold-1_contig-29.0_1_1139_[organism:Escherichia
# 3.3 Mbp        0.0%    0.2%       2.0    BCOI01000001.1 Salmonella enterica DNA, contig: contig000001, strain: NGUA10, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       2.0    APZC01000001.1 Escherichia coli P0302308.13 gecP030230813.contig.0, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       2.0    JHSU01000001.1 Escherichia coli strain 103338 103338_1, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       1.0    LHGD01000001.1 Salmonella enterica subsp. enterica serovar Livingstone strain CVM N45399 N45399_contig_1, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       1.0    AFCJ01002683.1 Salmonella enterica subsp. enterica serovar Alachua str. R6-377 Contig2683, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       1.0    MYZP01000010.1 Salmonella enterica subsp. enterica serovar Cubana strain 49-1017 NODE_10_length_180267_cov_3.89656, whole genome shotgun se...
# 3.3 Mbp        0.0%    0.2%       3.0    APZF01000001.1 Escherichia coli P0302308.3 gecP03023083.contig.0, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       1.0    AIBV01000001.1 Escherichia coli C807_09 C807_09_2, whole genome shotgun sequence
# 3.3 Mbp        0.0%    0.2%       3.0    LHLZ01000001.1 Salmonella enterica subsp. enterica serovar Liverpool strain CVM N50444 N50444_contig_1, whole genome shotgun sequence
# 3.2 Mbp        0.0%    0.2%       1.0    LHTC01000001.1 Salmonella enterica subsp. enterica serovar Milwaukee str. SA19950795 S_Milwaukee_mirahybrid1_c1, whole genome shotgun sequence
# 3.2 Mbp        0.0%    0.2%       1.0    MYCT01000010.1 Salmonella enterica subsp. enterica serovar Hadar strain BCW_2689 NODE_10_length_154511_cov_3.18601, whole genome shotgun se...
# 3.2 Mbp        0.0%    0.2%       1.0    APZN01000001.1 Escherichia coli MP020980.1 gecMP0209801.contig.0, whole genome shotgun sequence
# 3.2 Mbp        0.0%    0.1%       1.0    GL872452.1 Saccharomyces cerevisiae Y10 unplaced genomic scaffold Scfld0, whole genome shotgun sequence
# 3.2 Mbp        0.0%    0.2%       1.0    MYBC01000100.1 Salmonella enterica subsp. enterica serovar Monschaui strain BCW_2734 NODE_100_length_287_cov_1.825, whole genome shotgun se...
# 3.2 Mbp        0.0%    0.2%       1.0    MYLC01000100.1 Salmonella enterica strain BCW_5090 NODE_100_length_1030_cov_0.866384, whole genome shotgun sequence
# 3.2 Mbp        0.0%    0.2%       2.0    APZL01000001.1 Escherichia coli 2854350 gec2854350.contig.0, whole genome shotgun sequence
# 3.2 Mbp        0.0%    0.2%       3.0    CP013029.1 Escherichia coli strain 2012C-4227, complete genome
# 3.2 Mbp        0.0%    0.2%       1.0    AIHN01000076.1 Escherichia coli DEC14D gecDEC14D.contig.76_1, whole genome shotgun sequence
# 3.2 Mbp        0.0%    0.2%       1.0    JHNY01000001.1 Escherichia coli O28ac:NM str. 02-3404 contig1, whole genome shotgun sequence
# 3.2 Mbp        0.0%    0.2%       1.0    CP019410.1 Salmonella enterica subsp. enterica serovar Hillingdon str. N1529-D3, complete genome
# 3.1 Mbp        0.0%    0.2%       1.0    MXMX01000010.1 Salmonella enterica subsp. enterica serovar Tornow strain BCW_1594 NODE_10_length_199137_cov_2.7807, whole genome shotgun se...
# 3.1 Mbp        0.0%    0.2%       2.0    MZDZ01001000.1 Salmonella enterica subsp. enterica serovar Kiambu strain 84-0491 NODE_1000_length_252_cov_0.630303, whole genome shotgun se...
# 3.1 Mbp        0.0%    0.2%       1.0    CCUQ01000001.1 Salmonella enterica subsp. enterica serovar Paratyphi A genome assembly PA040, contig 1, whole genome shotgun sequence
# 3.0 Mbp        0.0%    0.2%       3.0    MZDK01000010.1 Salmonella enterica subsp. enterica serovar Telelkebir strain 2011K-0078 NODE_10_length_150307_cov_3.11692, whole genome sho...
# 3.0 Mbp        0.0%    0.2%       1.0    CJRX01000001.1 Salmonella enterica subsp. enterica serovar Typhi genome assembly 8616_4#58, scaffold ERS168236SCcontig000001, whole genome ...
# 3.0 Mbp        0.0%    0.2%       2.0    FTUB01000383.1 Shigella sonnei strain 2090STDY5461813 genome assembly, contig: ERS219801SCcontig000383, whole genome shotgun sequence
# 3.0 Mbp        0.0%    0.2%       1.0    CVKD01000001.1 Salmonella enterica subsp. enterica serovar Typhi genome assembly 9475_4#50, scaffold ERS206728SCcontig000001, whole genome ...
# 3.0 Mbp        0.0%    0.2%       1.0    JUCI01000001.1 Escherichia coli strain CVM N36609PS N36609PS_contig_1, whole genome shotgun sequence
# 3.0 Mbp        0.0%    0.2%       1.0    AFZA01000001.1 Salmonella enterica subsp. enterica serovar Paratyphi A str. ZJ98-53 contig01, whole genome shotgun sequence
# 3.0 Mbp        0.0%    0.2%       1.0    LFCE01000001.1 Salmonella enterica subsp. enterica strain F12M01383-1 FDLB_F12M01383-1_contig1, whole genome shotgun sequence
# 3.0 Mbp        0.0%    0.2%       1.0    JQVX01000001.1 Salmonella enterica strain CRJJGF_00093 SQ0141Contig0, whole genome shotgun sequence
# 2.9 Mbp        0.0%    0.2%       1.0    MYIB01000001.1 Salmonella enterica strain BCW_5867 NODE_1_length_5770_cov_4.24477, whole genome shotgun sequence
# 2.9 Mbp        0.0%    0.2%       1.0    JH958332.1 Escherichia coli 5905 genomic scaffold E5905.contig.1, whole genome shotgun sequence
# 2.9 Mbp        0.0%    0.2%       1.0    AERN01000110.1 Shigella boydii ATCC 9905 SGB.Contig217_1, whole genome shotgun sequence
# 2.9 Mbp        0.0%    0.2%       1.0    AFCW01002504.1 Salmonella enterica subsp. enterica serovar Urbana str. R8-2977 Contig2504, whole genome shotgun sequence
# 2.8 Mbp        0.0%    0.2%       4.0    AKNB01000001.1 Shigella boydii 4444-74 gss444474.contig.0, whole genome shotgun sequence
# 2.8 Mbp        0.0%    0.2%       1.0    ANAN01000001.1 Shigella flexneri S6162 S6162_10, whole genome shotgun sequence
# 2.8 Mbp        0.0%    0.2%       1.0    MSJZ02000001.1 Shigella flexneri strain BCW_4876 BCW_4876_1__paired__contig_1, whole genome shotgun sequence
# 2.8 Mbp        0.0%    0.2%       2.0    LGBB01000001.1 Escherichia coli strain STEC 994 contig_1, whole genome shotgun sequence
# 2.8 Mbp        0.0%    0.2%       1.0    LPTQ01000009.1 Shigella boydii strain 600690 600690_10, whole genome shotgun sequence
# 2.8 Mbp        0.0%    0.2%       1.0    AMKE01000009.1 Shigella boydii 09-0344 SB_08_0344_10, whole genome shotgun sequence
# 2.8 Mbp        0.0%    0.2%       1.0    JXWI01000001.1 Escherichia coli strain OLC-975 Cont0001, whole genome shotgun sequence
# 2.7 Mbp        0.0%    0.2%       1.0    MSJY02000001.1 Shigella flexneri strain BCW_4875 BCW_4875_1__paired__contig_1, whole genome shotgun sequence
# 2.6 Mbp        0.0%    0.3%       8.0    MJSX01000001.1 Listeria monocytogenes strain BCW_2370 PROKKA_contig000001, whole genome shotgun sequence
# 2.5 Mbp        0.0%    0.2%       2.0    KK736542.1 Escherichia coli UCI 57 genomic scaffold aedZc-supercont1.1, whole genome shotgun sequence
# 2.5 Mbp        0.0%    0.2%       3.0    KQ088929.1 Escherichia coli strain CHS199 genomic scaffold aeviY-supercont1.1, whole genome shotgun sequence
# 2.5 Mbp        0.0%    0.2%       1.0    JNPT01000001.1 Escherichia coli 1-392-07_S1_C1 e139207S1C1.contig.0_1, whole genome shotgun sequence
# 2.5 Mbp        0.0%    0.2%       1.0    AFAV01000027.1 Escherichia coli PCN079 contig1, whole genome shotgun sequence
# 2.4 Mbp        0.0%    0.3%       1.0    JOZQ01000001.1 Listeria monocytogenes strain FSL R8-6637 contig000001, whole genome shotgun sequence
# 2.3 Mbp        0.0%    0.2%       4.0    HE572566.1 Escherichia coli HM605 draft genomic scaffold, scaffold00001, whole genome shotgun sequence
# 2.3 Mbp        0.0%    0.2%       1.0    LQVK01000001.1 Escherichia coli strain GN02820 GCID_ECOLID_00120_NODE_1.ctg_1, whole genome shotgun sequence
# 2.2 Mbp        0.0%    0.2%       1.0    JSOO01000002.1 Escherichia coli strain upec-112 upec-112_ctg_195, whole genome shotgun sequence
# 2.2 Mbp        0.0%    0.2%       2.0    KE701208.1 Escherichia coli KOEGE 10 (25a) genomic scaffold acYzM-supercont1.1, whole genome shotgun sequence
# 2.2 Mbp        0.0%    0.2%       1.0    JQFS01000001.1 Escherichia coli K1 strain 937 contig001, whole genome shotgun sequence
# 2.2 Mbp        0.0%    0.2%       2.0    JSSD01000002.1 Escherichia coli strain blood-08-1463 blood-08-1463_ctg_2927, whole genome shotgun sequence
# 2.2 Mbp        0.0%    0.2%       1.0    LM996590.1 Escherichia coli genome assembly FHI74, scaffold scaffold-1_contig-22.0_1_1216_[organism:Escherichia
# 2.2 Mbp        0.0%    0.2%       1.0    CDRQ01000001.1 Escherichia coli genome assembly 57A_A8_assembly, contig MI.57A_A8_c1, whole genome shotgun sequence
# 2.2 Mbp        0.0%    0.2%       1.0    AIBH01000001.1 Escherichia coli C639_08 C639_08_2, whole genome shotgun sequence
# 2.2 Mbp        0.0%    0.2%       2.0    JSYS01000001.1 Escherichia coli strain GSK2528 contig001, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.2%       1.0    AOQM01000001.1 Escherichia coli TOP550-2, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.2%       1.0    JHRG01000001.1 Escherichia coli 302053 302053_1, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.2%       3.0    JYEE01000001.1 Escherichia coli strain Chronic_salp contig_1, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.4%       1.0    JEVU01000001.1 Staphylococcus aureus W89268 adHvw-supercont1.1, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.4%       1.0    KK072083.1 Staphylococcus aureus M42514 genomic scaffold adHuc-supercont1.1, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.4%       2.0    CSET01000001.1 Staphylococcus aureus genome assembly 7907_4#14, scaffold ERS092867SCcontig000001, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.4%       1.0    CYPM01000001.1 Staphylococcus aureus genome assembly 12673_2#71, scaffold ERS410903SCcontig000001, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.2%       5.0    KB733107.1 Escherichia coli KTE84 genomic scaffold acATv-supercont1.1, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.4%       1.0    KK097848.1 Staphylococcus aureus DAR3565 genomic scaffold adLXq-supercont1.1, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.4%       1.0    FKUT01000001.1 Staphylococcus aureus strain MRSA genome assembly, contig: ERS411154SCcontig000001, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.2%       1.0    LYRV01000001.1 Escherichia coli strain IMT26197 contig_1, whole genome shotgun sequence
# 2.1 Mbp        0.0%    0.4%       1.0    FKMK01000001.1 Staphylococcus aureus strain NA genome assembly, contig: ERS365786SCcontig000001, whole genome shotgun sequence
# 2.0 Mbp        0.0%    0.4%       2.0    FMQD01000035.1 Staphylococcus aureus strain GKP136-54 genome assembly, contig: ERS189617SCcontig000035, whole genome shotgun sequence
# 2.0 Mbp        0.0%    0.2%       2.0    KI532355.1 Escherichia coli 908519 genomic scaffold Scaffold0, whole genome shotgun sequence
# 1.9 Mbp        0.0%    0.4%       2.0    KI912926.1 Enterococcus faecalis B284 genomic scaffold adiGp-supercont1.1, whole genome shotgun sequence
# 1.9 Mbp        0.0%    0.2%       1.0    CUEG01000001.1 Staphylococcus aureus genome assembly 5364_2#11, scaffold ERS010883.5364_2_11.1, whole genome shotgun sequence
# 1.8 Mbp        0.0%    0.4%       1.0    AYOL01000001.1 Enterococcus faecalis NY9 contig00001, whole genome shotgun sequence
# 1.8 Mbp        0.0%    0.4%       1.0    FFHQ01000001.1 TPA: Listeria monocytogenes strain 2842STDY5753961 genome assembly, contig: ERS409707SCcontig000001, whole genome shotgun se...
# 1.8 Mbp        0.0%    0.2%       1.0    CAFD01000138.1 Salmonella enterica subsp. salamae str. 3588/07, WGS project CAFD00000000 data, contig SAL1.0.600824, whole genome shotgun s...
# 1.8 Mbp        0.0%    0.4%       2.0    FGZR01000001.1 Staphylococcus aureus strain st2968 genome assembly, contig: ERS073697SCcontig000001, whole genome shotgun sequence
# 1.8 Mbp        0.0%    0.4%       1.0    FHCO01000001.1 Staphylococcus aureus strain st1913 genome assembly, contig: ERS073835SCcontig000001, whole genome shotgun sequence
# 1.7 Mbp        0.0%    0.2%       1.0    JNPK01000001.1 Escherichia coli 2-156-04_S3_C2 e215604S3C2.contig.0_1, whole genome shotgun sequence
# 1.7 Mbp        0.0%    0.4%       1.0    KB944692.1 Enterococcus faecalis SS-7 genomic scaffold acAqA-supercont1.1, whole genome shotgun sequence
# 1.7 Mbp        0.0%    0.4%       1.0    FHBW01000001.1 Staphylococcus aureus strain st1515 genome assembly, contig: ERS073808SCcontig000001, whole genome shotgun sequence
# 1.7 Mbp        0.0%    0.4%       1.0    MJBW01000001.1 Enterococcus faecalis strain Enfs94 Enfs94_S4_L001_R1_001_(paired)_contig_1, whole genome shotgun sequence
# 1.7 Mbp        0.0%    0.2%       1.0    MXLJ01001000.1 Salmonella enterica subsp. salamae serovar 30:1,z28:z6 strain 175-85 NODE_1000_length_1581_cov_1.07497, whole genome shotgun...
# 1.6 Mbp        0.0%    0.4%       2.0    FMOW01000030.1 Staphylococcus aureus strain MSSA777507T genome assembly, contig: ERS391193SCcontig000030, whole genome shotgun sequence
# 1.6 Mbp        0.0%    0.4%       1.0    LRNB01000001.1 Staphylococcus aureus subsp. aureus strain 685 contig_1, whole genome shotgun sequence
# 1.5 Mbp        0.0%    0.4%       1.0    AYND01000001.1 Enterococcus faecalis KS19 contig00001, whole genome shotgun sequence
# 1.5 Mbp        0.0%    0.4%       1.0    MLQA01000010.1 Staphylococcus aureus strain 2016 NODE_10_length_89244_cov_47.1279_ID_19, whole genome shotgun sequence
# 1.4 Mbp        0.0%    0.5%      19.0    CP011536.1 Lactobacillus fermentum 3872, complete genome
# 1.4 Mbp        0.0%    0.5%       1.0    CP019030.1 Lactobacillus fermentum strain SNUV175, complete genome
# 1.4 Mbp        0.0%    0.4%       1.0    LNTF01000001.1 Staphylococcus aureus strain SA7112 SA7112_Contig1, whole genome shotgun sequence
# 1.4 Mbp        0.0%    0.4%       1.0    JIXO01000001.1 Staphylococcus aureus MSSA-47 adZzG-supercont1.1, whole genome shotgun sequence
# 1.4 Mbp        0.0%    0.5%       1.0    LAIK01000012.1 Lactobacillus fermentum strain LfQi6 plasmid unnamed scaffold12, whole genome shotgun sequence
# 1.2 Mbp        0.0%    0.5%       1.0    AVAB01000110.1 Lactobacillus fermentum MTCC 8711 plasmid pLF01, complete sequence, whole genome shotgun sequence
# 1.1 Mbp        0.0%    0.6%       1.0    LBDH01000089.1 Lactobacillus fermentum strain 90 TC-4 Contig_102, whole genome shotgun sequence
# 1.1 Mbp        0.0%    0.2%       1.0    MXOX01000100.1 Salmonella enterica subsp. arizonae serovar 48:z4,z23,z32:- strain BCW_1532 NODE_100_length_14960_cov_0.108135, whole genome...
# 1.0 Mbp        0.0%    0.4%       1.0    MKNG01000001.1 Listeria monocytogenes strain NRRL B-33360 PROKKA_contig000001, whole genome shotgun sequence
# 1.0 Mbp        0.0%    0.4%       1.0    MKNX01000001.1 Listeria monocytogenes strain NRRL B-33381 PROKKA_contig000001, whole genome shotgun sequence
# 1.0 Mbp        0.0%    0.4%       1.0    MDQY01000001.1 Listeria monocytogenes strain A51 contig000001, whole genome shotgun sequence
# 0.8 Mbp        0.0%    0.2%       2.0    AVPM01000001.1 Bacillus sp. EGD-AK10 contig1, whole genome shotgun sequence
# 0.7 Mbp        0.0%    0.4%       1.0    MJSZ01000001.1 Listeria monocytogenes strain BCW_2372 PROKKA_contig000001, whole genome shotgun sequence
# 0.7 Mbp        0.0%    0.2%       1.0    MOEO01000010.1 Bacillus sp. FMQ74 NODE_10_length_89456_cov_116.443, whole genome shotgun sequence
# 0.6 Mbp        0.0%    0.2%       1.0    CCBU01000001.1 Escherichia coli genome assembly Ec80B_L1 MIRA assembly, contig Ec80B_L1_021214_c1, whole genome shotgun sequence
# 0.6 Mbp        0.0%    0.2%       1.0    BBVJ01000001.1 Escherichia albertii DNA, contig: NIAH_Bird_16contig0001, strain: NIAH_Bird_16, whole genome shotgun sequence
# 0.6 Mbp        0.0%    0.2%       2.0    CH991901.1 Escherichia albertii TW07627 scf_1111898267896 genomic scaffold, whole genome shotgun sequence
# 0.5 Mbp        0.0%    0.2%       1.0    BBVO01000001.1 Escherichia albertii DNA, contig: NIAH_Bird_26contig0001, strain: NIAH_Bird_26, whole genome shotgun sequence
# 380.0 kbp      0.0%    0.1%       1.0    KE007191.1 Salinibacillus aidingensis MSP4 genomic scaffold scaffold00001, whole genome shotgun sequence
# 330.0 kbp      0.0%    0.4%       1.0    CCEO01000001.1 Staphylococcus schweitzeri genome assembly Sa_FSA090, contig FSA090001, whole genome shotgun sequence
# 230.0 kbp      0.0%    0.1%       1.0    AABY01000832.1 Saccharomyces paradoxus NRRL Y-17217 contig_531, whole genome shotgun sequence
# 190.0 kbp      0.0%    0.2%       2.0    MBTW01000001.1 Citrobacter koseri strain B1B contig_1, whole genome shotgun sequence
# 180.0 kbp      0.0%    0.2%       2.0    MCOM01000001.1 Citrobacter freundii strain CF4_ST90 NODE_100_length_277_cov_22.235_ID_1587, whole genome shotgun sequence
# 130.0 kbp      0.0%    0.2%       1.0    CP011602.1 Kluyvera intermedia strain CAV1151, complete genome
# 110.0 kbp      0.0%    0.2%       3.0    FTLF01000037.1 Raoultella ornithinolytica strain Marseille-P1025 genome assembly, contig: contig00037, whole genome shotgun sequence
# 100.0 kbp      0.0%    0.2%       1.0    CP005991.1 Enterobacter sp. R4-368, complete genome
# 100.0 kbp      0.0%    0.2%       4.0    BCYS01000001.1 Kluyvera intermedia NBRC 102594 = ATCC 33110 DNA, contig: KIN01S_CON0001_0001, whole genome shotgun sequence
# 90.0 kbp       0.0%    0.2%       1.0    LNII01000001.1 Enterobacter sp. 50793107 contig1, whole genome shotgun sequence
# 90.0 kbp       0.0%    0.2%       1.0    FABP01000001.1 Peptoclostridium difficile isolate VL_0125 genome assembly, contig: NODE_1_length_79558_cov_49.1838_ID_1, whole genome shotg...
# 60.0 kbp       0.0%    0.2%       1.0    LWHA01000001.1 Pseudomonas punonensis strain D1-6 D1-6A_S5.NODE_1, whole genome shotgun sequence
# 60.0 kbp       0.0%    0.3%       2.0    JTAK01000001.1 Pseudomonas flexibilis strain JCM 14085 PT85_1, whole genome shotgun sequence
# 60.0 kbp       0.0%    0.2%       1.0    BCNN01000001.1 Lelliottia amnigena NBRC 105700 DNA, contig: ENTAMN1_CON00110, whole genome shotgun sequence
# 50.0 kbp       0.0%    0.2%       1.0    FAFB01000001.1 Peptoclostridium difficile isolate VL_0218 genome assembly, contig: NODE_1_length_35060_cov_8.78634_ID_1, whole genome shotg...
# 50.0 kbp       0.0%    0.5%       1.0    BCAH01000001.1 Lactobacillus gorillae DNA, contig:KZ01T01, strain: KZ01, whole genome shotgun sequence
# 50.0 kbp       0.0%    0.2%       1.0    LFEJ01000003.1 Cronobacter sp. DJ34 Bact_A5_1, whole genome shotgun sequence
# 50.0 kbp       0.0%    0.4%       1.0    JXCF01000004.1 Staphylococcus gallinarum strain DSM 20610 contig_4, whole genome shotgun sequence
# 50.0 kbp       0.0%    0.2%       1.0    CP009459.1 Cedecea neteri strain ND14a, complete genome
# 50.0 kbp       0.0%    0.2%       1.0    FWGT01000058.1 Klebsiella pneumoniae strain VRCO0021 genome assembly, contig: ERS682822SCcontig000058, whole genome shotgun sequence
# 40.0 kbp       0.0%    0.2%       1.0    MKEP01000001.1 Pseudomonas sp. J237 Contig_1, whole genome shotgun sequence
# 40.0 kbp       0.0%    0.1%       1.0    JGYI01000001.1 Pseudomonas veronii 1YB2 Contig1, whole genome shotgun sequence
# 40.0 kbp       0.0%    0.2%       1.0    AMWJ01000001.1 Pseudomonas putida CSV86 contig00002, whole genome shotgun sequence
# 40.0 kbp       0.0%    0.2%       1.0    LMLK01000001.1 Erwinia sp. Leaf53 contig_1, whole genome shotgun sequence
# 30.0 kbp       0.0%    0.3%       1.0    CP015857.1 Lactobacillus plantarum strain LZ227, complete genome
# 30.0 kbp       0.0%    0.2%       1.0    AP014655.1 Pseudomonas sp. MT-1 DNA, nearly complete genome
# 30.0 kbp       0.0%    0.2%       1.0    CVLT01000001.1 Providencia rettgeri genome assembly Providencia rettgeri H1736, contig contig00001, whole genome shotgun sequence
# 30.0 kbp       0.0%    0.2%       1.0    ATKM01000001.1 Pseudomonas sp. P818 P818-scaffold1, whole genome shotgun sequence
# 30.0 kbp       0.0%    0.2%       1.0    HE610985.1 Bacillus timonensis genomic scaffold, scaffold00001, whole genome shotgun sequence
# 20.0 kbp       0.0%    0.3%       1.0    KV804027.1 Enterococcus sp. HMSC063H10 genomic scaffold Scaffold1, whole genome shotgun sequence
# 20.0 kbp       0.0%    0.2%       4.0    CVTO01000163.1 Pseudomonas sp. 24 E 13 isolate 24 E 13 genome assembly, contig: contig000163, whole genome shotgun sequence
# 20.0 kbp       0.0%    0.8%       1.0    CP012328.1 Spiroplasma turonicum strain Tab4c, complete genome
# 20.0 kbp       0.0%    0.4%       1.0    FHGX01000001.1 Streptococcus pneumoniae strain 2245STDY5609077 genome assembly, contig: ERS356448SCcontig000001, whole genome shotgun sequence
# 20.0 kbp       0.0%    0.4%       1.0    LXRZ01000001.1 Staphylococcus epidermidis strain LRKNS092 LRKNS092_contig_1, whole genome shotgun sequence
# 20.0 kbp       0.0%    0.7%       1.0    AZGC01000044.1 Lactobacillus equigenerosi DSM 18793 = JCM 14505 strain DSM 18793 NODE_100, whole genome shotgun sequence
# 20.0 kbp       0.0%    0.2%       1.0    LKKK01000001.1 Pseudomonas sp. TTU2014-080ASC scaffold_0, whole genome shotgun sequence
# 20.0 kbp       0.0%    0.2%       1.0    JMQC01000001.1 Bacillus mycoides strain BHP DJ93.Contig27, whole genome shotgun sequence
# 20.0 kbp       0.0%    0.2%       1.0    LT629746.1 Pseudomonas lini strain BS3782 genome assembly, chromosome: I
# 10.0 kbp       0.0%    0.3%       1.0    LSMC01000001.1 Enterococcus mundtii strain QAUEM2808 Contig_10, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.4%       2.0    CH672400.1 marine gamma proteobacterium HTCC2207 scf_1099215844476 genomic scaffold, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.2%       1.0    MKVY01000001.1 Rhizobiales bacterium 65-9 SCNpilot_cont_300_bf_scaffold_10078, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.3%       1.0    LYUK01000001.1 Lactobacillus plantarum strain SRCM101060 contig00001, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.2%       1.0    JRWU01000008.1 Rouxiella chamberiensis strain 130333 contig1, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.3%       1.0    LT629736.1 Pseudomonas xinjiangensis strain NRRL B-51270 genome assembly, chromosome: I
# 10.0 kbp       0.0%    0.1%       1.0    CAPV01000001.1 Candida nivariensis CBS 9983 WGS project CAPV00000000 data, contig 123, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.2%       1.0    JH947125.1 Pedobacter arcticus A12 genomic scaffold Scaffold1, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.4%       1.0    AAJO01000553.1 Streptococcus agalactiae 18RS21 s_agalactiae_18rs21_2957, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.2%       1.0    FNJL01000093.1 Acidovorax cattleyae strain DSM 17101 genome assembly, contig: Ga0074841_193, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.2%       1.0    LXWE01000001.1 Alteromonadaceae bacterium XY-R5 NODE_1, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.3%       1.0    FOWY01000091.1 Halanaerobium congolense strain WG1 genome assembly, contig: Ga0073276_191, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.2%       1.0    AXZL01000001.1 Shewanella decolorationis S12 Contig1, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.3%       1.0    MEGL01000001.1 Rhodanobacter sp. SCN 68-63 ABT19_C0001, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.5%       1.0    LWMI01000001.1 Candidatus Arsenophonus triatominarum strain ATi contig000001, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.2%       1.0    AXDT01000001.1 Photorhabdus temperata J3 Contig001, whole genome shotgun sequence
# 10.0 kbp       0.0%    0.3%       1.0    JNLN01000001.1 Clostridium algidicarnis strain B3 BV55DRAFT_scf7180000000002_quiver_dupTrim_6430.1_C, whole genome shotgun sequence
# found less than 100 bp  in common. => exiting

# found 419 matches total;
# the recovered matches hit 89.0% of the abundance-weighted query.
# the recovered matches hit 32.2% of the query k-mers (unweighted).


conda activate 'centrifuge'
export CENTRIFUGE_OUTPUT='SRR11207265_1'

centrifuge -x "$DATA_DIR/p_compressed+h+v/p_compressed+h+v" \
    -U $FASTQ_FILE \
    --report-file "$DATA_DIR/centrifuge/${CENTRIFUGE_OUTPUT}_report.txt" \
    -S "$DATA_DIR/centrifuge/${CENTRIFUGE_OUTPUT}_results.txt"
# report file data/leishmania_revisited/centrifuge/SRR11207265_1_report.txt
# Number of iterations in EM algorithm: 53
# Probability diff. (P - P_prev) in the last iteration: 7.22319e-11
# Calculating abundance: 00:00:00



# ----------------------------------------------------------------------
# Task 5

# rerun sourmash and centrifuge using different files again
conda activate "rnaseq"
export FASTQ_FILENAME="SRR12596172_1.fastq.gz"
export FASTQ_FILE="$DATA_DIR/fastq/$FASTQ_FILENAME"

sourmash sketch dna -p scaled=10000,k=31,abund $FASTQ_FILE --name-from-first
mv "$FASTQ_FILENAME.sig" "$DATA_DIR/$FASTQ_FILENAME.sig"

sourmash gather -k 31 "$DATA_DIR/$FASTQ_FILENAME.sig" $GENOME_FILE --threshold-bp 100
# == This is sourmash version 4.6.1. ==
# == Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

# selecting specified query k=31
# loaded query: SRR12596172.1 1/1... (k=31, DNA)
# --
# loaded 93249 total signatures from 1 locations.
# after selecting signatures compatible with search, 93249 remain.

# Starting prefetch sweep across databases.
# Found 4474 signatures via prefetch; now doing gather.

# overlap     p_query p_match avg_abund
# ---------   ------- ------- ---------
# 30.0 kbp       4.7%   16.7%       5.0    AUNV01000052.1 Comamonadaceae bacterium JGI 0001003-E1...
# 30.0 kbp       0.9%    0.5%       1.5    LFRI01000001.1 Aquabacterium parvum strain B6 contig00...
# 20.0 kbp      15.8%    4.5%      25.0    AAMH01000001.1 Uncultured human fecal virus RVHLib2_2,...
# 20.0 kbp       1.3%    0.5%       2.0    FMDR01000062.1 uncultured Clostridium sp. isolate 2789...
# 30.0 kbp       3.2%    0.2%      10.0    AMRP01000001.1 Acidovorax citrulli ZJU1106 AASCcontig1...
# 30.0 kbp       0.3%    0.2%       1.0    BADL01000631.1 Ideonella sp. B508-1 DNA, contig: IB508...
# 20.0 kbp       0.3%    0.2%       1.0    CUEG01000001.1 Staphylococcus aureus genome assembly 5...
# 20.0 kbp       0.3%    0.2%       1.0    AELB01000405.1 Escherichia coli TW10722 TW10722_rep_c8...
# 20.0 kbp       0.9%    0.2%       3.0    AJKJ01000001.1 Citreicella sp. 357 C357_001, whole gen...
# 20.0 kbp       0.3%    0.2%       1.0    LSVE01000001.1 Pseudomonas stutzeri strain ODKF13 cont...
# 10.0 kbp       0.3%    0.2%       1.0    LFHV01000001.1 Escherichia coli strain E-D 371 EC1678_...
# 10.0 kbp       0.3%    0.2%       1.0    CP000749.1 Marinomonas sp. MWYL1, complete genome
# 10.0 kbp       6.3%    0.2%      20.0    JTBG01000002.1 Aeromonas caviae strain FDAARGOS_76 deg...
# 10.0 kbp       0.3%    0.2%       1.0    JYJT01000108.1 Vibrio parahaemolyticus strain 10-4255 ...
# 10.0 kbp       0.6%    0.2%       2.0    LZEY01000001.1 Morganella psychrotolerans strain GCSL-...
# 10.0 kbp       0.3%    0.2%       1.0    LGYY01000235.1 Shewanella sp. Sh95 contig_1, whole gen...
# 10.0 kbp       0.3%    0.2%       1.0    AIFR01000073.1 Escherichia coli DEC5B gecDEC5B.contig....
# 10.0 kbp       0.6%    0.3%       2.0    LT629780.1 Pseudomonas guangdongensis strain CCTCC 201...
# 10.0 kbp       0.6%    0.3%       2.0    GG703878.1 Prevotella copri DSM 18205 genomic scaffold...
# 10.0 kbp       0.3%    0.3%       1.0    MRCJ01000001.1 Neptunomonas phycophila strain 3CM2.5 C...
# 10.0 kbp       3.5%    0.3%      11.0    LN898212.1 Vitreoscilla sp. SN6 genome assembly Vitreo...
# found less than 100 bp  in common. => exiting

# found 21 matches total;
# the recovered matches hit 41.6% of the abundance-weighted query.
# the recovered matches hit 24.1% of the query k-mers (unweighted).


conda activate 'centrifuge'
export CENTRIFUGE_OUTPUT='SRR12596172_1'

centrifuge -x "$DATA_DIR/hvc/hvc" \
    -U $FASTQ_FILE \
    --report-file "$DATA_DIR/centrifuge/${CENTRIFUGE_OUTPUT}_report.txt" \
    -S "$DATA_DIR/centrifuge/${CENTRIFUGE_OUTPUT}_results.txt"
# report file data/leishmania_revisited/centrifuge/SRR12596172_1_report.txt
# Number of iterations in EM algorithm: 14
# Probability diff. (P - P_prev) in the last iteration: 0
# Calculating abundance: 00:00:00
