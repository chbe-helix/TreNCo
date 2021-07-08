# Topologically Associated Domain (TAD) aware Regulatory Network Construction (TReNCo)

## About

TReNCo is designed to be an all-in-one Transcription Factor gene-by-gene regulatory network builder. TReNCo utilizes common assay’s, ChIP-seq, RNA-seq, and TAD boundaries as a hard cutoff, instead of distance based, to efficiently create context-specific TF-gene regulatory networks. There are three parts to the network construction: i) Enhancer-gene mapping ii) TF-enhancer mapping and iii) TAD Enhancer-gene weighing.

__This is a fork of the original TReNCo found at https://git.biohpc.swmed.edu/BICF/Software/trenco for which you need an account__

## Release

May 1, 2021 - v1.0.0 - Initial release of core TReNCo code

## Installation

### Requirements
python    3.7  <
conda     4.7  <
bcftools  1.4  <
samtools  1.4  <
bedops    2.4  <
bedtools  2.29 <
meme      5.0  <

### Basic Steps on Linux

Clone TReNCo repo into local repository and add to path
(Assuming .bashrc file and repo location in ~/ home directory on Linux)
```
git clone https://git.biohpc.swmed.edu/BICF/Software/trenco.git
echo "export PATH=~/trenco:$PATH" >> ~/.bashrc
echo "export PYHONPATH=~/trenco_modules:$PYTHONPATH" >> ~/.bashrc
```

Activate TReNCo conda environment
```
conda env create -f trenco_env.yml
conda activate trenco_env
```

## Usage
```
trenco --design [DESIGN FILE (txt)] --alignment [ALIGNMENT FILES (tsv)] --expression [EXPRESSION FILES (bam)] --peaks [PEAK FILES (bed)] -g [GENOME VERSION] [OPTIONS]
```
TAD aware Regulatory Network Construction (TReNCo)

### Main Arguments
```
  -h, --help                              |show this help message and exit. 
------------------------------------------|----------------------------------------------------------
                                          |
  --design DESIGN                         |Design file containing link information to samples.   
------------------------------------------|----------------------------------------------------------
                                          |
  --alignment ALIGNMENT [ALIGNMENT ...]   |Full path to ChIP alingment files in bam format
------------------------------------------|----------------------------------------------------------
                                          |
  --expression EXPRESSION [EXPRESSION ...]|Full path to transcript expression table in tsv format
------------------------------------------|----------------------------------------------------------
                                          |
  --enhMarks TARGET                       |Mark for enchancers: Default H3K27ac
------------------------------------------|----------------------------------------------------------
                                          |
  --tadBED TADBED                         |TAD file: Default - mm10 TAD in common_files
------------------------------------------|----------------------------------------------------------
                                          |
  -p, --peaks PEAKS [PEAKS ...].          |The full path to peak files in bed format
------------------------------------------|----------------------------------------------------------
                                          |
  -g, --genome-version GVERS.             |Version of genome to use. Default is newest
                                          |
```
### Optional Arguments:
```
  -r, --region REGION                     |The number of bases pairs to exclude around TSS 
                                          |  (Default: 2500)
------------------------------------------|----------------------------------------------------------
                                          |
  -q, --promoter-range PROMOTER_RANGE     |Range of numbers before TSS and after TSS to consider as Promoter 
                                          |  (Default: 1000-200)
------------------------------------------|----------------------------------------------------------
                                          |
  -d, --distance DISTANCE                 |The number of bases pairs to merge between adjacent peaks 
                                          |  (Default: 150)
------------------------------------------|----------------------------------------------------------
                                          |
  --annotation-file ANNOTFNAME            |Genode annotations file in gtf format 
                                          |  (overwrites --annotation-version and --organism)
------------------------------------------|----------------------------------------------------------
                                          |
  --annotation-version ANNOTATIONS        |The Gencode annotations file in gtf format. (Default: vM4) 
                                          |  (WARNING: All entries are indexed to this version)
------------------------------------------|----------------------------------------------------------
                                          |
  --organism ORGANISM                     |Organism gencode to download (Default: mouse)
------------------------------------------|----------------------------------------------------------
                                          |
  -b, --biotypes BIOTYPES [BIOTYPES ...]  |The biotypes to get transcript TSS. (default: protein)
------------------------------------------|----------------------------------------------------------
                                          |
  --meme-db MEME_DB                       |MEME database to use (Default: cis-bp)
------------------------------------------|----------------------------------------------------------
                                          |
  --db DB                                 |Motif database name if different from species 
                                          |  (ex JASPER CORE 2014 for jasper)
------------------------------------------|----------------------------------------------------------
                                          |
  --bed BED                               |ChIP and Promoter bed file for getting motifs 
                                          |  (ex enh.bed,promoter.bed)
```
