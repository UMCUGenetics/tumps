# tumps
tumor purity simulation

will only work on an HPC system. Can be easily adapted

## Basic usage
Required parameters are: 
- Tumor BAM file
- Normal BAM file
- Tumor Coverage
- Normal Coverage

Example basic usage:

```
./tumps.sh -t tumor.bam -n normal.bam -tc 90 -nc 30
```
Purities are specified providing the -p | --purities argument with a comma separated list. By default, purities are 0,10,20,25,50,75,100. 0 corresponds to the normal BAM file and 100 to the tumor BAM file.
Somatic SV calling is performed using GRIDSS by default. Other supported options are NanoSV or PBSV. SV calling can be disabled with --nosv.
I.e. for a run to get tumor purities of 30% and 60% and no somatic SV calling:

```
./tumps.sh -t tumor.bam -n normal.bam -tc 90 -nc 30 -p 30,60 --nosv
```

## Advanced usage

```
./tumps.sh --help

Parameters:
    -t|--tumor                    Path to tumor BAM file, must be indexed (Required)
    -n|--normal                    Path to normal BAM file, must be indexed (Required)
    -tc|--tumor_coverage                    Average tumor coverage (Required)
    -nc|--normal_coverage                    Average normal coverage for mixing (Required)
    -nsv|--normal_sv                    Path to normal BAM file to use for SV calling (if different than normal above), must be indexed [equal to --normal]
    -p|--purities                    Merge Tumor and Normal with the following tumor purities (comma separated): [0,10,20,25,50,75,100]
    -d|--output_dir                    OUTPUT DIRECTORY [.]
    --nosv                     Don't perform SV calling and stop after mixing the BAMs [FALSE]
    --nooverlap                     Don't perform truth and stop after SV calling [FALSE]
    -m|--mode                    SV calling mode, choose from: gridss, nanopore, pbsv [gridss]
    -g|--truthset                    Path to TRUTHSET VCF file for comparison [files/COLO829.somatic.vcf]
    --bed                    BED file for coverage calculations [files/human_hg19.bed]
    --sambamba                    Path to SAMBAMBA [/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba]
    --samtools                    Path to SAMTOOLS [/hpc/local/CentOS7/cog_bioinf/samtools-1.7/samtools]
    --mail                    Mail for HPC jobs [jespejov@umcutrecht.nl]
    --libgridss                    Path to LIBGRIDSS directory [scripts/libgridss/]
    --gridss_pon                    Path to GRIDSS PON directory (too big to share here) [/hpc/cog_bioinf/cuppen/project_data/Roel_pipeline_validation/gridss]
    --sniffles                    Path to Sniffles executable [/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/tools_kloosterman/Sniffles-1.0.8/bin/sniffles-core-1.0.8/sniffles]
    --survivor                    Path to SURVIVOR executable (for nanopore and pbsv) [/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/tools_kloosterman/SURVIVOR-1.0.6/Debug/SURVIVOR]
    --pbsv                    Path to PBSV executable [/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/bin/miniconda3/bin/pbsv]
    --ref                    Path to reference genome (for PBSV) [/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta]
    --nanosv_venv                    Path to NanoSV VENV folder [/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/bin/NanoSV/bin/NanoSV]
    --nanosv_config                    Path to NanoSV config file [files/config_COLO829_NGMLR.ini]

```
