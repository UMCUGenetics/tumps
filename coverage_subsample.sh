#! /bin/bash


#Jose Espejo Valle-Inclan
#2020
#jespejov@umcutrecht.nl/jaesvi@gmail.com

usage() {
echo "
Parameters:
    -b|--bam                    Path to BAM file for subsetting, must be indexed (Required)
    -bc|--bamcov                    Coverage of BAM file (Required)
    -n|--normal                    Path to normal BAM file for somatic SV calling, must be indexed (Required)
    -c|--coverages                    Subset BAM to following coverages (comma separated): [1,5,10,30,50]
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

"
exit
}

echo "0"
POSITIONAL=()

##DEFAULT ARGUMENTS

OUTDIR=.
MODE=gridss
TRUTHSET=files/COLO829.somatic.vcf
COVERAGES="1,5,10,30,50"
BED=files/human_hg19.bed
SAMBAMBA=/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba
SAMTOOLS=/hpc/local/CentOS7/cog_bioinf/samtools-1.7/samtools
MAIL=jespejov@umcutrecht.nl
LIBGRIDSS=script/libgridss/
GRIDSS_PON=/hpc/cog_bioinf/cuppen/project_data/Roel_pipeline_validation/gridss
SNIFFLES=/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/tools_kloosterman/Sniffles-1.0.8/bin/sniffles-core-1.0.8/sniffles
SURVIVOR=/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/tools_kloosterman/SURVIVOR-1.0.6/Debug/SURVIVOR
NANOSV_VENV=/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/bin/NanoSV/
NANOSV_CONFIG=files/config_COLO829_NGMLR.ini
PBSV=/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/bin/miniconda3/bin/pbsv
REF=/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta
echo "00"
##READ PARAMETERS
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
    -h|--help)
    usage
    shift
    shift# past argument
    ;;
    -b|--bam)
    BAM="$2"
    shift # past argument
    shift # past value
    ;;
    -bc|--bamcov)
    BAMCOV="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--normal)
    NORMAL="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--output_dir)
    OUTDIR="$2"
    shift # past argument
    shift # past value
    ;;
    --nosv)
    NOSV_FLAG=true
    shift # past argument
    shift # past value
    ;;
    --nooverlap)
    NOOLAP_FLAG=true
    shift # past argument
    shift # past value
    ;;
    -m|--mode)
    MODE="$2"
    shift # past argument
    shift # past value
    ;;
    -g|--truthset)
    TRUTHSET="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--COVERAGES)
    COVERAGES="$2"
    shift # past argument
    shift # past value
    ;;
    --bed)
    BED="$2"
    shift # past argument
    shift # past value
    ;;
    --sambamba)
    SAMBAMBA="$2"
    shift # past argument
    shift # past value
    ;;
    --samtools)
    SAMTOOLS="$2"
    shift # past argument
    shift # past value
    ;;
    --mail)
    MAIL="$2"
    shift # past argument
    shift # past value
    ;;
     --libgridss)
    LIBGRIDSS="$2"
    shift # past argument
    shift # past value
    ;;
     --gridss_pon)
    GRIDSS_PON="$2"
    shift # past argument
    shift # past value
    ;;
     --sniffles)
    SNIFFLES="$2"
    shift # past argument
    shift # past value
    ;;
    --survivor)
    SURVIVOR="$2"
    shift # past argument
    shift # past value
    ;;
    --pbsv)
    PBSV="$2"
    shift # past argument
    shift # past value
    ;;
    --ref)
    REF="$2"
    shift # past argument
    shift # past value
    ;;
     --nanosv_venv)
    NANOSV_VENV="$2"
    shift # past argument
    shift # past value
    ;;
     --nanosv_config)
    NANOSV_CONFIG="$2"
    shift # past argument
    shift # past value
    ;;
    esac
done
echo "HERE"
set -- "${POSITIONAL[@]}" # restore positional parameters
echo "2"
#Check required parameters
if [ -z $BAM ]; then
  echo "Missing -b|--bam parameter"
  usage
elif [ -z $BAMCOV ]; then
  echo "Missing -bc|--bamcov parameter"
  usage
elif [ -z $NORMAL ]; then
  echo "Missing -n|--normal parameter"
  usage
fi
echo "3"
SCRIPT_DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
GRIDSS_SCRIPT=${SCRIPT_DIR}/scripts/GRIDSS.sh
OVERLAP_SCRIPT=${SCRIPT_DIR}/scripts/annotate_sv_vcf_file_without_ori.py

RAND=$(cat /dev/urandom | tr -cd 'a-zA-Z0-9' | head -c 7)
echo "Random code for this run is $RAND"

OUTDIR=`realpath $OUTDIR`
if [ ! -d $OUTDIR ]; then
    echo "Output directory ($OUTDIR) does not exist, creating"
    mkdir $OUTDIR
fi

echo "Moving to $OUTDIR"
cd $OUTDIR

if [ ! -d logs ]; then
    echo "Logs directory (${OUTDIR}/logs) does not exist, creating"
    mkdir logs
fi

if [ ! -f ${BAM}.bai ] || [ ! -f ${NORMAL}.bai ] ; then
  echo "BAM files must be indexed!"
  exit
fi

COVS=`echo "${COVERAGES}" | sed "s:,: :g"`

for COV in $COV
do

  COVBAM=cov${COV}.bam
  echo ""
  echo "###COVERAGE ${COV}###"
  #Calculate percentage
  SUBSAMPLE_TUMOR=$(awk -v coverage="$COV" -v read_depth="$BAMCOV" 'BEGIN{print coverage/read_depth}')
  echo "Subsampling tumor by $SUBSAMPLE_TUMOR"
  echo "$SAMBAMBA view -f bam -s $SUBSAMPLE_TUMOR -o $COVBAM $BAM" |\
  qsub -cwd -l h_rt=12:0:0 -l h_vmem=30G -N subsample_${COV}_${RAND} -e logs/subsample_${COV}_${RAND}.err -o logs/subsample_${COV}_${RAND}.o
    if [ $NOSV_FLAG ]; then
      echo "Skipping SV calling"
      continue; fi
    
  if [ $MODE == "gridss" ]; then
    RAWVCF=${COVBAM/.bam/.raw.vcf}
    SOMVCF=${RAWVCF/.raw.vcf/.somatic.vcf}
    SOMFULLVCF=${RAWVCF/.raw.vcf/.somatic.full.vcf}
    
    if  [ ! -f $RAWVCF ]; then
    echo "GRIDSS calling"
    qsub -cwd -hold_jid subsample_${COV}_${RAND} -l h_rt=48:0:0 -l h_vmem=70G -pe threaded 4 -N GRIDSS_calling_${COV}_${RAND} -e logs/GRIDSS_calling_${COV}_${RAND}.err -o logs/GRIDSS_calling_${COV}_${RAND}.o $GRIDSS_SCRIPT -t $COVBAM -n $NORMAL_SV -o $RAWVCF
    else echo "$RAWVCF exists, skipping to SV filtering"
    fi
    
    if  [ ! -f $SOMVCF ]; then
    echo "GRIDSS filter"
    echo "#! /bin/sh" > logs/GRIDSS_filter_${RAND}.sh
    echo "guixr load-profile /hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/bin/hmf-pipeline -- <<EOF" >> logs/GRIDSS_filter_${RAND}.sh
    echo "Rscript ${LIBGRIDSS}/gridss_somatic_filter.R --normalordinal 2 -p $GRIDSS_PON -i ${RAWVCF} -o ${SOMVCF} -f ${SOMFULLVCF} -s $LIBGRIDSS" >> logs/GRIDSS_filter_${RAND}.sh
    echo "mv ${SOMVCF}.bgz ${SOMVCF}.gz" >> logs/GRIDSS_filter_${RAND}.sh
    echo "mv ${SOMFULLVCF}.bgz ${SOMFULLVCF}.gz" >> logs/GRIDSS_filter_${RAND}.sh
    echo "gunzip ${SOMVCF}.gz" >> logs/GRIDSS_filter_${RAND}.sh
    echo "gunzip ${SOMFULLVCF}.gz" >> logs/GRIDSS_filter_${RAND}.sh
    echo "EOF" >> logs/GRIDSS_filter_${RAND}.sh
    qsub -cwd -l h_rt=12:0:0 -l h_vmem=20G -hold_jid GRIDSS_calling_${COV}_${RAND} -N GRIDSS_filter_${COV}_${RAND} -e logs/GRIDSS_filter_${COV}_${RAND}.err -o logs/GRIDSS_filter_${COV}_${RAND}.o logs/GRIDSS_filter_${RAND}.sh
    else echo "$SOMVCF exists, skipping to cleanup"
    fi
    
    if [ ! -f logs/GRIDSS_cleanup_${COV}.done ]; then
    echo "Cleanup GRIDSS"
    echo "if [ -f ${SOMVCF} ]; then rm -r ${SOMVCF}.bgz.tbi ${SOMFULLVCF}.bgz.tbi ${RAWVCF}.idx purity${COV}.gridss.assembly* GRIDSS_wdir_purity${COV}; touch logs/GRIDSS_cleanup_${COV}.done; fi" | qsub -cwd -hold_jid GRIDSS_filter_${COV}_${RAND} -l h_rt=0:10:0 -l h_vmem=1G -e logs/GRIDSS_cleanup_${COV}_${RAND}.err -o logs/GRIDSS_cleanup_${COV}_${RAND}.o -N GRIDSS_cleanup_${COV}_${RAND}
    else echo "Directory looks clean, skipping to overlap"
    fi
    
    if [ $NOOLAP_FLAG ]; then
      echo "Skipping overlap with truthset"
      continue; fi
        
    if [ ! -f logs/truth_overlap_${COV}.done ]; then
    echo "Truthset overlap"
    RAWTRUTH=${RAWVCF/.vcf/.truth.vcf}
    SOMTRUTH=${SOMVCF/.vcf/.truth.vcf}
    SOMFULLTRUTH=${SOMFULLVCF/.vcf/.truth.vcf}
    TRUTHPUR=truth.purity${COV}.vcf

    echo "python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $RAWVCF --file2 $TRUTHSET > $RAWTRUTH; \
    python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $SOMVCF --file2 $TRUTHSET > $SOMTRUTH; \
    python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $SOMFULLVCF --file2 $TRUTHSET > $SOMFULLTRUTH; \
    python $OVERLAP_SCRIPT --distance 100 --annotation ${COV}_RAW --input $TRUTHSET --file2 $RAWVCF | python $OVERLAP_SCRIPT --distance 100 --annotation ${COV}_SOM --input /dev/stdin --file2 $SOMVCF | \
    python $OVERLAP_SCRIPT --distance 100 --annotation ${COV}_FULL --input /dev/stdin --file2 $SOMFULLVCF > $TRUTHPUR;
    touch logs/truth_overlap_${COV}.done" | \
    qsub -cwd -hold_jid GRIDSS_cleanup_${COV}_${RAND} -N truthset_overlap_${COV}_${RAND} -l h_rt=2:0:0 -e logs/truthset_overlap_${COV}_${RAND}.err -o logs/truthset_overlap_${COV}_${RAND}.o
    else echo "Overlap files exist, skipping"
    fi

  elif [ $MODE == "nanopore" ]; then
    NSV_TUMVCF=${COVBAM/.bam/.tumor.raw.nanosv.vcf}
    SNF_TUMVCF=${COVBAM/.bam/.tumor.raw.sniffles.vcf}
    MERGEDVCF=${COVBAM/.bam/.merged.vcf}
    SOMVCF=${COVBAM/.bam/.somatic.vcf}
    NORSVNAME=`basename $NORMAL_SV`
    NSV_NORVCF=${NORSVNAME/.bam/.nanosv.vcf}
    SNF_NORVCF=${NORSVNAME/.bam/.sniffles.vcf}
    
    if [ ! -f logs/normal.run ]; then
        echo "$SNIFFLES -m $NORMAL_SV -v $SNF_NORVCF -t 4 --report_BND --genotype" | \
        qsub -cwd -pe threaded 4 -l h_rt=24:0:0 -l h_vmem=20G -N SNIFFLES_NORMAL_${RAND} -e logs/SNIFFLES_NORMAL_${RAND}.err -o logs/SNIFFLES_NORMAL_${RAND}.o
        echo ". $NANOSV_VENV/bin/activate; $NANOSV_VENV/bin/NanoSV -t 4 -c $NANOSV_CONFIG -b $NANOSV_VENV/bin/human_hg19.bed -o $NSV_NORVCF $NORMAL_SV" | \
        qsub -cwd -pe threaded 4 -l h_rt=96:0:0 -l h_vmem=80G -N NANOSV_NORMAL_${RAND} -e logs/NANOSV_NORMAL_${RAND}.err -o logs/NANOSV_NORMAL_${RAND}.o
        touch logs/normal.run
    else echo "Normal VCF exists, or the job is running. Skipping Normal SV calling"
    fi
    
    if [ ! -f $SNF_TUMVCF ]; then
        echo "$SNIFFLES -m $COVBAM -v $SNF_TUMVCF -t 4 --report_BND --genotype" | \
        qsub -hold_jid subsample_${COV}_${RAND} -cwd -pe threaded 4 -l h_rt=24:0:0 -l h_vmem=20G -N SNIFFLES_${COV}_${RAND} -e logs/SNIFFLES_${COV}_${RAND}.err -o logs/SNIFFLES_${COV}_${RAND}.o
    else echo "$SNF_TUMVCF exists, skipping SV calling"
    fi
    
    if [ ! -f $NSV_TUMVCF ]; then
        echo ". $NANOSV_VENV/bin/activate; $NANOSV_VENV/bin/NanoSV -t 4 -c $NANOSV_CONFIG -b $NANOSV_VENV/bin/human_hg19.bed -o $NSV_TUMVCF $COVBAM" | \
        qsub -hold_jid subsample_${COV}_${RAND} -cwd -pe threaded 4 -l h_rt=96:0:0 -l h_vmem=80G -N NANOSV_${COV}_${RAND} -e logs/NANOSV_${COV}_${RAND}.err -o logs/NANOSV_${COV}_${RAND}.o
    else echo "$NSV_TUMVCF exists, skipping SV calling"
    fi
    
    
    if [ ! -f $SOMVCF ]; then
        echo "$NSV_NORVCF" > purity${COV}.nanopore.list
        echo "$SNF_NORVCF" >> purity${COV}.nanopore.list
        echo "$NSV_TUMVCF" >> purity${COV}.nanopore.list
        echo "$SNF_TUMVCF" >> purity${COV}.nanopore.list
        echo "$SURVIVOR merge purity${COV}.nanopore.list 200 1 0 0 0 0 $MERGEDVCF; grep -P '^#|VEC=00[01][01];' $MERGEDVCF > ${SOMVCF}" | \
        qsub -hold_jid SNIFFLES_${COV}_${RAND} -hold_jid SNIFFLES_NORMAL_${RAND} -hold_jid NANOSV_${COV}_${RAND} -hold_jid NANOSV_NORMAL_${RAND} -N NANOPORE_filter_${COV}_${RAND} -cwd -l h_rt=2:0:0 -e logs/NANOPORE_filter_${COV}_${RAND}.err -o logs/NANOPORE_filter_${COV}_${RAND}.o
    else echo "$SOMVCF exists, skipping SV filtering"
    fi
    
    if [ $NOOLAP_FLAG ]; then
      echo "Skipping overlap with truthset"
      continue; fi
    
    if [ ! -f logs/truth_overlap_${COV}.done ]; then
    RAWTRUTH=${MERGEDVCF/.vcf/.truth.vcf}
    SOMTRUTH=${SOMVCF/.vcf/.truth.vcf}
    TRUTHPUR=truth.purity${COV}.vcf
    
    echo "python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $MERGEDVCF --file2 $TRUTHSET > $RAWTRUTH; \
    python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $SOMVCF --file2 $TRUTHSET > $SOMTRUTH; \
    python $OVERLAP_SCRIPT --distance 100 --annotation ${COV}_RAW --input $TRUTHSET --file2 $MERGEDVCF | python $OVERLAP_SCRIPT --distance 100 --annotation ${COV}_SOM --input /dev/stdin --file2 $SOMVCF > $TRUTHPUR; \
    touch logs/truth_overlap_${COV}.done" | \
    qsub -cwd -hold_jid  NANOPORE_filter_${COV}_${RAND} -N truthset_overlap_${COV}_${RAND} -l h_rt=2:0:0 -e logs/truthset_overlap_${COV}_${RAND}.err -o logs/truthset_overlap_${COV}_${RAND}.o
    else echo "Overlap files exist, skipping"
    fi


  elif [ $MODE == "pbsv" ]; then
    GET_SOMATIC=${SCRIPT_DIR}/get_somatic_pacbio.py
    TUMSIG=${COVBAM/.bam/.tumor.svsig.gz}
    NORNAME=`basename $NORMAL_SV`
    NORSIG=${NORNAME/.bam/.normal.svsig.gz}
    MERGVCF=${COVBAM/.bam/.merged.vcf}
    SOMVCF=${COVBAM/.bam/.somatic.vcf}
    if [ ! -f logs/pbsv_normal.run ]; then
        echo "$PBSV discover $NORMAL_SV $NORSIG" | \
        qsub -cwd -l h_rt=24:0:0 -l h_vmem=60G -N PBSV_DISC_NORMAL_${RAND} -e logs/PBSV_DISC_NORMAL_${RAND}.err -o logs/PBSV_DISC_NORMAL_${RAND}.o
        touch logs/pbsv_normal.run
    else echo"$NORSIG exists, or the job is running. Skipping Normal SV discovery"
    fi
    
    if [ ! -f $TUMSIG ]; then
        echo "$PBSV discover $COVBAM $TUMSIG" | \
        qsub -hold_jid subsample_${COV}_${RAND} -cwd -l h_rt=24:0:0 -l h_vmem=60G -N PBSV_DISC_${COV}_${RAND} -e logs/PBSV_DISC_${COV}_${RAND}.err -o logs/PBSV_DISC_${COV}_${RAND}.o
    else echo "$TUMSIG exists, skipping SV discovery"
    fi
    
    if [ ! -f $MERGVCF ]; then
        echo "$PBSV call -j 4 $REF $NORSIG $TUMSIG $MERGVCF" | \
        qsub -cwd -l h_rt=24:0:0 -l h_vmem=60G -pe threaded 4 -hold_jid PBSV_DISC_${COV}_${RAND} -hold_jid PBSV_DISC_NORMAL_${RAND} -N PBSV_CALL_${COV}_${RAND} -e logs/PBSV_CALL_${COV}_${RAND}.err -o logs/PBSV_CALL_${COV}_${RAND}.o
    else echo "$MERGVCF exists, skipping SV calling"
    fi
    
    if [ ! -f $SOMVCF ]; then
        echo "python $GET_SOMATIC $MERGVCF | awk -F $'\t' '{OFS = FS}{\$10=\"\"; print \$0}' > $SOMVCF" | \
        qsub -l h_rt=1:0:0 -cwd -hold_jid PBSV_CALL_${COV}_${RAND} -N PBSV_FILTER_${COV}_${RAND} -e logs/PBSV_FILTER_${COV}_${RAND}.err -o logs/PBSV_FILTER_${COV}_${RAND}.o
    else echo "$SOMVCF exists, skipping SV filtering"
    fi
    
    if [ $NOOLAP_FLAG ]; then
      echo "Skipping overlap with truthset"
        continue; fi
    
    if [ ! -f logs/truth_overlap_${COV}.done ]; then
    RAWTRUTH=${MERGVCF/.vcf/.truth.vcf}
    SOMTRUTH=${SOMVCF/.vcf/.truth.vcf}
    TRUTHPUR=truth.purity${COV}.vcf
    
    echo "python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $MERGVCF --file2 $TRUTHSET > $RAWTRUTH; \
    python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $SOMVCF --file2 $TRUTHSET > $SOMTRUTH; \
    python $OVERLAP_SCRIPT --distance 100 --annotation ${COV}_RAW --input $TRUTHSET --file2 $MERGVCF | python $OVERLAP_SCRIPT --distance 100 --annotation ${COV}_SOM --input /dev/stdin --file2 $SOMVCF > $TRUTHPUR; \
    touch logs/truth_overlap_${COV}.done" | \
    qsub -cwd -hold_jid PBSV_FILTER_${COV}_${RAND} -N truthset_overlap_${COV}_${RAND} -l h_rt=2:0:0 -e logs/truthset_overlap_${COV}_${RAND}.err -o logs/truthset_overlap_${COV}_${RAND}.o
    else echo "Overlap files exist, skipping"
    fi
  fi
done

done