#! /bin/bash

#Jose Espejo Valle-Inclan
#2020
#jespejov@umcutrecht.nl/jaesvi@gmail.com

usage() {
echo "
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

"
exit
}

POSITIONAL=()

  
SCRIPT_DIR=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
GRIDSS_SCRIPT=${SCRIPT_DIR}/scripts/GRIDSS.sh
OVERLAP_SCRIPT=${SCRIPT_DIR}/scripts/annotate_sv_vcf_file_without_ori.py

##DEFAULT ARGUMENTS
TUMOR_COVERAGE=0
NORMAL_COVERAGE=0
OUTDIR=.
MODE=gridss
TRUTHSET=${SCRIPT_DIR}/files/truthset_somaticSVs_COLO829.vcf
PURITIES="0,10,20,25,50,100"
BED=${SCRIPT_DIR}/files/human_hg19.bed
SAMBAMBA=/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba
SAMTOOLS=/hpc/local/CentOS7/cog_bioinf/samtools-1.7/samtools
MAIL=jespejov@umcutrecht.nl
LIBGRIDSS=${SCRIPT_DIR}/script/libgridss/
GRIDSS_PON=/hpc/cog_bioinf/cuppen/project_data/Roel_pipeline_validation/gridss
SNIFFLES=/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/tools_kloosterman/Sniffles-1.0.8/bin/sniffles-core-1.0.8/sniffles
SURVIVOR=/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/tools_kloosterman/SURVIVOR-1.0.6/Debug/SURVIVOR
NANOSV_VENV=/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/bin/NanoSV/
NANOSV_CONFIG=${SCRIPT_DIR}/files/config_COLO829_NGMLR.ini
PBSV=/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/bin/miniconda3/bin/pbsv
REF=/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta


##READ PARAMETERS
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
    -h|--help)
    usage
    shift # past argument
    ;;
    -t|--tumor)
    TUMOR="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--normal)
    NORMAL="$2"
    shift # past argument
    shift # past value
    ;;
    -tc|--tumor_coverage)
    TUMOR_COVERAGE="$2"
    shift # past argument
    shift # past value
    ;;
    -nc|--normal_coverage)
    NORMAL_COVERAGE="$2"
    shift # past argument
    shift # past value
    ;;
    -nsv|--normal_sv)
    NORMAL_SV="$2"
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
    -p|--purities)
    PURITIES="$2"
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
set -- "${POSITIONAL[@]}" # restore positional parameters

#Check required parameters
if [ -z $TUMOR ]; then
  echo "Missing -t|--tumor parameter"
  usage
elif [ -z $NORMAL ]; then
  echo "Missing -n|--normal parameter"
  usage
elif [ -z $TUMOR_COVERAGE ]; then
  echo "Missing -tc|--tumor_coverage parameter"
  usage
elif [ -z $NORMAL_COVERAGE ]; then
  echo "Missing -nc|--normal_coverage parameter"
  usage
elif [ -z $NORMAL_SV ]; then
  echo "Using $NORMAL also as normal for SV"
  NORMAL_SV=$NORMAL
fi


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

if [ ! -f ${TUMOR}.bai ] || [ ! -f ${NORMAL}.bai ] ; then
  echo "BAM files must be indexed!"
  exit
fi

PURS=`echo $PURITIES | sed "s:,: :g"`

echo "Tumor percentage purities: $PURS"
echo "Tumor coverage: $TUMOR_COVERAGE"
echo "Normal coverage: $NORMAL_COVERAGE"


for PURITY in $PURS

do
    PURBAM=purity${PURITY}.bam
    echo ""
    echo "###PURITY ${PURITY}###"
    
    
    ###SET THE CONTROLS
    if [ $PURITY -eq "100" ]; then
        ln -sf $TUMOR purity100.bam
    elif [ $PURITY -eq "0" ]; then
        ln -sf $NORMAL purity0.bam
    else
    ###CREATE THE MIXED BAM FILES
        if  [ ! -f logs/${PURBAM}.done ]; then
            echo "Calculating ${PURITY}% purity"
            RATIO=`echo "${PURITY} / 100" | bc -l`
            TCOV=`echo "${TUMOR_COVERAGE} * ${RATIO}" | bc -l`
            echo "Subsample tumor BAM to a ${PURITY}%. That is ${TCOV}x coverage"
            echo "$SAMBAMBA view -f bam -s $RATIO -o tumorOnly${PURITY}.bam.tmp $TUMOR" | \
            qsub -cwd -l h_rt=12:0:0 -l h_vmem=30G -N tumor_subsample_${PURITY}_${RAND} -e logs/tumor_subsample_${PURITY}_${RAND}.err -o logs/tumor_subsample_${PURITY}_${RAND}.o
            
            NCOV=`echo "${TUMOR_COVERAGE} - ${TCOV}" | bc -l`
            NPUR=`echo "${NCOV} / ${NORMAL_COVERAGE}" | bc -l`
            echo "Subsample ${NCOV}x coverage from normal BAM. That is $NPUR of the original"
            echo "$SAMBAMBA view -f bam -s $NPUR -o normalOnly${PURITY}.bam.tmp $NORMAL" | \
            qsub -cwd -l h_rt=12:0:0 -l h_vmem=30G -N normal_subsample_${PURITY}_${RAND} -e logs/normal_subsample_${PURITY}_${RAND}.err -o logs/normal_subsample_${PURITY}_${RAND}.o
            
            echo "Merge both"
            echo "$SAMBAMBA merge $PURBAM tumorOnly${PURITY}.bam.tmp normalOnly${PURITY}.bam.tmp" | \
            qsub -cwd -hold_jid tumor_subsample_${PURITY}_${RAND} -hold_jid normal_subsample_${PURITY}_${RAND} -l h_rt=12:0:0 -l h_vmem=30G -N merge_${PURITY}_${RAND} -e logs/merge_${PURITY}_${RAND}.err -o logs/merge_${PURITY}_${RAND}.o
            
            echo "Clean up for purity $PURITY"
            echo "if [ -f purity$PURITY.bam ]; then rm tumorOnly${PURITY}.bam.tmp normalOnly${PURITY}.bam.tmp; touch logs/${PURBAM}.done; fi" | \
            qsub -cwd -hold_jid merge_${PURITY}_${RAND} -N cleanup_${PURITY}_${RAND} -e logs/cleanup_${PURITY}_${RAND}.err -o logs/cleanup_${PURITY}_${RAND}.o
        else echo "${PURBAM}.done exists, skipping to SV calling"
        fi
    fi
    if [ $NOSV_FLAG ]; then
      echo "Skipping SV calling"
      continue; fi
    
    if [ $MODE == "gridss" ]; then
        RAWVCF=${PURBAM/.bam/.raw.vcf}
        SOMVCF=${RAWVCF/.raw.vcf/.somatic.vcf}
        SOMFULLVCF=${RAWVCF/.raw.vcf/.somatic.full.vcf}
        
        if  [ ! -f $RAWVCF ]; then
        echo "GRIDSS calling"
        qsub -cwd -hold_jid merge_${PURITY}_${RAND} -l h_rt=48:0:0 -l h_vmem=70G -pe threaded 4 -N GRIDSS_calling_${PURITY}_${RAND} -e logs/GRIDSS_calling_${PURITY}_${RAND}.err -o logs/GRIDSS_calling_${PURITY}_${RAND}.o $GRIDSS_SCRIPT -t $PURBAM -n $NORMAL_SV -o $RAWVCF
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
        qsub -cwd -l h_rt=12:0:0 -l h_vmem=20G -hold_jid GRIDSS_calling_${PURITY}_${RAND} -N GRIDSS_filter_${PURITY}_${RAND} -e logs/GRIDSS_filter_${PURITY}_${RAND}.err -o logs/GRIDSS_filter_${PURITY}_${RAND}.o logs/GRIDSS_filter_${RAND}.sh
        else echo "$SOMVCF exists, skipping to cleanup"
        fi
        
        if [ ! -f logs/GRIDSS_cleanup_${PURITY}.done ]; then
        echo "Cleanup GRIDSS"
        echo "if [ -f ${SOMVCF} ]; then rm -r ${SOMVCF}.bgz.tbi ${SOMFULLVCF}.bgz.tbi ${RAWVCF}.idx purity${PURITY}.gridss.assembly* GRIDSS_wdir_purity${PURITY}; touch logs/GRIDSS_cleanup_${PURITY}.done; fi" | qsub -cwd -hold_jid GRIDSS_filter_${PURITY}_${RAND} -l h_rt=0:10:0 -l h_vmem=1G -e logs/GRIDSS_cleanup_${PURITY}_${RAND}.err -o logs/GRIDSS_cleanup_${PURITY}_${RAND}.o -N GRIDSS_cleanup_${PURITY}_${RAND}
        else echo "Directory looks clean, skipping to overlap"
        fi
        
        if [ $NOOLAP_FLAG ]; then
          echo "Skipping overlap with truthset"
          continue; fi
            
        if [ ! -f logs/truth_overlap_${PURITY}.done ]; then
        echo "Truthset overlap"
        RAWTRUTH=${RAWVCF/.vcf/.truth.vcf}
        SOMTRUTH=${SOMVCF/.vcf/.truth.vcf}
        SOMFULLTRUTH=${SOMFULLVCF/.vcf/.truth.vcf}
        TRUTHPUR=truth.purity${PURITY}.vcf

        echo "python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $RAWVCF --file2 $TRUTHSET > $RAWTRUTH; \
        python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $SOMVCF --file2 $TRUTHSET > $SOMTRUTH; \
        python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $SOMFULLVCF --file2 $TRUTHSET > $SOMFULLTRUTH; \
        python $OVERLAP_SCRIPT --distance 100 --annotation ${PURITY}_RAW --input $TRUTHSET --file2 $RAWVCF | python $OVERLAP_SCRIPT --distance 100 --annotation ${PURITY}_SOM --input /dev/stdin --file2 $SOMVCF | \
        python $OVERLAP_SCRIPT --distance 100 --annotation ${PURITY}_FULL --input /dev/stdin --file2 $SOMFULLVCF > $TRUTHPUR;
        touch logs/truth_overlap_${PURITY}.done" | \
        qsub -cwd -hold_jid GRIDSS_cleanup_${PURITY}_${RAND} -N truthset_overlap_${PURITY}_${RAND} -l h_rt=2:0:0 -e logs/truthset_overlap_${PURITY}_${RAND}.err -o logs/truthset_overlap_${PURITY}_${RAND}.o
        else echo "Overlap files exist, skipping"
        fi
        
        
        
        
    elif [ $MODE == "nanopore" ]; then
        NSV_TUMVCF=${PURBAM/.bam/.tumor.raw.nanosv.vcf}
        SNF_TUMVCF=${PURBAM/.bam/.tumor.raw.sniffles.vcf}
        MERGEDVCF=${PURBAM/.bam/.merged.vcf}
        SOMVCF=${PURBAM/.bam/.somatic.vcf}
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
            echo "$SNIFFLES -m $PURBAM -v $SNF_TUMVCF -t 4 --report_BND --genotype" | \
            qsub -hold_jid merge_${PURITY}_${RAND} -cwd -pe threaded 4 -l h_rt=24:0:0 -l h_vmem=20G -N SNIFFLES_${PURITY}_${RAND} -e logs/SNIFFLES_${PURITY}_${RAND}.err -o logs/SNIFFLES_${PURITY}_${RAND}.o
        else echo "$SNF_TUMVCF exists, skipping SV calling"
        fi
        
        if [ ! -f $NSV_TUMVCF ]; then
            echo ". $NANOSV_VENV/bin/activate; $NANOSV_VENV/bin/NanoSV -t 4 -c $NANOSV_CONFIG -b $NANOSV_VENV/bin/human_hg19.bed -o $NSV_TUMVCF $PURBAM" | \
            qsub -hold_jid merge_${PURITY}_${RAND} -cwd -pe threaded 4 -l h_rt=96:0:0 -l h_vmem=80G -N NANOSV_${PURITY}_${RAND} -e logs/NANOSV_${PURITY}_${RAND}.err -o logs/NANOSV_${PURITY}_${RAND}.o
        else echo "$NSV_TUMVCF exists, skipping SV calling"
        fi
        
        
        if [ ! -f $SOMVCF ]; then
            echo "$NSV_NORVCF" > purity${PURITY}.nanopore.list
            echo "$SNF_NORVCF" >> purity${PURITY}.nanopore.list
            echo "$NSV_TUMVCF" >> purity${PURITY}.nanopore.list
            echo "$SNF_TUMVCF" >> purity${PURITY}.nanopore.list
            echo "$SURVIVOR merge purity${PURITY}.nanopore.list 200 1 0 0 0 0 $MERGEDVCF; grep -P '^#|VEC=00[01][01];' $MERGEDVCF > ${SOMVCF}" | \
            qsub -hold_jid SNIFFLES_${PURITY}_${RAND} -hold_jid SNIFFLES_NORMAL_${RAND} -hold_jid NANOSV_${PURITY}_${RAND} -hold_jid NANOSV_NORMAL_${RAND} -N NANOPORE_filter_${PURITY}_${RAND} -cwd -l h_rt=2:0:0 -e logs/NANOPORE_filter_${PURITY}_${RAND}.err -o logs/NANOPORE_filter_${PURITY}_${RAND}.o
        else echo "$SOMVCF exists, skipping SV filtering"
        fi
        
        if [ $NOOLAP_FLAG ]; then
          echo "Skipping overlap with truthset"
          continue; fi
        
        if [ ! -f logs/truth_overlap_${PURITY}.done ]; then
        RAWTRUTH=${MERGEDVCF/.vcf/.truth.vcf}
        SOMTRUTH=${SOMVCF/.vcf/.truth.vcf}
        TRUTHPUR=truth.purity${PURITY}.vcf
        
        echo "python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $MERGEDVCF --file2 $TRUTHSET > $RAWTRUTH; \
        python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $SOMVCF --file2 $TRUTHSET > $SOMTRUTH; \
        python $OVERLAP_SCRIPT --distance 100 --annotation ${PURITY}_RAW --input $TRUTHSET --file2 $MERGEDVCF | python $OVERLAP_SCRIPT --distance 100 --annotation ${PURITY}_SOM --input /dev/stdin --file2 $SOMVCF > $TRUTHPUR; \
        touch logs/truth_overlap_${PURITY}.done" | \
        qsub -cwd -hold_jid  NANOPORE_filter_${PURITY}_${RAND} -N truthset_overlap_${PURITY}_${RAND} -l h_rt=2:0:0 -e logs/truthset_overlap_${PURITY}_${RAND}.err -o logs/truthset_overlap_${PURITY}_${RAND}.o
        else echo "Overlap files exist, skipping"
        fi
    
    
    
    
    elif [ $MODE == "pbsv" ]; then
        GET_SOMATIC=${SCRIPT_DIR}/get_somatic_pacbio.py
        TUMSIG=${PURBAM/.bam/.tumor.svsig.gz}
        NORNAME=`basename $NORMAL_SV`
        NORSIG=${NORNAME/.bam/.normal.svsig.gz}
        MERGVCF=${PURBAM/.bam/.merged.vcf}
        SOMVCF=${PURBAM/.bam/.somatic.vcf}
        if [ ! -f logs/pbsv_normal.run ]; then
            echo "$PBSV discover $NORMAL_SV $NORSIG" | \
            qsub -cwd -l h_rt=24:0:0 -l h_vmem=60G -N PBSV_DISC_NORMAL_${RAND} -e logs/PBSV_DISC_NORMAL_${RAND}.err -o logs/PBSV_DISC_NORMAL_${RAND}.o
            touch logs/pbsv_normal.run
        else echo"$NORSIG exists, or the job is running. Skipping Normal SV discovery"
        fi
        
        if [ ! -f $TUMSIG ]; then
            echo "$PBSV discover $PURBAM $TUMSIG" | \
            qsub -hold_jid merge_${PURITY}_${RAND} -cwd -l h_rt=24:0:0 -l h_vmem=60G -N PBSV_DISC_${PURITY}_${RAND} -e logs/PBSV_DISC_${PURITY}_${RAND}.err -o logs/PBSV_DISC_${PURITY}_${RAND}.o
        else echo "$TUMSIG exists, skipping SV discovery"
        fi
        
        if [ ! -f $MERGVCF ]; then
            echo "$PBSV call -j 4 $REF $NORSIG $TUMSIG $MERGVCF" | \
            qsub -cwd -l h_rt=24:0:0 -l h_vmem=60G -pe threaded 4 -hold_jid PBSV_DISC_${PURITY}_${RAND} -hold_jid PBSV_DISC_NORMAL_${RAND} -N PBSV_CALL_${PURITY}_${RAND} -e logs/PBSV_CALL_${PURITY}_${RAND}.err -o logs/PBSV_CALL_${PURITY}_${RAND}.o
        else echo "$MERGVCF exists, skipping SV calling"
        fi
        
        if [ ! -f $SOMVCF ]; then
            echo "python $GET_SOMATIC $MERGVCF | awk -F $'\t' '{OFS = FS}{\$10=\"\"; print \$0}' > $SOMVCF" | \
            qsub -l h_rt=1:0:0 -cwd -hold_jid PBSV_CALL_${PURITY}_${RAND} -N PBSV_FILTER_${PURITY}_${RAND} -e logs/PBSV_FILTER_${PURITY}_${RAND}.err -o logs/PBSV_FILTER_${PURITY}_${RAND}.o
        else echo "$SOMVCF exists, skipping SV filtering"
        fi
        
        if [ $NOOLAP_FLAG ]; then
          echo "Skipping overlap with truthset"
            continue; fi
        
        if [ ! -f logs/truth_overlap_${PURITY}.done ]; then
        RAWTRUTH=${MERGVCF/.vcf/.truth.vcf}
        SOMTRUTH=${SOMVCF/.vcf/.truth.vcf}
        TRUTHPUR=truth.purity${PURITY}.vcf
        
        echo "python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $MERGVCF --file2 $TRUTHSET > $RAWTRUTH; \
        python $OVERLAP_SCRIPT --distance 100 --annotation TRUTHSET --input $SOMVCF --file2 $TRUTHSET > $SOMTRUTH; \
        python $OVERLAP_SCRIPT --distance 100 --annotation ${PURITY}_RAW --input $TRUTHSET --file2 $MERGVCF | python $OVERLAP_SCRIPT --distance 100 --annotation ${PURITY}_SOM --input /dev/stdin --file2 $SOMVCF > $TRUTHPUR; \
        touch logs/truth_overlap_${PURITY}.done" | \
        qsub -cwd -hold_jid PBSV_FILTER_${PURITY}_${RAND} -N truthset_overlap_${PURITY}_${RAND} -l h_rt=2:0:0 -e logs/truthset_overlap_${PURITY}_${RAND}.err -o logs/truthset_overlap_${PURITY}_${RAND}.o
        else echo "Overlap files exist, skipping"
        fi
    fi
done


