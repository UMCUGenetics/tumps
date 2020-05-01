#!/bin/bash

usage() {
echo "
Required parameters:
    -t|--tumor_bam		    
    -n|--normal_bam		

Optional parameters:
    -r|--reference  Reference
    -bl|--blacklist  Blacklist
    -g|--gridss_jar GRIDSS jar
    -a|--assembly   Assembly
    -o|--output     Output


"
}

POSITIONAL=()

#DEFAULT
#BLACKLIST=/hpc/cog_bioinf/kloosterman/users/sdeblank/data/NextGen_Blacklisted_Regions_ENCFF001TDO_chr.bed
REFERENCE=/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta
GRIDSS_JAR=/gnu/store/iyw3yd1gl36ga89rb1sbn4cpz90jplf2-gridss-1.8.0/share/java/gridss/gridss.jar
OUTDIR=.

while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
    -h|--help)
    usage
    shift # past argument
    ;;
    -t|--tumor_bam)
    TUMOR_BAM="$2"
    shift # past argument
    shift # past value
    ;;
    -n|normal_bam)
    NORMAL_BAM="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--reference)
    REFERENCE="$2"
    shift # past argument
    shift # past value
    ;;
    -bl|--blacklist)
    BLACKLIST="$2"
    shift # past argument
    shift # past value
    ;;
    -g|--gridss_jar)
    GRIDSS_JAR="$2"
    shift # past argument
    shift # past value
    ;;
    -a|--assembly)
    ASSEMBLY="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z $TUMOR_BAM ]; then
    echo "Missing -t|--tumor_bam parameter"
    usage
    exit
fi

if [ -z $NORMAL_BAM ]; then
    echo "Missing -n|--normal_bam parameter"
    usage
    exit
fi

WDIR=./GRIDSS_wdir_${TUMOR_BAM/.bam/}
if [ ! -d $WDIR ]; then
    mkdir $WDIR
fi

TMPDIR=$WDIR/GRIDSS_tmp_${TUMOR_BAM/.bam/}

if [ ! -d $TMPDIR ]; then
    mkdir $TMPDIR
fi

if [ -z $OUTPUT ]; then
  OUTPUT=${TUMOR_BAM/.bam/.sv.vcf}
fi

if [ -z $ASSEMBLY ]; then
  ASSEMBLY=${TUMOR_BAM/.bam/.gridss.assembly.bam}
fi


if ! which bwa >/dev/null 2>&1 ; then
	echo "Missing bwa. Please add to PATH"
	exit 1
fi
if [[ ! -f "$REFERENCE" ]] ; then
	echo "Missing reference genome $REFERENCE. Update the REFERENCE variable in the shell script to your hg19 location"
	echo "For the example file chr12.1527326.DEL1024.bam, ReorderSam can be used to match to your version of hg19. In the case of this example, only \"chr12\" is required to exist, and difference in alternate contigs can be ignored (using ALLOW_INCOMPLETE_DICT_CONCORDANCE=true)."
	echo "For real data, please ensure that all BAM files are aligned to the same reference, and the reference supplied to GRIDSS matches that used for alignment."
	exit 1
fi
if [[ ! -f "$REFERENCE.bwt" ]] ; then
	echo "Missing bwa index for $REFERENCE. Could not find $REFERENCE.bwt. Create a bwa index (using \"bwa index $REFERENCE\") or symlink the index files to the expected file names."
	exit 1
fi
if [[ ! -f $GRIDSS_JAR ]] ; then
	echo "Missing $GRIDSS_JAR. Update the GRIDSS_JAR variable in the shell script to your location"
	exit 1
fi
if ! which java >/dev/null 2>&1 ; then
	echo "Missing java. Please add java 1.8 or later to PATH"
	exit 1
fi
JAVA_VERSION="$(java -version 2>&1 | grep version )"
if [[ ! "$JAVA_VERSION" =~ "\"1.8" ]] ; then
	echo "Detected $JAVA_VERSION. GRIDSS requires Java 1.8 or later."
	exit 1
fi
if ! which Rscript >/dev/null 2>&1 ; then
	echo "Missing R installation. Please add Rscript to PATH"
	exit 1
fi

ulimit -n $(ulimit -Hn) # Reduce likelihood of running out of open file handles
unset DISPLAY # Prevents errors attempting to connecting to an X server when starting the R plotting device

# -Dreference_fasta is only required for CRAM input files
# -Dgridss.gridss.output_to_temp_file=true allows GRIDSS to continue where it left off without data errors due to truncated files
# -Dsamjdk.create_index=true is required for multi-threaded operation
# -Dsamjdk.use_async_io allow for async read/write on background threads which improves BAM I/O performancce
#-Dreference_fasta=$REFERENCE \

/gnu/store/1hykmyl04mhvrwd5qrz88ymamj7nhc1p-icedtea-3.7.0/bin/java -ea -Xmx64g \
	-Dsamjdk.create_index=true \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.compression_level=1 \
	-Dgridss.gridss.output_to_temp_file=true \
	-cp $GRIDSS_JAR gridss.CallVariants \
	WORKER_THREADS=4 \
	TMP_DIR=$TMPDIR \
	WORKING_DIR=$WDIR \
	REFERENCE_SEQUENCE=$REFERENCE \
	INPUT=$TUMOR_BAM \
	INPUT=$NORMAL_BAM \
	OUTPUT=$OUTPUT \
	ASSEMBLY=$ASSEMBLY \
# 	#BLACKLIST="$BLACKLIST" \
2>&1 | tee -a gridss.$HOSTNAME.$$.log

echo "Done"
