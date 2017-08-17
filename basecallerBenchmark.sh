#!/bin/bash

# Load config file
source basecallerBenchmark.config
ALIGN=true; INDEX_REF=true; USECLUSTER=false

# Override defaults with given options
argCount=0
while getopts "b:si:r:t:h" opt; do
case $opt in
	b) BAM_FILE=$OPTARG; ALIGN=false; ((argCount+=2));; # skip alignment to reference, bam-file per directory must be available
	s) INDEX_REF=false; ((argCount++));; # Skip reference indexing, fai-file must be available
	i) INT=$OPTARG; ((argCount+=2));; # Location of intermediate folder
	r) REF_FILE=$OPTARG; ((argCount+=2));;
	t) NTHREADS=$OPTARG; ((argCount+=2));;
	q) USECLUSTER=true; ((argCount++));; # Use cluster through qsub for analysis of pileup file. Recommended for large genomes!
	h)
	echo "

Basecaller benchmark: output MinION basecalling and assembly quality metrics

Usage:bash basecallerBenchmark.sh [ -options ] readsDirectory1 [ readsDirectory2 ...]

NOTE: reads Directories must have meaningful and different names!

options:
-i : location at which intermediate files are saved [mandatory].
-r : reference file [mandatory].
-t : number of threads (standard 1).
-b : Provide alignment bam-file and skip alignment to reference.
-s : skip reference fasta indexing by BWA-MEM. Index files must be available.

"
	exit 0
	;;
esac
done
shift $argCount

# Create log file
mkdir -p logFiles
LOG_FILE="logFiles/basecallerBenchmark_run_$(date +%Y-%m-%d_%H:%M).log"
exec 3>&1 1>>${LOG_FILE} 2>&1


# Create subdirectories for intermediate results, if required
mkdir -p $INT
mkdir -p $INT/basecallerMetrics
mkdir -p $INT/assemblerMetrics
mkdir -p $INT/assemblerMetrics/quast
mkdir -p $INT/denovoAssemblies
mkdir -p $INT/images

# REFERENCE PRE-PROCESSING -------------------------------------------------------------------------
if $INDEX_REF;then 
	echo "Indexing ref file..." | tee /dev/fd/3
	$BWA index $REF_FILE
else 
	echo "Skipping reference indexing..." | tee /dev/fd/3
fi

$JELLYFISH count -m 5 -s 100M -t $NTHREADS -C $REF_FILE -o $INT/assemblerMetrics/kmerCount_reference.jf
$JELLYFISH dump -c $INT/assemblerMetrics/kmerCount_reference.jf > $INT/assemblerMetrics/kmerCount_reference.txt

# START BENCHMARK ----------------------------------------------------------------------------------
while [ "$1" ]; do
READ_DIR=$1; RD=$(basename $READ_DIR)
INT_RD=$INT/$RD
mkdir -p $INT_RD
echo "$(date +%Y-%m-%d_%H:%M) Processing read directory $RD" | tee /dev/fd/3

# READ PRE-PROCESSING-------------------------------------------------------------------------------

# If file supplied, check if type is fastq
if  [ -f "$READ_DIR" ]; then
	checkFastq=`head -n3 $READ_DIR | grep -Pq "@.+\n.+\n\+"`
	if grep -Pq "^@" $READ_DIR ;then # TODO needs a more robust fastq-check
		echo "(Multi-)fastq file with reads supplied." | tee /dev/fd/3
		READS_FILE=$READ_DIR
	else
		echo "Supplied reads file does not seem to be of fastq format. Skipping this readsfile." | tee /dev/fd/3
		shift; continue
	fi
fi

if [ -d "$READ_DIR" ]; then
	# Convert files in directory to multifastq
	READS_FILE=$INT_RD/reads_$RD.fq
	python3 extract_fastq.py $READ_DIR $READS_FILE
fi

# ALIGNMENT TO REFERENCE----------------------------------------------------------------------------
if $ALIGN; then
	echo "$(date +%Y-%m-%d_%H:%M) Aligning reads to reference genome..." | tee /dev/fd/3
	BAM_FILE=$INT_RD/bwaAlignment_${RD}.bam
	# Align to reference with bwa mem
	$BWA mem -t $NTHREADS -k $MIN_SEED_LENGTH -w $BANDWIDTH -d $Z_DROPOFF -r $SEED_SPLIT_RATIO -O $GAP_OPEN_PENALTY -E $GAP_EXTENSION_PENALTY -c $MAX_OCCURRENCE -A $MATCHING_SCORE -B $MISMATCH_PENALTY -L $CLIPPING_PENALTY -U $UNPAIR_PENALTY -x $ONT_FLAG  $REF_FILE $READS_FILE | $SAMTOOLS view -b > $BAM_FILE
fi

# EVALUATE BASECALLING QUALITY----------------------------------------------------------------------
# Create pileup file
echo "$(date +%Y-%m-%d_%H:%M) Creating pileup-file..." | tee /dev/fd/3
PILEUP_FILE=${BAM_FILE}.pileup
$SAMTOOLS sort -@ $NTHREADS -o $BAM_FILE $BAM_FILE 
$SAMTOOLS mpileup -f $REF_FILE $BAM_FILE > $PILEUP_FILE

# Generate basecaller benchmark statistics
echo "$(date +%Y-%m-%d_%H:%M) Generating basecaller quality metrics..." | tee /dev/fd/3
if $USECLUSTER; then
	echo "Preparing to send job to cluster..."
	CLUSTER_DIR=$INT/clusterfiles/$RD
	CLUSTER_PILEUPS=$CLUSTER_DIR/pileups
	CLUSTER_SCRIPTS=$CLUSTER_DIR/scripts
	CLUSTER_OUTPUT=$CLUSTER_DIR/output
	mkdir -p $CLUSTER_PILEUPS
	mkdir -p $CLUSTER_SCRIPTS
	mkdir -p $CLUSTER_OUTPUT

	# Split pileup into non-overlapping subfiles
	split --numeric-suffixes=1 --suffix-length=9 -l $REPEATS_NBASES $PILEUP_FILE $CLUSTER_PILEUPS/
	# Create scripts
	COUNTER=1
	for nRepeats in $(seq $REPEATS_MIN $REPEATS_MAX); do
		for pb in `find $CLUSTER_PILEUPS`; do
			FILE=`readlink -f $pb`
			printf "#!/bin/Rscript
pileupName =  \"$pb\";
rFileNb=$COUNTER;
i = $nRepeats;
u = $REPEATS_UNITS;
outDir=\"$CLUSTER_COUNTS\";
" > $CLUSTER_SCRIPTS/${COUNTER}.R	
			cat rScripts/basecallerBenchmarkCluster.R >> $CLUSTER_SCRIPTS/${COUNTER}.R
			chmod 777 $CLUSTER_SCRIPTS/${COUNTER}.R
			let COUNTER=COUNTER+1
		done
	done
	numPileups=`ls $CLUSTER_SCRIPTS | wc -w`
	
	printf "#!/bin/bash\n$CUSTER_SCRIPTS/\"\$SGE_TASK_ID\".sh" > $CLUSTER_DIR/qsubJob.sh
	# qsub -wd $INT/clusterfiles/$RD -tc $CLUSTER_NCORES -t 1-$numPileups -e $CLUSTER_DIR/stderr -o $CLUSTER_DIR/stdout -S /bin/bash $CLUSTER_DIR/qsubJob.sh
else
	Rscript rScripts/basecallerBenchmarkMaster.R $PILEUP_FILE $INT $NTHREADS
fi
# DE NOVO ALIGNMENT---------------------------------------------------------------------------------
echo "$(date +%Y-%m-%d_%H:%M) Aligning reads de novo..." | tee /dev/fd/3

MAPPING_PAF=$INT_RD/minimap_${RD}.paf.gz
ASSEMBLY_GFA=$INT_RD/miniasm_${RD}.gfa
ASSEMBLY_PAF=$INT_RD/miniasm_${RD}.paf
ASSEMBLY_GFAMOD=$INT_RD/miniasm_${RD}.gfamod
ASSEMBLY_PAFMOD=$INT_RD/miniasm_${RD}.pafmod
ASSEMBLY_FA=$INT_RD/miniasm_${RD}.fa

# Assemble with Minimap/Miniasm
$MINIMAP -x ava10k $READS_FILE $READS_FILE | gzip -1 > $MAPPING_PAF
$MINIASM -f $READS_FILE $MAPPING_PAF > $ASSEMBLY_GFA
$MINIASM -p paf -f $READS_FILE $MAPPING_PAF > $ASSEMBLY_PAF
# Miniasm adds clipping positions to read names. Turn into extra columns to avoid changing read name.
sed -r 's/:/\t/g' $ASSEMBLY_GFA | perl -ne 's/(?<=[0-9])-(?=[0-9]+    )/      /g; print;' > $ASSEMBLY_GFAMOD
sed -r 's/:/\t/g' $ASSEMBLY_PAF | perl -ne 's/(?<=[0-9])-(?=[0-9]+    )/      /g; print;' > $ASSEMBLY_PAFMOD
awk '/^S/{print ">"$2"\n"$3}' $ASSEMBLY_GFA | fold > $ASSEMBLY_FA # Extract fasta from miniasm's gfa

# If Miniasm returns empty file, no overlaps were found
linesAssembly=` wc -l $ASSEMBLY_FA | cut -f1 -d ' '`
if [ "$linesAssembly" -eq 0 ]; then 
	echo "No overlaps found, assembly failed." | tee /dev/fd/3
	shift; continue
fi
cp $ASSEMBLY_FA $INT/denovoAssemblies/.

# EVALUATE DE NOVO ASSEMBLY QUALITY ----------------------------------------------------------------
echo "$(date +%Y-%m-%d_%H:%M) calculating k-mer frequencies..." | tee /dev/fd/3
# Count k-mers with Jellyfish
DENOVO_OVERLAPS=$INT/denovoAssemblies/miniasm_${RD}.ovl
REF_OVERLAPS=${BAM_FILE}.ovl
REF_OVERLAPS_REV=${BAM_FILE}_rev.ovl

$JELLYFISH count -m 5 -s 100M -t $NTHREADS -C $ASSEMBLY_FA -o $INT/assemblerMetrics/kmerCount_${RD}.jf
$JELLYFISH dump -c $INT/assemblerMetrics/kmerCount_${RD}.jf > $INT/assemblerMetrics/kmerCount_${RD}.txt

# End of loop, assembly quality metrics continue
shift
done

#NOTE requires matplotlib for graphs

# If de novo assemblies were created succesfully, run QUAST
if [[ `find $INT/denovoAssemblies -type f | wc -l` -gt 0 ]];then
	echo "$(date +%Y-%m-%d_%H:%M)  Running QUAST assembly benchmark..."
	$QUAST -o $INT/assemblerMetrics/quast -R $REF_FILE $INT/denovoAssemblies/*.fa
fi

# Draw graphs
Rscript rScripts/drawGraphs.R $INT
