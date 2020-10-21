#!/bin/bash
#set -o errexit
#set -o pipefail

################################################################################
# Set default values for mapping parameters
################################################################################

nrMM="2"
mapLen="36"
gap="1000"
USE=6
KEEPLOG="1"
FILTER="0"

################################################################################
# Requirements
################################################################################

# Programs:
# * fastx-toolkit
# * bowtie
# * Jim Kent's utils

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
echo >&2 "
$(basename $0) - Mapping FASTQ file to reference genome and returning BigBed file with coordinates of mapped reads (or fragments for paired end data)

USAGE: $(basename $0) -i <FASTQ file(s)> -n <outfile basename> -s <species> -l <read lengthn> -m <max mismatches> -o <output directory> [OPTIONS]
-i    FASTQ (gzipped) file(s) read 1                                                        [required]
-I    FASTQ (gzipped) file read 2 (if PE)                                                   [required]
-t    Type (SE or PE)                                                                       [required]
-n    Basename of the output file                                                           [required]
-a    Genome assembly (hg38, hg19, dm3, etc.)                                               [required]
-l    Length to which sequenced reads will be trimmed before mapping                        [default: $mapLen]
                      (50 for human, 36 for fly)
-m    Maximal number of mismatches for mapping (usually 3 for human, 2 for fly)             [default: $nrMM]
-g    Maximal gap (in base pairs) between paired end reads                                  [default: $gap]
-f    Filters                                                                               [default: $FILTER]
            (0) Report all mapped reads
            (1) Report fragments collapsed on identical chromosome, start, end,
                  and strand information
            (A) Report both all mapped reads AND unique fragments,
                  i.e. report results from both filters 0 and 1
-r    Report mapping statistics (0/1)                                                      [default: $KEEPLOG]
-u    Number of processes to use (will be overwritten by NSLOTS on SGE, see below)         [default: $USE]
-o    Output directory                                                                     [required]

NOTES:
Note that you need the bowtie index files for your genome assembly.
In addition the programs bowtie, fastx-toolkit and Jim Kent's utils must be installed on the machine.
"
exit 1
fi


################################################################################
# Parse input and error checking
################################################################################

while getopts "i:I:t:n:l:m:g:a:f:r:u:o:" op
do
    case "$op" in
        i)  fastqFile1="$OPTARG";;
        I)  fastqFile2="$OPTARG";;
        t)  type="$OPTARG";;
        n)  sample="$OPTARG";;
        l)  mapLen="$OPTARG";;
        m)  nrMM="$OPTARG";;
        g)  gap="$OPTARG";;
        f)  FILTER="$OPTARG";;
        r)  KEEPLOG="$OPTARG";;
        u)  USE="$OPTARG";;
        a)  assembly="$OPTARG";;
        o)  outDir="$OPTARG";;
        \?) exit 1;;
    esac
done

if [ -z "$fastqFile1" ]; then
    echo >&2 "ERROR: Specify input FASTQ file R1!"
    exit 1
fi

if [ -z "$fastqFile2" ] && [ "$type" == "PE" ]; then
    echo >&2 "ERROR: Specify input FASTQ file R2 for PE data!"
    exit 1
fi

if [ -z "$mapLen" ]; then
    echo >&2 "ERROR: Specify desired length of the read for mapping!"
    exit 1
fi

if [ -z "$sample" ]; then
    echo >&2 "ERROR: Specify basename of the output file!"
    exit 1
fi

if [ -z "$assembly" ]; then
    echo >&2 "ERROR: Specify genome assembly!"
    exit 1
fi

if [ -z "$outDir" ]; then
    echo >&2 "ERROR: Specify output directory!"
    exit 1
fi

echo "input: $fastqFile1 $fastqFile2; type: $type; sampleID: $sample; mapping length: $mapLen; num mismatches: $nrMM; gap for PE reads: $gap, Filter: $FILTER; KeepLog: $KEEPLOG; Use: $USE;  genome assembly: $assembly; Output Directory: $outDir"

################################################################################
# Set index and chromosome sizes
################################################################################

INDEX=/home/groups/CEDAR/anno/indices/bowtie/${assembly}/${assembly}
SIZES=/home/groups/CEDAR/anno/chromsizes/${assembly}.chrom.sizes
export PATH=$PATH:/home/groups/CEDAR/tools/kentUtils/


echo ${INDEX} ${SIZES} 

[ -e "${INDEX}.1.ebwt" ] || echo >&2 "ERROR: No bowtie index files found for genome assembly ${assembly}!"
[ -e "${SIZES}" ] || echo >&2 "ERROR: No chromosome sizes file found for genome assembly ${assembly}!"

################################################################################
# Create tmp directory and start log file
################################################################################

# A temporary directory
TMP=$(mktemp -d)

# Generating temp log file to store bowtie stats
# logs
(
  echo -n "Start: "
  date
  echo "mapReads.sh"
  bowtie --version | head -n 1
  echo "$0 $@"
  echo
) > $TMP/LOG.log


################################################################################
# Set processor usage
# Either assigned by SGE or half of the available cores
# (at least 3 as we anyway need 3 later on)
################################################################################
# what is the equivalent in slurm I think it is $SLURM_CPUS_PER_TASK

#if [ -z "$NSLOTS" ]; then
#    USE=$(grep processor /proc/cpuinfo | wc -l | awk '{if($1>6){print int($1/2)}else{print 3}}')
#else
#    USE=$NSLOTS
#fi


################################################################################
# Determine type of input - single end or paired end
################################################################################

echo $type
# For NextSeq, you receive 4 fastq files, need to combine them prior to analysis
# This separation of se and pe would need to be done at the trimming step, extraction step


################################################################################
# Trim reads to specified length
################################################################################
(echo -n "Running fastx_trimmer: "; date) >> $TMP/LOG.log

if [ "$type" == "SE" ]; then
    zcat $fastqFile1 | fastx_trimmer -Q33 -l $mapLen -o ${TMP}/${sample}.trimmed.fq
elif [ "$type" == "PE" ]; then
    zcat $fastqFile1 | fastx_trimmer -Q33 -l $mapLen -o ${TMP}/${sample}.read_1.trimmed.fq
    zcat $fastqFile2 | fastx_trimmer -Q33 -l $mapLen -o ${TMP}/${sample}.read_2.trimmed.fq
fi


################################################################################
# Map with bowtie
################################################################################
(echo -n "Running bowtie: "; date) >> $TMP/LOG.log

if [ "$type" == "SE" ]; then
    bowtie -p $USE -q -v $nrMM -m 1 --best --strata $INDEX ${TMP}/${sample}.trimmed.fq > ${TMP}/${sample}.reads.bowtie 2>> $TMP/LOG.log
    rm ${TMP}/${sample}.trimmed.fq
elif [ "$type" == "PE" ]; then
    bowtie -p $USE -q -X $gap -v $nrMM -m 1 --best --strata $INDEX -1 ${TMP}/${sample}.read_1.trimmed.fq -2 ${TMP}/${sample}.read_2.trimmed.fq > ${TMP}/${sample}.reads.bowtie 2>> $TMP/LOG.log
    # new: save output as sam file, copy output to samples/bams, sort, save as bam
    bowtie -S -p $USE -q -X $gap -v $nrMM -m 1 --best --strata $INDEX -1 ${TMP}/${sample}.read_1.trimmed.fq -2 ${TMP}/${sample}.read_2.trimmed.fq > ${TMP}/${sample}.sam 2>> $TMP/LOG.log
    cp ${TMP}/${sample}.sam samples/bams
    samtools view -@ ${USE} -bS samples/bams/${sample}.sam > samples/bams/${sample}.bam
    samtools sort samples/bams/${sample}.bam -@ 4 -m 4G > samples/bams/${sample}.sorted.bam
    samtools index samples/bams/${sample}.sorted.bam
    # rm files in tmp space + unsorted bam
    rm ${TMP}/${sample}.sam
    rm samples/bams/${sample}.bam
    rm samples/bams/${sample}.sam
    # end new
    rm ${TMP}/${sample}.read_1.trimmed.fq ${TMP}/${sample}.read_2.trimmed.fq
fi
# -q ensures fastq input/output

################################################################################
# Convert output to BED and sort
################################################################################
(echo -n "Parsing & sorting output: "; date) >> $TMP/LOG.log

if [ "$type" == "SE" ]; then
    
    cat ${TMP}/${sample}.reads.bowtie | \
     awk -vFS="\t" -vOFS="\t" -vR=$mapLen '
        {
            print $3,$4,$4+R,".",0,$2
        }' | \
     sort -k1,1 -k2,2n -k6,6 -S 1G --parallel $(echo $USE | awk '{print $1-2}') > ${TMP}/${sample}.reads.sorted.bed

elif [ "$type" == "PE" ]; then

    cat ${TMP}/${sample}.reads.bowtie | \
     awk -vFS="\t" -vOFS="\t" '{
        if(NR%2==1)
        {
            split($1,T,"/")
            R[T[2]]=T[1]

            if (T[2]==1)
            {
                S="+"
            }
            else
            {
                S="-"
            }

            B=$4

        }
        else
        {
            print $3,B,$4+length($5),".",0,S

        }
    }' | \
    sort -k1,1 -k2,2n -k6,6 -S 1G --parallel $(echo $USE | awk '{print $1-2}') > ${TMP}/${sample}.reads.sorted.bed

fi

rm ${TMP}/${sample}.reads.bowtie

################################################################################
# Apply filters to collapse reads and/or remove information:
#
#   (0) Report all mapped reads
#   (1) Report fragments collapsed on identical chromosome, start, end,
#       and strand information
#   (A) Report both all mapped reads AND unique fragments,
#       i.e. report results from both filters 0 and 1
#
################################################################################
(echo -n "Applying filters: "; date) >> $TMP/LOG.log

if [ $FILTER == "0" ]; then
    bedToBigBed ${TMP}/${sample}.reads.sorted.bed $SIZES ${TMP}/${sample}.all.bb 2>> ${TMP}/LOG.log
    echo >> ${TMP}/LOG.log
    mv ${TMP}/${sample}.all.bb ${outDir}/${sample}.all.bb
elif [ $FILTER == "1" ]; then
    awk -vOFS='\t' '(!x[$1" "$2" "$3" "$6]++){print $0}' ${TMP}/${sample}.reads.sorted.bed > ${TMP}/${sample}.reads.filtered.bed
    bedToBigBed ${TMP}/${sample}.reads.filtered.bed $SIZES ${TMP}/${sample}.bb 2>> ${TMP}/LOG.log
    echo >> ${TMP}/LOG.log
    mv ${TMP}/${sample}.bb ${outDir}/${sample}.bb
elif [ $FILTER == "A" ]; then
    bedToBigBed ${TMP}/${sample}.reads.sorted.bed $SIZES ${TMP}/${sample}.all.bb 2>> ${TMP}/LOG.log
    awk -vOFS='\t' '(!x[$1" "$2" "$3" "$6]++){print $0}' ${TMP}/${sample}.reads.sorted.bed > ${TMP}/${sample}.reads.filtered.bed
    bedToBigBed ${TMP}/${sample}.reads.filtered.bed $SIZES ${TMP}/${sample}.bb 2>> ${TMP}/LOG.log
    echo >> ${TMP}/LOG.log
    mv ${TMP}/${sample}.bb ${outDir}/${sample}.bb
    mv ${TMP}/${sample}.all.bb ${outDir}/${sample}.all.bb
fi    
# filter=1 line !x[$1" "$2" "$3" "$6]++ means: do not keep reads with the same chr, start and strand --> remove duplicates

(echo -n "End: "
date ) >> $TMP/LOG.log

# clean up
if [ $FILTER == "0" ]; then
    rm ${TMP}/${sample}.reads.sorted.bed
elif [ $FILTER == "1" ]; then
    rm ${TMP}/${sample}.reads.sorted.bed
    rm ${TMP}/${sample}.reads.filtered.bed
elif [ $FILTER == "A" ]; then
    rm ${TMP}/${sample}.reads.sorted.bed
    rm ${TMP}/${sample}.reads.filtered.bed
fi    


if [ $KEEPLOG == "1" ]; then
    mv ${TMP}/LOG.log ${outDir}/${sample}.log
fi

################################################################################
# Clean temp folder and exit
################################################################################

rm -rf ${TMP}

exit 0
