#!/bin/bash

conda activate pytorch

# Basic parameter definitions
work_dir=~/workspace/RNA-Seq
# genome=hg38
# anot=refGene
genome_file=GRCh38.primary_assembly.genome.fa
anot_file=gencode.v47.basic.annotation.gtf
read_length=150
trim_length=10
thread=32
max_io_jobs=8
min_SJ_reads=5


# Paths
genome_fa=${work_dir}/${genome_file}
anot_gtf=${work_dir}/${anot_file}
fastq_dir=${work_dir}/fastq
align_dir=${work_dir}/align
expr_dir=${work_dir}/expr
tmp_dir=${work_dir}/tmp
rsem_reference=${work_dir}/rsem_reference/
log_dir=${work_dir}/log
fastq_suffix=fq.gz
fastqc_dir=${work_dir}/fastqc
fastq_trim_dir=${work_dir}/fastq_trim
# Create directories
mkdir -p ${align_dir} ${tmp_dir} ${expr_dir} ${log_dir} ${rsem_reference} ${fastqc_dir} ${fastq_trim_dir}


# Quality control of raw sequencing data
for i in ${fastq_dir}/*.${fastq_suffix}; do
    (
        sample=$(basename "$i")
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Performing FastQC quality check on: $sample..."
        fastqc -t ${thread} -o ${fastqc_dir} ${fastq_dir}/${sample} > ${log_dir}/fastqc_${sample}.log 2>&1
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished performing FastQC quality check on: $sample..."
    ) &
    while [[ $(jobs -r -p | wc -l) -ge $max_io_jobs ]]; do
        sleep 1
    done
done
wait


# Trim adaptor and index
for i in ${fastq_dir}/*_1.${fastq_suffix}; do
    (
        sample=$(basename "$i")
        sample=${sample%"_1.$fastq_suffix"}
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Trimming adaptor and indexing on: $sample..."
        trim_galore --output_dir ${fastq_trim_dir} --gzip --clip_R1 $trim_length --clip_R2 $trim_length \
            --fastqc --fastqc_args "-t ${thread} -o ${fastqc_dir}" \
            --paired ${fastq_dir}/${sample}_1.$fastq_suffix ${fastq_dir}/${sample}_2.$fastq_suffix \
            > ${log_dir}/trim_${sample}.log 2>&1
        mv ${fastq_trim_dir}/${sample}_1_val_1.fq.gz ${fastq_trim_dir}/${sample}_1.fq.gz
        mv ${fastq_trim_dir}/${sample}_2_val_2.fq.gz ${fastq_trim_dir}/${sample}_2.fq.gz
        date=$(date "+%Y-%m-%d %H:%M:%S")
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished trimming adaptor and indexing on $sample."
    ) &
    while [[ $(jobs -r -p | wc -l) -ge $max_io_jobs ]]; do
        sleep 1
    done
done
wait
fastqc_dir=${fastq_trim_dir}


# Generate initial STAR index (approximately 2 hours and 10 minutes)
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Indexing genome with STAR (First Indexing)..."
STAR \
    --runMode genomeGenerate \
    --genomeDir ${work_dir}/star_genome_index \
    --genomeFastaFiles ${genome_fa} \
    --sjdbOverhang $((read_length - trim_length - 1)) \
    --sjdbGTFfile ${anot_gtf} \
    --runThreadN ${thread} \
    > ${log_dir}/index_1.log 2>&1
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished indexing genome with STAR (First Indexing)..."


# First alignment (approximately 5 minutes per sample)
if [[ ${fastq_suffix} == *gz ]]; then
    read_method="zcat"
else
    read_method="cat"
fi
> ${align_dir}/combined_SJ.out.tab

for i in ${fastq_dir}/*_1.${fastq_suffix}; do
    sample=$(basename "$i")
    sample=${sample%"_1.$fastq_suffix"}
    # Check if input files exist
    if [ ! -f "${fastq_dir}/${sample}_1.${fastq_suffix}" ] || [ ! -f "${fastq_dir}/${sample}_2.${fastq_suffix}" ]; then
        echo "Error: Missing FASTQ files for sample ${sample}" >&2
        continue
    fi
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] First Aligning sample: $sample..."
    STAR \
        --genomeDir ${work_dir}/star_genome_index \
        --readFilesIn ${fastq_dir}/${sample}_1.${fastq_suffix} ${fastq_dir}/${sample}_2.${fastq_suffix} \
        --runThreadN ${thread} \
        --outFileNamePrefix ${tmp_dir}/${sample}. \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 1 \
        --genomeLoad NoSharedMemory \
        --readFilesCommand $read_method \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --sjdbOverhang $((read_length - trim_length - 1)) \
        --outSAMstrandField intronMotif \
        --outSAMtype None \
        --outSAMmode None \
        > ${log_dir}/alignment_1_${sample}.log 2>&1
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished first Aligning sample: $sample..."
    # Merge SJ files
    cat ${tmp_dir}/${sample}.SJ.out.tab >> ${align_dir}/combined_SJ.out.tab
done


# Filter and sort SJ files
awk -v min_SJ_reads="$min_SJ_reads" '$5 >= min_SJ_reads && $6 > 0 {print $0}' ${align_dir}/combined_SJ.out.tab > ${align_dir}/filtered_combined_SJ.out.tab
sort -k1,1 -k2,2n -k3,3n ${align_dir}/filtered_combined_SJ.out.tab | uniq > ${align_dir}/unique_combined_SJ.out.tab


# Second STAR indexing (approximately 2 hours and 10 minutes)
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Indexing genome with STAR (Second Indexing)..."
STAR \
    --runMode genomeGenerate \
    --genomeDir ${work_dir}/star_index \
    --genomeFastaFiles $genome_fa \
    --sjdbOverhang $((read_length - trim_length - 1)) \
    --runThreadN $thread \
    --sjdbFileChrStartEnd ${align_dir}/unique_combined_SJ.out.tab \
    > ${log_dir}/index_2.log 2>&1
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished indexing genome with STAR (Second Indexing)..."


# Second alignment (approximately 10 minutes per sample)
for i in ${fastq_dir}/*_1.${fastq_suffix}; do
    sample=$(basename "$i")
    sample=${sample%"_1.$fastq_suffix"}
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Second Aligning sample: $sample..."
    STAR \
        --genomeDir ${work_dir}/star_index \
        --readFilesIn ${fastq_dir}/${sample}_1.${fastq_suffix} ${fastq_dir}/${sample}_2.${fastq_suffix} \
        --runThreadN ${thread} \
        --outFileNamePrefix ${align_dir}/${sample}. \
        --sjdbGTFfile ${anot_gtf} \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 1 \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM 0 \
        --readFilesCommand $read_method \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --sjdbOverhang $((read_length - trim_length - 1)) \
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI NM MD AS XS \
        --outSAMunmapped Within \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMheaderHD @HD VN:1.4 \
        --outSAMattrRGline ID:$sample SM:$sample PL:ILLUMINA LB:Lib_$sample \
        --quantMode GeneCounts TranscriptomeSAM \
        > ${log_dir}/alignment_2_${sample}.log 2>&1
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished second Aligning sample: $sample..."
done


# RSEM expression calculation
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Indexing genome with RSEM ..."
rsem-prepare-reference -p $thread --gtf ${anot_gtf} ${genome_fa} ${rsem_reference} \
    > ${log_dir}/rsem_index.log 2>&1

for i in ${fastq_dir}/*_1.${fastq_suffix}; do
    (
        sample=$(basename "$i")
        sample=${sample%"_1.$fastq_suffix"}
        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Calculating expression level of sample: $sample..."
        rsem-calculate-expression -p $thread --calc-ci --bam --paired-end \
            ${align_dir}/${sample}.Aligned.toTranscriptome.out.bam \
            ${rsem_reference} \
            ${expr_dir}/${sample} \
            > ${log_dir}/rsem_${sample}.log 2>&1

        echo "[$(date "+%Y-%m-%d %H:%M:%S")] Finished calculating expression level for $sample."
    ) &
    while [[ $(jobs -r -p | wc -l) -ge $max_io_jobs ]]; do
        sleep 1
    done
done
wait


# rm -rf ${tmp_dir}
echo "[$(date "+%Y-%m-%d %H:%M:%S")] Alignment and expression level calculation complete!"
