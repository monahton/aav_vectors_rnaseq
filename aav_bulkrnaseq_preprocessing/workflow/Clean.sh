#!/bin/bash

# RNA-seq mini-pipeline: QC → alignment → counting

# -------------------- Setup --------------------


# Define project structure relative to current location
PROJECT_DIR="./data_pre_processing"
mkdir -p "${PROJECT_DIR}"/{raw,fastq,aligned,counts,logs,qc,hisat2_index}


cd "${PROJECT_DIR}/raw"


# Group SRA run IDs by biological sample 
UT_rep1=(SRR27673943)   
UT_rep2=(SRR27673944)   
UT_rep3=(SRR27673945)  
UT_rep4=(SRR27673949)
AAV9_rep1=(SRR27673950)   
AAV9_rep2=(SRR27673951)   
AAV9_rep3=(SRR27673956)

# -------------------- Download & Convert --------------------

# Download .sra files
for file in "${UT_rep1[@]}" "${UT_rep2[@]}" "${UT_rep3[@]}" "${UT_rep4[@]}" "${AAV9_rep1[@]}" "${AAV9_rep2[@]}" "${AAV9_rep3[@]}"; do
  prefetch "$file"
done

# Convert to gzipped FASTQ
for file in "${UT_rep1[@]}" "${UT_rep2[@]}" "${UT_rep3[@]}" "${UT_rep4[@]}" "${AAV9_rep1[@]}" "${AAV9_rep2[@]}" "${AAV9_rep3[@]}"; do
  fasterq-dump -e 16 -p -O . "$file"
  gzip -f "${file}.fastq"
done

# Concatenate per-sample FASTQs
cat "${UT_rep1[@]/%/.fastq.gz}"  > UT1.fastq.gz
cat "${UT_rep2[@]/%/.fastq.gz}"  > UT2.fastq.gz
cat "${UT_rep3[@]/%/.fastq.gz}" > UT3.fastq.gz
cat "${UT_rep4[@]/%/.fastq.gz}" > UT4.fastq.gz
cat "${AAV9_rep1[@]/%/.fastq.gz}" > AAV9_1.fastq.gz
cat "${AAV9_rep2[@]/%/.fastq.gz}" > AAV9_2.fastq.gz
cat "${AAV9_rep3[@]/%/.fastq.gz}" > AAV9_3.fastq.gz

# Move to fastq/ folder
mv UT*.fastq.gz AAV9*.fastq.gz ../fastq/

# -------------------- QC --------------------

cd ../fastq
fastqc UT*.fastq.gz AAV9*.fastq.gz -o ../qc --threads 16

# -------------------- Alignment (hisat) --------------------

curl -O ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg38.tar.gz

mkdir hisat2_index
cd hisat2_index
tar -xzf hg38.tar.gz

cd hisat2_index
SAMPLES=(UT1, UT2, UT3, UT4, AAV9_1, AAV9_2, AAV9_3)

# Align each sample
for sample in "${SAMPLES[@]}"
do
  echo "Aligning ${sample}..."
  hisat2 -p 4 \
    -x hisat2_index/hg38/genome \
    -U fastq/${sample}.fastq.gz \
    2> logs/${sample}_hisat2.log | \
    samtools sort -@ 4 -o aligned/${sample}.bam
  samtools index aligned/${sample}.bam
  echo "${sample} alignment done."
done


# -------------------- Quantification (featureCounts) --------------------

cd ..
curl -L -o gencode.v35.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz
gunzip -f gencode.v35.annotation.gtf.gz

featureCounts -T 16 -t exon -g gene_name \
  -a gencode.v35.annotation.gtf \
  -o counts/raw_counts_gene_sym.txt aligned/*.bam \
  &> logs/featureCounts_gene_sym.log

# Format counts matrix
{ printf "GeneSymbol\t"; head -n 2 counts/raw_counts_gene_sym.txt | tail -n 1 | cut -f7-; } > counts/final_counts_symbols.tsv
tail -n +3 counts/raw_counts_gene_sym.txt | \
  awk -v OFS="\t" '{ out=$1; for(i=7;i<=NF;i++) out=out OFS $i; print out }' >> Normalized_counts/final_counts_symbols.tsv

sed -i '' '1 s|aligned/||g; 1 s|\.bam||g' counts/final_counts_symbols.tsv

# Done
echo "Pipeline complete. Output saved in: ${PROJECT_DIR}/counts/final_counts_symbols.tsv"