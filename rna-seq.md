1.TPM values for tuber flesh transcriptome expression were obtained from StringTie and then merged.

hisat2 -x /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/03_genome/DM8.1_genome.hisat -1 C001R_1.filter.fq.gz -2 C001R_2.filter.fq.gz -S C001R.sam --summary-file C001R.hisat.summary

stringtie C001R.sorted.bam -e -b C001R.c -A C001R.tpm.txt -G /public/agis/yanjianbin_group/shangyi/01_potato/06_potato_data/03_genome/DM8.1_gene.gff3 -o C001R.c/C001R.gtf-p 5 

#merge tpm value with make_tpm_matrix_from_list.py



2.Transcriptome analysis of developing tubers was performed using a bash script with Salmon.


```
#!/bin/bash

# Initialization

WORK_DIR="/path/way/to/02_tuber_development_rna_seq"
REF_CDS="/path/way/to/DM8.1_gene.cds.fasta"
SRA_DIR="${WORK_DIR}"
FASTQ_DIR="${WORK_DIR}/fastq"
INDEX_DIR="${WORK_DIR}/salmon_index"
QUANT_DIR="${WORK_DIR}/salmon_quant"
MATRIX_DIR="${WORK_DIR}/salmon_matrix"

THREADS=8 # Number of threads per sample for quantification
LIB_TYPE="A" # Library type: auto-detect. Change to IU/ISF/ISR if known
SRA_TOOLKIT_PATH="" # Specify full path to fasterq-dump if sra-tools is not in PATH

# ===================================================

cd "${WORK_DIR}"
mkdir -p "${FASTQ_DIR}" "${INDEX_DIR}" "${QUANT_DIR}" "${MATRIX_DIR}"

# ==================== Step 0: SRA to FASTQ ====================

echo "[$(date)] Step 0: Converting SRA to FASTQ ..."
for sra_dir in "${SRA_DIR}"/SRR*; do
 [ -d "${sra_dir}" ] || continue
 sra_file=$(ls "${sra_dir}"/*.sra 2>/dev/null | head -1)
 [ -z "${sra_file}" ] && continue
 srr_id=$(basename "${sra_dir}")
 out_fq1="${FASTQ_DIR}/${srr_id}_1.fastq.gz"
 out_fq2="${FASTQ_DIR}/${srr_id}_2.fastq.gz"
 [ -f "${out_fq1}" ] && [ -f "${out_fq2}" ] && echo " SKIP ${srr_id}" && continue
 echo " Converting ${srr_id} ..."
 ${SRA_TOOLKIT_PATH}fasterq-dump --split-files --threads 4 "${sra_file}" -O "${FASTQ_DIR}"
 # Compress R1 and R2 in parallel to speed up
 gzip -f "${FASTQ_DIR}/${srr_id}_1.fastq" &
 gzip -f "${FASTQ_DIR}/${srr_id}_2.fastq" &
 wait
done
echo "[$(date)] Step 0 done."

# ==================== Step 1: Build Salmon index ====================

echo "[$(date)] Step 1: Building Salmon index ..."
if [ ! -f "${INDEX_DIR}/reflengths.bin" ]; then
 salmon index -t "${REF_CDS}" -i "${INDEX_DIR}" -p "${THREADS}" --kmerLen 31
else
 echo " Index already exists, skip."
fi
echo "[$(date)] Step 1 done."

# ==================== Step 2: Batch Salmon quantification ====================

echo "[$(date)] Step 2: Running Salmon quantification ..."
for fq1 in "${FASTQ_DIR}"/*_1.fastq.gz; do
 [ -f "${fq1}" ] || continue
 base=$(basename "${fq1}" _1.fastq.gz)
 fq2="${FASTQ_DIR}/${base}_2.fastq.gz"
 [ -f "${fq2}" ] || { echo " WARN: missing R2 for ${base}, skip"; continue; }
 out="${QUANT_DIR}/${base}.quant"
 [ -f "${out}/quant.sf" ] && echo " SKIP ${base}" && continue
 echo " Quantifying ${base} ..."
 salmon quant \
 -i "${INDEX_DIR}" \
 -l "${LIB_TYPE}" \
 -1 "${fq1}" -2 "${fq2}" \
 -p "${THREADS}" \
 --validateMappings --gcBias --seqBias \
 -o "${out}" \
 > "${out}.log" 2>&1
done
echo "[$(date)] Step 2 done."

# ==================== Step 3: Merge expression matrices ====================

echo "[$(date)] Step 3: Merging expression matrices ..."

# Merge TPM matrix

salmon quantmerge --quants "${QUANT_DIR}"/*.quant --column TPM \
 -o "${MATRIX_DIR}/gene.TPM.tsv"

# Merge NumReads count matrix

salmon quantmerge --quants "${QUANT_DIR}"/*.quant --column NumReads \
 -o "${MATRIX_DIR}/gene.Count.tsv"

# Remove .quant suffix from column names

sed -i '1s/\.quant//g' "${MATRIX_DIR}"/gene.*
echo "[$(date)] Step 3 done."
echo "Pipeline finished at $(date)."
```
