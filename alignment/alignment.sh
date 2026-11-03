
#!/bin/bash

set -Eeuo pipefail

# ============================================================
# MICB 405 RNA-seq pipeline
# Downloads FASTQs, runs FastQC, STAR, and MultiQC
# ============================================================

# Run this script from the project root directory.

THREADS=8

# 1. Directory paths

FASTQ_DIR="fastq"
RAW_QC_DIR="qc/raw"
MULTIQC_DIR="qc/multiqc"
REFERENCE_DIR="reference"
STAR_INDEX="$REFERENCE_DIR/star_index"
BAM_DIR="results/bam"
STAR_DIR="results/star_align"
LOG_DIR="logs"


# 2. Sample accessions

SAMPLES=(
    SRR2938958
    SRR2938964
    SRR2938974
    SRR2938976
    SRR2938981
    SRR2938960
    SRR2938961
    SRR2938967
    SRR2938971
    SRR2938980
    SRR2938963
    SRR2938973
    SRR2938977
    SRR2938985
    SRR2938962
    SRR2938968
    SRR2938969
    SRR2938972
)

# 3. Check required programs and reference

for program in curl gzip fastqc STAR samtools multiqc; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: $program is not installed or not in PATH."
        exit 1
    fi
done

if [[ ! -s "$REFERENCE_DIR/GRCm39.primary_assembly.genome.fa" ]]; then
    echo "ERROR: Reference genome FASTA not found."
    exit 1
fi

if [[ ! -s "$REFERENCE_DIR/gencode.vM39.primary_assembly.annotation.gtf" ]]; then
    echo "ERROR: GTF annotation not found."
    exit 1
fi

if [[ ! -s "$STAR_INDEX/Genome" ]]; then
    echo "ERROR: STAR index not found at $STAR_INDEX."
    exit 1
fi

# Reset the failure list for this run.
FAILED_FILE="$LOG_DIR/failed_samples.txt"
: > "$FAILED_FILE"

# 4. Process one sample

process_sample() {

    local SRR="$1"

    local FASTQ="$FASTQ_DIR/${SRR}.fastq.gz"

    local QC_HTML="$RAW_QC_DIR/${SRR}_fastqc.html"
    local QC_ZIP="$RAW_QC_DIR/${SRR}_fastqc.zip"

    local STAR_PREFIX="$STAR_DIR/${SRR}."
    local STAR_BAM="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"
    local STAR_LOG="${STAR_PREFIX}Log.final.out"

    local BAM="$BAM_DIR/${SRR}.sorted.bam"
    local BAI="${BAM}.bai"

    echo ""
    echo "================================================"
    echo "Processing $SRR"
    echo "================================================"

    # --------------------------------------------------------
    # STEP 1: Download FASTQ from ENA
    # --------------------------------------------------------

    if [[ -s "$FASTQ" ]] && gzip -t "$FASTQ" 2>/dev/null; then

        echo "Valid FASTQ already exists. Skipping download."

    else

        echo "Downloading $SRR..."

        rm -f "$FASTQ" "${FASTQ}.part"

        local FASTQ_URL

        FASTQ_URL=$(curl -fsSL --retry 3 --retry-delay 5 \
            "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp&format=tsv" \
            | awk -F '\t' 'NR == 2 {print $2}')

        if [[ -z "$FASTQ_URL" ]]; then
            echo "ERROR: Could not retrieve FASTQ URL for $SRR."
            return 1
        fi

        # These samples are expected to be single-end.
        if [[ "$FASTQ_URL" == *";"* ]]; then
            echo "ERROR: Multiple FASTQ files returned for $SRR."
            echo "Check whether this sample is paired-end."
            return 1
        fi

        # Convert the ENA FTP address to HTTPS.
        FASTQ_URL="${FASTQ_URL#ftp://}"
        FASTQ_URL="${FASTQ_URL#https://}"
        FASTQ_URL="https://${FASTQ_URL}"

        curl -fL --retry 3 --retry-delay 5 \
            --progress-bar "$FASTQ_URL" \
            -o "${FASTQ}.part"

        if ! gzip -t "${FASTQ}.part"; then
            echo "ERROR: Downloaded FASTQ is corrupted."
            rm -f "${FASTQ}.part"
            return 1
        fi

        mv "${FASTQ}.part" "$FASTQ"

        echo "FASTQ download complete."
    fi

    # --------------------------------------------------------
    # STEP 2: FastQC
    # --------------------------------------------------------

    if [[ -s "$QC_HTML" && -s "$QC_ZIP" ]]; then

        echo "FastQC reports already exist. Skipping."

    else

        echo "Running FastQC..."

        fastqc -t "$THREADS" \
            -o "$RAW_QC_DIR" \
            "$FASTQ"

        if [[ ! -s "$QC_HTML" || ! -s "$QC_ZIP" ]]; then
            echo "ERROR: FastQC output missing."
            return 1
        fi

        echo "FastQC complete."
    fi

    # --------------------------------------------------------
    # STEP 3: STAR alignment
    # --------------------------------------------------------

    # Check whether a valid BAM and STAR summary already exist.
    if [[ -s "$BAM" && -s "$STAR_LOG" ]] \
        && samtools quickcheck -q "$BAM"; then

        echo "Valid BAM and STAR log already exist. Skipping alignment."

    else

        echo "Running STAR alignment..."

        rm -f "$STAR_BAM"

        STAR \
            --runThreadN "$THREADS" \
            --genomeDir "$STAR_INDEX" \
            --readFilesIn "$FASTQ" \
            --readFilesCommand zcat \
            --outFileNamePrefix "$STAR_PREFIX" \
            --outSAMtype BAM SortedByCoordinate  || { echo "ERROR: STAR failed."; return 1; }

        if [[ ! -s "$STAR_BAM" ]]; then
            echo "ERROR: STAR BAM output missing."
            return 1
        fi

        if ! samtools quickcheck -q "$STAR_BAM"; then
            echo "ERROR: STAR produced an invalid BAM."
            return 1
        fi

        # Move BAM to the centralized BAM directory.
        mv "$STAR_BAM" "$BAM"

        echo "STAR alignment complete."
    fi

    # --------------------------------------------------------
    # STEP 4: Index BAM
    # --------------------------------------------------------

    if [[ -s "$BAI" ]] && samtools quickcheck -q "$BAM"; then

        echo "BAM index already exists. Skipping."

    else

        echo "Indexing BAM..."

        rm -f "$BAI"

        samtools index -@ "$THREADS" "$BAM"

        if [[ ! -s "$BAI" ]]; then
            echo "ERROR: BAM index was not created."
            return 1
        fi

        echo "BAM indexing complete."
    fi

    echo "Finished processing $SRR."

}

# ------------------------------------------------------------
# 5. Process all samples
# ------------------------------------------------------------

for SRR in "${SAMPLES[@]}"; do

    SAMPLE_LOG="$LOG_DIR/${SRR}.log"

    if process_sample "$SRR" 2>&1 | tee "$SAMPLE_LOG"; then

        echo "$SRR completed successfully."

    else

        echo "WARNING: $SRR failed. Continuing with next sample."

        echo "$SRR" >> "$FAILED_FILE"

    fi

done

# ------------------------------------------------------------
# 6. Generate MultiQC report
# ------------------------------------------------------------

echo ""
echo "================================================"
echo "Trying Next Step: MultiQC"
echo "================================================"

if [[ -s "$FAILED_FILE" ]]; then

    echo "WARNING: Some samples failed. Skipping MultiQC."
    echo "Fix the failures and rerun; completed samples will be skipped."

elif compgen -G "$RAW_QC_DIR/*_fastqc.zip" > /dev/null; then

    multiqc "$RAW_QC_DIR" \
        --outdir "$MULTIQC_DIR" \
        --force

    echo "MultiQC report generated."

else

    echo "WARNING: No FastQC reports found. Skipping MultiQC."

fi

# ------------------------------------------------------------
# 7. Final summary
# ------------------------------------------------------------

echo ""
echo "================================================"
echo "PIPELINE SUMMARY"
echo "================================================"

if [[ -s "$FAILED_FILE" ]]; then

    echo "Some samples failed:"
    cat "$FAILED_FILE"

    echo ""
    echo "Check individual logs in $LOG_DIR."

else

    echo "All samples completed successfully!"

fi

echo ""
echo "Pipeline finished."
