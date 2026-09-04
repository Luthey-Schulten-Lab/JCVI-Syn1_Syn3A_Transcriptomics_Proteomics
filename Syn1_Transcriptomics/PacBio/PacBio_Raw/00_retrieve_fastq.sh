# 1) Define the list of accessions of three PacBio datasets of Syn1
SRR_ACCESSIONS=("SRR36012643" "SRR36012642" "SRR36012641")

FASTERQ_PATH=$(which fasterq-dump 2>/dev/null)
if [ -z "$FASTERQ_PATH" ]; then
    echo "Error: fasterq-dump not found in the active conda environment. Please activate the appropriate environment."
    exit 1
fi

# Check if steps 2 and 3 can be skipped
SRR_PRESENT=true
for ID in "${SRR_ACCESSIONS[@]}"; do
    [ -f "${ID}.fastq" ] || { SRR_PRESENT=false; break; }
done

REP_PRESENT=true
for i in 1 2 3; do
    [ -f "PacBio.hifi_reads.rep${i}.fastq" ] || [ -f "PacBio.hifi_reads.rep${i}.fastq.gz" ] || { REP_PRESENT=false; break; }
done

if $REP_PRESENT; then
    echo "Rep files already present, skipping steps 2 and 3."
elif $SRR_PRESENT; then
    echo "SRR fastq files already present, skipping download (step 2). Running rename (step 3)."
    mv SRR36012643.fastq PacBio.hifi_reads.rep1.fastq
    mv SRR36012642.fastq PacBio.hifi_reads.rep2.fastq
    mv SRR36012641.fastq PacBio.hifi_reads.rep3.fastq
else
    # 2) Loop through each accession
    for ID in "${SRR_ACCESSIONS[@]}"; do
        echo "Processing $ID..."

        # Run fasterq-dump (using -e 16 for multi-threading speed)
        "$FASTERQ_PATH" --split-files --progress -e 16 "$ID"

        # Optional: Compress immediately to save disk space
        # gzip "${ID}"*.fastq
    done

    # 3) Rename the files
    # Assuming your files are named SRR36012643.fastq, etc.
    mv SRR36012643.fastq PacBio.hifi_reads.rep1.fastq
    mv SRR36012642.fastq PacBio.hifi_reads.rep2.fastq
    mv SRR36012641.fastq PacBio.hifi_reads.rep3.fastq
fi

# 4) Compress them into fastq.gz
# -p 16: use 16 threads
# -k: keep original files (optional, remove if you want to save space)
GZ_PRESENT=true
for i in 1 2 3; do
    [ -f "PacBio.hifi_reads.rep${i}.fastq.gz" ] || { GZ_PRESENT=false; break; }
done

if $GZ_PRESENT; then
    echo "Compressed rep files already present, skipping step 4."
else
    pigz -p 16 PacBio.hifi_reads.rep1.fastq
    pigz -p 16 PacBio.hifi_reads.rep2.fastq
    pigz -p 16 PacBio.hifi_reads.rep3.fastq
fi

# 5) Concatenate them into one merged file
cat PacBio.hifi_reads.rep1.fastq.gz \
    PacBio.hifi_reads.rep2.fastq.gz \
    PacBio.hifi_reads.rep3.fastq.gz > merged.hifi_reads.fastq.gz