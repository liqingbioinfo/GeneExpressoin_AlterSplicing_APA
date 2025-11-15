#!/usr/bin/env bash

module load StdEnv/2023
module load StdEnvACCRE/2023

module load bedops/2.4.41
module load kent_tools/486
module load bedtools/2.31.0

set -euo pipefail

###############################################################
#                       ARGUMENT PARSER                       
###############################################################
usage() {
    echo ""
    echo "Usage:"
    echo "  $0 --input sample_sheet.csv --gtf annotation.gtf.gz --outdir /abspath/to/outdir --rnaseq nf-core/rnaseq"
    echo ""
    echo "Sample sheet format:"
    echo "  sample,bam"
    echo "  S001,path/to/nf-core/rnaseq/star_salmon/S001.bam"
    echo "  S002,path/to/nf-core/rnaseq/star_salmon/S002.bam"
    echo ""
    echo "Optional flags:"
    echo "  --wait slurm      Automatically wait for SLURM jobs (default)"
    echo "  --wait manual     User manually checks and continues"
    echo "  --resume none     Run full pipeline (default)"
    echo "  --resume wig      Skip steps before WIG (resume after wig generation)"
    echo "  --addchr yes|no   Add 'chr' prefix to WIG files (default=no). Match chromosomes in wig files and 3UTR_annotation.bed"
    exit 1
}

# Default parameters
THREADS=8
SLURM_CPUS=4
SLURM_MEM="40G"
SLURM_TIME="1-00:00:00"
DAPARS="/data/sbcs/GuoLab/backup/liq17/biosoft/DaPars2/src/"
WAIT_MODE="slurm"
RESUME_MODE="none"
ADDCHR="no"    ### ADDCHR DEFAULT

while [[ $# -gt 0 ]]; do
    case $1 in
        --input) SAMPLE_SHEET="$2"; shift 2 ;;
        --gtf) GTF_FILE="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --rnaseq) RNASEQ_OUT="$2"; shift 2 ;;
        --wait) WAIT_MODE="$2"; shift 2 ;;
        --resume) RESUME_MODE="$2"; shift 2 ;;
        --addchr) ADDCHR="$2"; shift 2 ;;   ### ADDCHR PARSER
        *) echo "Unknown option: $1"; usage ;;
    esac
done

[[ -z "${SAMPLE_SHEET:-}" || -z "${GTF_FILE:-}" || -z "${OUTDIR:-}" || -z "${RNASEQ_OUT:-}" ]] && usage

###############################################################
#                       DIRECTORY SETUP                       
###############################################################
mkdir -p "$OUTDIR"/{annotation,wig,slurm_logs,dapars_output,merged}

echo "======================================================"
echo "                DaPars2 APA Pipeline"
echo " Samplesheet : $SAMPLE_SHEET"
echo " GTF         : $GTF_FILE"
echo " Output dir  : $OUTDIR"
echo " nf-core/rnaseq dir  : $RNASEQ_OUT"
echo " Wait mode   : $WAIT_MODE"
echo " Resume mode : $RESUME_MODE"
echo " Add chr?    : $ADDCHR"
echo "======================================================"

# ###############################################################
# # STEP 1a — Convert GTF → BED12
# ###############################################################
# if [[ "$RESUME_MODE" == "wig" ]]; then
    # echo "[RESUME] Skipping GTF → BED12"
# else
    # echo "[1a] Converting GTF → BED12..."

    # if [[ "$GTF_FILE" == *.gz ]]; then
        # GTF_UNZIPPED="${GTF_FILE%.gz}"
        # gunzip -c "$GTF_FILE" > "$GTF_UNZIPPED"
        # GTF_FILE="$GTF_UNZIPPED"
    # fi

    # BED12="$OUTDIR/annotation/${GTF_FILE%.gtf}.bed12"

    # if [[ -f "$BED12" ]]; then
        # echo "BED12 already exists: $BED12 — skipping."
    # else
        # gtfToGenePred "$GTF_FILE" stdout \
            # | genePredToBed /dev/stdin "$BED12.tmp"

        # awk 'BEGIN{OFS="\t"} { $1="chr"$1; print }' "$BED12.tmp" > "$BED12"
        # rm "$BED12.tmp"
        # echo "Created BED12: $BED12"
    # fi
# fi

# ###############################################################
# # STEP 1b — Build transcript→gene_name ID mapping
# ###############################################################
# IDMAP="$OUTDIR/annotation/IDmapping.bed"

# if [[ "$RESUME_MODE" == "wig" ]]; then
    # echo "[RESUME] Skipping ID mapping"
# else
    # if [[ -f "$IDMAP" ]]; then
        # echo "IDMAP exists: $IDMAP — skipping."
    # else
        # echo "[1b] Creating ID mapping..."
        # echo -e "#name\tname2" > "$IDMAP"
        # awk '
            # BEGIN{OFS="\t"}
            # $3=="transcript"{
                # tx=""; gene="";
                # if(match($0,/transcript_id "([^"]+)"/,a)) tx=a[1]
                # if(match($0,/gene_name "([^"]+)"/,b)) gene=b[1]
                # if(tx!="" && gene!="") print tx, gene
            # }' "$GTF_FILE" >> "$IDMAP"
        # echo "Created IDMAP: $IDMAP"
    # fi
# fi

# ###############################################################
# # STEP 2 — DaPars 3’UTR annotation
# ###############################################################
# THREEUTR="$OUTDIR/annotation/3UTR_annotation.bed"

# if [[ "$RESUME_MODE" == "wig" ]]; then
    # echo "[RESUME] Skipping 3’UTR annotation"
# else
    # if [[ -f "$THREEUTR" ]]; then
        # echo "3’UTR annotation exists — skipping."
    # else
        # echo "[2] Running DaPars_Extract_Anno.py..."
        # python $DAPARS/DaPars_Extract_Anno.py -b "$BED12" -s "$IDMAP" -o "$THREEUTR"
    # fi
# fi

# ###############################################################
# # STEP 3a — Submit WIG SLURM jobs
# ###############################################################
# if [[ "$RESUME_MODE" == "wig" ]]; then
    # echo "[RESUME] Skipping WIG generation."
# else
    # echo "[3a] Submitting WIG jobs…"

    # submitted_jobs=()

    # read header < "$SAMPLE_SHEET"
    # while IFS=',' read -r sample bam; do
        # [[ -z "$sample" ]] && continue

        # job="$sample"
        # slurm="$OUTDIR/slurm_logs/${job}.slurm"

        # cat > "$slurm" <<EOF
# #!/bin/bash
# #SBATCH --job-name=$job
# #SBATCH --cpus-per-task=$SLURM_CPUS
# #SBATCH --mem=$SLURM_MEM
# #SBATCH --time=$SLURM_TIME
# #SBATCH --output=$OUTDIR/slurm_logs/${job}.log

# bedtools genomecov -ibam $bam -bga -split -trackline \
    # > $OUTDIR/wig/${sample}.wig
# EOF

        # jid=$(sbatch "$slurm" | awk '{print $4}')
        # submitted_jobs+=("$jid")
        # echo " Submitted $job (JobID=$jid)"

    # done < "$SAMPLE_SHEET"

    # if [[ "$WAIT_MODE" == "manual" ]]; then
        # echo ""
        # echo "WIG jobs submitted."
        # echo "Check files in: $OUTDIR/wig/"
        # read -p "Press ENTER to continue once all WIG files exist..."
    # else
        # echo "[WAIT] Waiting for SLURM jobs..."
        # while true; do
            # running=0
            # for jid in "${submitted_jobs[@]}"; do
                # state=$(squeue -j "$jid" -h -o "%T" 2>/dev/null || echo "COMPLETED")
                # if [[ "$state" != "COMPLETED" && "$state" != "FAILED" ]]; then
                    # running=$((running+1))
                # fi
            # done
            # [[ $running -eq 0 ]] && break
            # echo " $running jobs running..."
            # sleep 30
        # done
    # fi

    # echo "[OK] WIG generation complete."
    # echo "Pipeline exiting. Re-run with: --resume wig"
    # exit 0
# fi

###############################################################
# STEP 3b — Optional Add 'chr' prefix to wig files
###############################################################
if [[ "$ADDCHR" == "yes" ]]; then
    echo "[3b] Adding 'chr' prefix to WIG files..."

    for wig in "$OUTDIR"/wig/*.wig; do
        {
            read header
            echo "$header"
            awk 'BEGIN{OFS="\t"} NR>1 { $1="chr"$1; print }'
        } < "$wig" > "${wig}.tmp"

        mv "${wig}.tmp" "$wig"
        echo " Updated: $(basename "$wig")"
    done

    echo "[OK] chr prefixes added."
else
    echo "[3b] Skipping chr-prefix modification (ADDCHR=$ADDCHR)"
fi

###############################################################
# STEP 4 — Build sequencing depth file
###############################################################

THREEUTR="$OUTDIR/annotation/3UTR_annotation.bed"

echo "[4] Extracting read depths from salmon logs..."
MAP_DEPTH="$OUTDIR/mapping_wig_location_with_depth.txt"
> "$MAP_DEPTH"

for log in $RNASEQ_OUT/*/logs/salmon_quant.log; do
    sample=$(basename "$(dirname "$(dirname "$log")")")
    depth=$(grep "# of uniquely mapped reads" "$log" | awk -F':' '{print $2}' | tr -d ' ')
    echo -e "$OUTDIR/wig/${sample}.wig\t$depth" >> "$MAP_DEPTH"
done

###############################################################
# STEP 5 — Build DaPars config
# Output_directory: prefix for the folder build under current folder
# Output_result_file: prefix for the output file 
###############################################################
CONFIG="$OUTDIR/DaPars_configure.txt"
: > "$CONFIG"
echo "[5] Writing DaPars2 config…"

wig_list=$(ls "$OUTDIR/wig/"*.wig | tr '\n' ',' | sed 's/,$//')

cat >> "$CONFIG" <<EOF
Annotated_3UTR=$THREEUTR
Aligned_Wig_files=$wig_list
Output_directory=$(basename "$OUTDIR")/dapars_output/dapars2
Output_result_file=dapars2
Coverage_threshold=10
Num_Threads=$THREADS
sequencing_depth_file=$MAP_DEPTH
EOF

###############################################################
# STEP 6 — Run DaPars2 per chromosome (1..22)
###############################################################
###############################################################
# STEP 6 — Submit DaPars2 per chromosome (SLURM, no wait)
###############################################################
echo "[6] Submitting DaPars2 SLURM jobs (no waiting)..."

for chr in {1..22}; do
    jobname="DaPars_chr${chr}"
    slurmfile="$OUTDIR/slurm_logs/dapars_chr${chr}.slurm"

    cat > "$slurmfile" <<EOF
#!/bin/bash
#SBATCH --job-name=$jobname
#SBATCH --cpus-per-task=$SLURM_CPUS
#SBATCH --mem=$SLURM_MEM
#SBATCH --time=$SLURM_TIME
#SBATCH --output=$OUTDIR/slurm_logs/dapars_chr${chr}.log

~/.conda/envs/polyfun/bin/python \
    "$DAPARS"/Dapars2_Multi_Sample.py \
    "$CONFIG" chr${chr}
EOF

    jid=$(sbatch "$slurmfile" | awk '{print $4}')
    echo " Submitted chr${chr} as JobID $jid"
done

echo "[6] All DaPars2 chromosome jobs submitted."
echo "======================================================"
echo "Finished analysis. APA results present in chromosome folders"
echo "======================================================"