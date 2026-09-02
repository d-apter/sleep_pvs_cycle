#!/bin/bash
# fMRIPrep anatomical (T1w) preprocessing + FreeSurfer recon-all.
#
# Requires: BIDS_DIR, OUT_DIR, WORK_DIR, FS_OUTPUT_DIR, FMRIPREP_SIF,
# FS_LICENSE set below.

SUBJECT=$(echo "$1" | sed 's/^sub-//')
SESSION=$(echo "$2" | sed 's/^ses-//')

BIDS_DIR="/path/to/bids"
OUT_DIR="${BIDS_DIR}/derivatives/fmriprep"
WORK_DIR="/tmp/fmriprep_work/sub-${SUBJECT}_ses-${SESSION}"
FS_OUTPUT_DIR="${BIDS_DIR}/derivatives/freesurfer_subjects"
FMRIPREP_SIF="/path/to/fmriprep_25.1.3.sif"
FS_LICENSE="/path/to/freesurfer/license.txt"

mkdir -p "$OUT_DIR" "$WORK_DIR" "$FS_OUTPUT_DIR"

# T1w-only BIDS filter
cat > "${WORK_DIR}/filter.json" << EOF
{
  "t1w": {
    "datatype": "anat",
    "session": "${SESSION}",
    "suffix": "T1w"
  }
}
EOF

singularity run \
    -B "${BIDS_DIR}":/bids:ro \
    -B "${OUT_DIR}":/out \
    -B "${FS_OUTPUT_DIR}":/fs_output \
    -B "${WORK_DIR}":/work \
    -B "${FS_LICENSE}":/opt/freesurfer/license.txt \
    "${FMRIPREP_SIF}" \
    /bids /out participant \
    --participant-label "${SUBJECT}" \
    --bids-filter-file /work/filter.json \
    --output-spaces MNI152NLin2009cAsym T1w \
    --fs-subjects-dir /fs_output \
    --work-dir /work \
    --anat-only \
    --nprocs 8 \
    --omp-nthreads 8 \
    --notrack \
    --skip_bids_validation
