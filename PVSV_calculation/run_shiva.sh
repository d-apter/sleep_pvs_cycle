#!/bin/bash
# Runs the SHiVAi PVS segmentation pipeline via Apptainer/Singularity.
# Requires: SHiVAi container (shiva_0.4.2.sif), a patched image.py interface
# (see patched_image.py), FreeSurfer LUT, and a prepared config file.
# Expects derivatives/shiva/inputs/ to already contain the T1w inputs in
# SHiVAi's "standard" input structure (see SHiVAi README).

singularity exec --nv \
  -B $(pwd) \
  -B code/shiva/patched_image.py:/usr/local/lib/python3.8/dist-packages/shivai/interfaces/image.py \
  code/containers/shiva_0.4.2.sif \
  shiva \
    --in derivatives/shiva/inputs \
    --out derivatives/shiva/outputs \
    --config code/shiva/shiva_config.yml \
    --prediction PVS \
    --input_type standard \
    --brain_seg custom \
    --custom_LUT code/shiva/FreeSurferColorLUT.txt
