# Warps each subject/session's FreeSurfer aparc+aseg parcellation (native
# space) into MNI space, to serve as SHiVAi's --custom_LUT brain
# segmentation for region-wise PVS statistics.
#
# Requires:
#   subjects              - list of subject IDs
#   dataset               - DataFrame with columns subject_id, session
#   get_fmriprep_data()   - returns path to an fMRIPrep output file
#   get_freesurfer_path() - returns path to a subject/session's FreeSurfer dir
#   root, work, output    - project path variables
#   fmriprep               - fMRIPrep container image path
#   freesurfer_license     - path to the local FreeSurfer license.txt
#   singularity            - Singularity/Apptainer executable
#   submit()               - cluster job submission wrapper

import os

for sub in subjects:
    sessions = list(dataset[dataset.subject_id == sub]["session"])
    for ses in sessions:
        ref_path = get_fmriprep_data(sub, ses, 'desc-preproc_T1w.nii.gz')
        freesurfer_path = get_freesurfer_path(sub, ses)

        fsnative_to_t1 = f"{root}/derivatives/fmriprep/{sub}/{ses}/anat/{sub}_{ses}_from-fsnative_to-T1w_mode-image_xfm.txt"
        t1_to_mni = f"{root}/derivatives/fmriprep/{sub}/{ses}/anat/{sub}_{ses}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5"

        out_path = ref_path.replace("space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz", "aparc+aseg.nii.gz")

        # Single-step resampling: reads native-space .mgz, applies both
        # transforms, writes MNI-space .nii.gz directly.
        ants_command = (
            "antsApplyTransforms "
            "-d 3 "
            f"-i {freesurfer_path}/mri/aparc+aseg.mgz "
            f"-r {ref_path} "
            f"-t {t1_to_mni} "
            f"-t {fsnative_to_t1} "
            "-n NearestNeighbor "
            f"-o {out_path} "
            "--verbose"
        )

        fmriprep_run = (
            "exec "
            f"-B {work},{root},{output},{root}/code:/code,{freesurfer_license}:/opt/freesurfer/license.txt "
            f"{fmriprep} "
            f"{ants_command} "
        )

        if not os.path.exists(out_path):
            submit("fsnative_to_mni", singularity, fmriprep_run)
