# Aggregates SHiVAi's per-region PVS volumes (pvs_stats.csv) into the
# region groupings used in the statistical analysis (CSO, basal ganglia,
# hippocampus, total).
#
# Requires:
#   tsubs    - list of subject IDs
#   dataset  - DataFrame with columns subject_id, session
#   root     - project path variable

import pandas as pd

for sub in tsubs:
    sessions = list(dataset[dataset.subject_id == sub]["session"])
    df = pd.DataFrame({"subject_id": sub, "session": sessions})

    for ses in sessions:
        shiva_path = f"{root}/derivatives/shiva/outputs/results/segmentations/pvs_segmentation/{sub}_{ses}"
        pvs = pd.read_csv(f"{shiva_path}/pvs_stats.csv")

        def region_vol(region):
            return pvs.loc[pvs.Region == region, "Total Biomarker volume"].item()

        pvs_cs = region_vol("Left-Cerebral-White-Matter") + region_vol("Right-Cerebral-White-Matter")
        df.loc[df.session == ses, "pvs_volume_centrumsemiovale"] = pvs_cs

        pvs_bg = (region_vol("Left-Caudate")  + region_vol("Right-Caudate") +
                  region_vol("Left-Putamen")   + region_vol("Right-Putamen") +
                  region_vol("Left-Pallidum")  + region_vol("Right-Pallidum"))
        df.loc[df.session == ses, "pvs_volume_basalganglia"] = pvs_bg

        pvs_hc = region_vol("Left-Hippocampus") + region_vol("Right-Hippocampus")
        df.loc[df.session == ses, "pvs_volume_hippocampus"] = pvs_hc

        total_pvs = pvs_cs + pvs_bg + pvs_hc
        df.loc[df.session == ses, "total_pvs_volume"] = total_pvs

    df.to_csv(f"{root}/derivatives/fmriprep/{sub}/{sub}_pvs-stat.csv")
