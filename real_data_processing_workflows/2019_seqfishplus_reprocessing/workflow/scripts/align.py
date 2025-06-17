import pandas as pd
import copy

from FMAligner import FMAligner

# Datafiles are not in repository because of space
ref_init = pd.read_csv(snakemake.input[0])
ref_final = pd.read_csv(snakemake.input[1])
ro = pd.read_csv(snakemake.input[2])

ref_final.drop("s", axis=1, inplace = True)
ref_init.drop("s", axis=1, inplace = True)
ro.drop("s", axis=1, inplace = True)

ref_final.insert(2, 'z', 0)
ref_init.insert(2, 'z', 0)
ro.insert(3,'z', 0)

ref_final.rename(columns={'w':'int'}, inplace=True)
ref_init.rename(columns={'w':'int'}, inplace=True)
ro.rename(columns={'w':'int'}, inplace=True)


ref_final_ro_form = copy.deepcopy(ref_final)
ref_final_ro_form["hyb"] = 1
ro.set_index("hyb", inplace=True)

dal = FMAligner(ro, ref_init, ref_final)
dal.load_saved_parameters(snakemake.input[3])#"params.csv")

# set manually provided parameters
print(snakemake.params['min_fm_hyb_matches'])
if snakemake.params['min_fm_hyb_matches'] != 'default':
    dal.set_min_fm_hyb_matches(snakemake.params['min_fm_hyb_matches'])

print(snakemake.params['outlier_sd_thresh'])
if snakemake.params['outlier_sd_thresh'] != 'default':
    dal.set_outlier_sd_thresh(snakemake.params['outlier_sd_thresh'])

print(snakemake.params['max_lat_offset'])
if snakemake.params['max_lat_offset'] != "none":
    dal.set_max_lat_offset(snakemake.params['max_lat_offset'])

print(snakemake.params['set_xy_search_error'])
if snakemake.params['set_xy_search_error']:
    dal.set_xy_search_error(snakemake.params['set_xy_search_error'])


dal.set_max_lat_offset(50)

dal.align()
dal.save_offsets(snakemake.output[0])#"offsets_all.csv")
dal.save_matches(snakemake.output[1])#"matches_all.csv")
dal.save_loocv_errors(snakemake.output[2])#"loov_errors_all.csv")
#dal.save_ro_wout_fm(snakemake.output[3])#"pnts_no_fm_all.csv")
