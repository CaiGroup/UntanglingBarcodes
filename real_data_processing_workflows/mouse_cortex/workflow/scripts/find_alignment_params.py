import pandas as pd
import copy

from seqfish_fm_match import FMMatcher
#from FMAligner import FMAligner

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

dal = FMMatcher(ro, ref_init, ref_final)
opt_matches, mm_lt, mm_zt, x_offset, y_offset, z_offset, lat_offset = dal.auto_set_params()
dal.save_params(snakemake.output[0])#"params.csv")
