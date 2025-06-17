import pandas as pd
#from FMAligner import FMAligner
from seqfish_fmalign import FMAligner
from scipy.spatial import KDTree

ch_fids_ = [pd.read_csv(file) for file in snakemake.input]

ch_fids = []
ind_ch = {}
for i, fids in enumerate(ch_fids_):
    fids["int"] = fids["w"]
    fids["hyb"] = i
    tree = KDTree(fids.loc[:,["x", "y"]])
    fids = fids.loc[[len(tree.query_ball_point(row[["x","y"]], 4)) < 2 for i, row in fids.iterrows()]]
    ch_fids.append(fids[["hyb", "x", "y", "z", "int", "ch"]])
    ind_ch[i] = fids.ch.iloc[0]

ros = pd.concat(ch_fids)
#ros["hyb"] = ros["ch"]
ros.set_index("hyb", inplace=True)

matches = []
for i, fids in enumerate(ch_fids):
    aligner = FMAligner(ros, fids)

    aligner.set_xy_search_error(10.0)
    aligner.set_z_search_error(1.0)
    aligner.set_min_bright_prop(0.01)
    aligner.set_max_bright_prop(100)
    aligner.set_max_lat_offset(50)

    aligner.align()

    print(aligner.matchesDF)

    fid_matches = aligner.matchesDF
    fid_matches["ref_ch"] = fids.loc[0,"ch"]
    #fid_matches["comp_ch"] = fid_matches["hyb"]

    matches.append(fid_matches)

all_matches = pd.concat(matches)

all_matches["comp_ch"] = [ind_ch[h] for h in all_matches.hyb]

#print()
#print("len(all_matches): ", len(all_matches))
all_matches = all_matches.loc[all_matches.ref_ch != all_matches.comp_ch]
#print("len(all_matches): ", len(all_matches))


all_matches.to_csv(snakemake.output[0])
