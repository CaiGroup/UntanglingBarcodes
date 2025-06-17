import pandas as pd 
import numpy as np 
import cv2
import skimage.io
import copy

ref_matches = pd.read_csv(snakemake.input[0])
dots = pd.read_csv(snakemake.input[1])


matches488 = pd.read_csv(snakemake.input[2]) #"results/alignment/matches_ch488_pos_0.csv")
matches561 = pd.read_csv(snakemake.input[3]) #"results/alignment/matches_ch561_pos_0.csv")
matches640 = pd.read_csv(snakemake.input[4]) #"results/alignment/matches_ch640_pos_0.csv")

matches488.loc[:, "hyb"] += 1
matches561.loc[:, "hyb"] += 1
matches640.loc[:, "hyb"] += 1

#matchesDict = {488 : matches488, 561 : matches561, 640 : matches640}
matchesDict = {2 : matches488, 1 : matches561, 0 : matches640}

ref_hyb_matches = pd.concat([matches488, matches561, matches640])

offsets488 = pd.read_csv(snakemake.input[5]) #"results/alignment/offsets_ch488_pos_0.csv")
offsets561 = pd.read_csv(snakemake.input[6]) #"results/alignment/offsets_ch561_pos_0.csv")
offsets640 = pd.read_csv(snakemake.input[7]) #"results/alignment/offsets_ch640_pos_0.csv")

#offsetsDict = {488 : offsets488, 561 : offsets561, 640 : offsets561}
offsetsDict = {2 : offsets488, 1 : offsets561, 0 : offsets561}


ch_n_matches = ref_matches.groupby("ref_ch").size()
ref_ch = ch_n_matches.index[ch_n_matches.argmax()]

ref_ch_offsets = offsetsDict[ref_ch]

dots['hyb'] = dots['hyb'] + 1
ref_ch_offsets['hyb'] = ref_ch_offsets['hyb'] + 1

dots.set_index(["ch", "hyb"], inplace = True)
ref_ch_offsets.set_index("hyb", inplace = True)

print("register reference channel (", ref_ch, ")")

for hyb in ref_ch_offsets.index:
    #dots.loc[hyb,"x":"y"] -= offsets.loc[hyb,'x_offset':'y_offset']
    if hyb in ref_ch_offsets.index and (ref_ch, hyb) in dots.index:
        dots.loc[(ref_ch, hyb),"x"] = dots.loc[(ref_ch, hyb),"x"] - ref_ch_offsets.loc[hyb,'x']
        dots.loc[(ref_ch, hyb),"y"] = dots.loc[(ref_ch, hyb),"y"] - ref_ch_offsets.loc[hyb,'y']
        #dots.loc[(ref_ch, hyb),"z"] = dots.loc[(ref_ch, hyb),"z"] - ref_ch_offsets.loc[hyb,'z']
    else:
        pass


ref_matches = ref_matches.loc[ref_matches.ref_ch == ref_ch]

for ch in matchesDict.keys():
    if ch != ref_ch:
        print("register ch: ", ch)
        merged = pd.merge(ref_matches.loc[ref_matches.comp_ch == ch], matchesDict[ch], right_on = ["ref_x", "ref_y", "ref_z"], left_on = ["comp_x", "comp_y", "comp_z"])
        for hyb in np.unique(merged.hyb_y):
            if (ch, hyb) in dots.index:
                hyb_merged = merged.loc[merged.hyb_y == hyb]
                tform, throwaway = cv2.estimateAffine2D(np.array(hyb_merged.loc[:, ["comp_x_y", "comp_y_y"]]), np.array(hyb_merged.loc[:, ["ref_x_x", "ref_y_x"]]))
                dots.loc[(ch, hyb), ["x", "y"]] = np.transpose(np.matmul(tform[:, :2],np.transpose(np.array(dots.loc[(ch, hyb), ["x", "y"]]))))
                dots.loc[(ch, hyb), "x"] += tform[0,2]
                dots.loc[(ch, hyb), "y"] += tform[1,2]

pc_round_table = pd.read_csv(snakemake.input[8])

#ch_ind_2_wavelength = {1 : 640, 2 : 561, 3 : 488}
ch_ind_2_wavelength = {1 : 0, 2 : 1, 3 : 2}

pc_dict = {(ch_ind_2_wavelength[row.channel], row["readout hyb"]) : row.pseudocolor for i, row in pc_round_table.iterrows()}
round_dict = {(ch_ind_2_wavelength[row.channel], row["readout hyb"]) : row["block"] for i, row in pc_round_table.iterrows()}

pcs = []
rounds = []
cellid = []
for i, row in dots.iterrows():
    if i in pc_dict:
        pcs.append(pc_dict[i])
        rounds.append(round_dict[i])
    else:
        pcs.append(None)
        rounds.append(None)

dots["pseudocolor"] = pcs
dots["block"] = rounds

labeled_img = skimage.io.imread(snakemake.input[9])
print("np.shape(labeled_img): ", np.shape(labeled_img))
dots = dots.loc[(dots.x < 2048) & (dots.y < 2048) & (dots.x > 0) & (dots.x > 0)]
#dots["cellid"] = [labeled_img[int(row.y)-1, int(row.x)-1] for i, row in dots.iterrows()]
dots["cellid"] = [labeled_img[int(row.y)-1, int(row.x)-1] for i, row in dots.iterrows()]



dots.to_csv(snakemake.output[0])

'''
cross_channel_register
params
    - ref_ch: (Int) the channel used as reference to which all channels are registered to
    - ref_ch_offsets: DataFame containing the x, y, and z translation offsets found between each readout
         image and the reference image in the reference channel. Indexed by hyb.
    - ref_matches: DataFrame containing the matches of dots in the reference images of each channel to each other
    - _dots_to_register: DataFrame columns include x, y ,z and indexed by channel and hyb
    - matchesDict: dictionary with keys for each channel, and value dataFrame of fiducial marker matches from each
        readout image in that channel to fiducial markers in that channel's reference image.
'''
def cross_channel_register_loov(ref_ch, ref_ch_offsets, ref_matches, _dots_to_register, matchesDict, ref_hyb_matches):
    dots_to_register = copy.deepcopy(_dots_to_register)
    for hyb in ref_ch_offsets.index:
        #dots.loc[hyb,"x":"y"] -= offsets.loc[hyb,'x_offset':'y_offset']
        if hyb in ref_ch_offsets.index and (ref_ch, hyb) in dots.index:
            dots_to_register.loc[(ref_ch, hyb),"x"] = dots_to_register.loc[(ref_ch, hyb),"x"] - ref_ch_offsets.loc[hyb,'x']
            dots_to_register.loc[(ref_ch, hyb),"y"] = dots_to_register.loc[(ref_ch, hyb),"y"] - ref_ch_offsets.loc[hyb,'y']
            #dots_to_register.loc[(ref_ch, hyb),"z"] = dots_to_register.loc[(ref_ch, hyb),"z"] - ref_ch_offsets.loc[hyb,'z']
        else:
            pass

    ref_ch_matches = ref_matches.loc[ref_matches.ref_ch == ref_ch]
    ref_ch_matches.reset_index(drop=True, inplace=True)
    loocv_alignments = pd.DataFrame(columns=["ref_x_refch","ref_y_refch","ref_z_refch", "hyb_targetch",
                                                "loocv_align_x", "loocv_align_y", "loocv_align_z",
                                               "x_diff", "y_diff", "z_diff", 'comp_ch'])
    for ch in matchesDict.keys():
        if ch != ref_ch:
            print("register ch: ", ch)
            merged = pd.merge(ref_ch_matches.loc[ref_ch_matches.comp_ch == ch], matchesDict[ch], left_on = ["comp_x", "comp_y", "comp_z"],
                              right_on = ["ref_x", "ref_y", "ref_z"], suffixes=('_refch', '_targetch'))
            
            print(merged.columns)
            for hyb in np.unique(merged.hyb_targetch):
                if (ch, hyb) in dots.index:
                    hyb_merged = merged.loc[merged.hyb_targetch == hyb]
                    tform, inliers = cv2.estimateAffine2D(np.array(hyb_merged.loc[:, ["comp_x_targetch", "comp_y_targetch"]]), np.array(hyb_merged.loc[:, ["ref_x_refch", "ref_y_refch"]]))
                    for fm_ref in np.unique(hyb_merged.index):
                        drop_hyb_merged = hyb_merged.drop(fm_ref)
                        
                        tform, inliers = cv2.estimateAffine2D(np.array(drop_hyb_merged.loc[:, ["comp_x_targetch", "comp_y_targetch"]]), np.array(drop_hyb_merged.loc[:, ["ref_x_refch", "ref_y_refch"]]), ransacReprojThreshold=np.inf)
                        assert(all(inliers == 1))
                        z_trans = np.mean(np.array(drop_hyb_merged.loc[:, "ref_z_refch"]) -np.array(drop_hyb_merged.loc[:, "comp_z_targetch"]))
                        #print("z_trans",z_trans)
                        #print(drop_hyb_merged.index)
                        #print(drop_hyb_merged.columns)
                        #np.round("sdfal")
                        loo_samp = copy.deepcopy(drop_hyb_merged.loc[:, ["ref_x_refch","ref_y_refch","ref_z_refch","comp_x_targetch", "comp_y_targetch", "comp_z_targetch", "hyb_targetch", "comp_ch"]])
                        loo_samp.loc[:, ["loocv_align_x", "loocv_align_y"]] = np.transpose(np.matmul(tform[:, :2],np.transpose(np.array(loo_samp.loc[:, ["comp_x_targetch", "comp_y_targetch"]]))))
                        loo_samp.loc[:, "loocv_align_x"] += tform[0,2]
                        loo_samp.loc[:, "loocv_align_y"] += tform[1,2]
                        loo_samp.loc[:, "loocv_align_z"] = z_trans + loo_samp.comp_z_targetch

                        #print(loo_samp)
                        loo_samp.loc[:, "x_diff"] = loo_samp["loocv_align_x"] - loo_samp["ref_x_refch"]
                        loo_samp.loc[:, "y_diff"] = loo_samp["loocv_align_y"] - loo_samp["ref_y_refch"]
                        loo_samp.loc[:, "z_diff"] = loo_samp["loocv_align_z"] - loo_samp["ref_z_refch"]
                        
                        loocv_alignments = pd.concat([loocv_alignments, loo_samp])#[["ref_x", "ref_y", "ref_z", "ref_int", "hyb",
                                                                #"loocv_align_x", "loocv_align_y", "loocv_align_z",
                                                                #"comp_int", "x_diff", "y_diff", "z_diff"]]])

    return loocv_alignments

cc_loov = cross_channel_register_loov(ref_ch, ref_ch_offsets, ref_matches, dots, matchesDict, ref_hyb_matches)

cc_loov.to_csv(snakemake.output[1])
