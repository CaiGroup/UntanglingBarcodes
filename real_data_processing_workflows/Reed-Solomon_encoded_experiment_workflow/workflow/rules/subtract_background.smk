rule subtract_background:
    input:
        "resources/RS_seqfish_half_pool/HybCycle_{hyb}/HybCycle_{hyb}_MMStack_Pos{pos}.ome.tif",
        "resources/RS_seqfish_half_pool/chromatic_aberration/chromatic_aberration_MMStack_Pos{pos}.ome.tif",
        "results/alignment/offsets_ch{ch}_pos_{pos}.csv",
        "results/Shifted_Labeled_Images/hyb_{hyb}_ch{ch}_pos_{pos}.tif"
    params: r_med_filt = config['r_med_filt'], r_ball = config['rb_radius'], bg_sub_multiplier=config['bg_sub_multiplier'], bead_thresh=config['bead_thresh'], bead_dilation=config['bead_dilation_radius']
    output: "results/ims_bg_sub/HybCycle_{hyb}_ch{ch}_pos_{pos}.png"
    wildcard_constraints: hyb = "\d+"
    group: "fit_dots"
    resources: mem_mb = 5000
    conda: "../envs/bgsub.yaml"
    script: "../scripts/subtract_bg.py"

rule shift_label_images:
    input: "results/filtered_labeled_imgs/lbld_img_pos{pos}.png", #"results/labeled_stacks/labeled_stack_pos{pos}.tif",
           "results/alignment/offsets_ch{ch}_pos_{pos}.csv"
    output: "results/Shifted_Labeled_Images/hyb_{hyb}_ch{ch}_pos_{pos}.tif"
    conda: "../envs/bgsub.yaml"
    group: "fit_dots"
    resources: mem_mb = 3000
    script: "../scripts/shift_label_images.py"
