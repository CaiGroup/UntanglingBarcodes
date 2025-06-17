rule subtract_background:
    input:
        "results/rep{rep}/ch_pngs/HybCycle_{hyb}_ch_{ch}_pos_{pos}.png",
        "results/rep{rep}/ch_pngs/beads_initial_ch_{ch}_pos_{pos}.png",
        "results/rep{rep}/alignment/offsets_ch_{ch}_pos_{pos}.csv",
        "results/rep{rep}/Shifted_Labeled_Images/hyb_{hyb}_ch_{ch}_pos_{pos}.png"
    params: r_med_filt = config['r_med_filt'], r_ball = config['rb_radius']
    output: "results/rep{rep}/ims_bg_sub/HybCycle_{hyb}_ch_{ch}_pos_{pos}.png"
    wildcard_constraints: hyb = "\d+"
    group: "fit_dots"
    resources: mem_mb = 5000
    conda: "../envs/bgsub.yaml"
    script: "../scripts/subtract_bg.py"

rule convert_contours_to_labeled_ims:
    input: "resources/replicates/roi/RoiRep{rep}/RoiSet_Pos{pos}.zip"
    output: "results/rep{rep}/Labeled_Images/MMStack_Pos{pos}.png"
    resources: mem_mb = 3000
    conda: "../envs/contour2labeledim.yaml"
    script: "../scripts/roi_contours_2_labeledim.py"

rule shift_label_images:
    input: "results/rep{rep}/Labeled_Images/MMStack_Pos{pos}.png",
           "results/rep{rep}/alignment/offsets_ch_{ch}_pos_{pos}.csv"#"dapi_align/offsets_hyb_{hyb}_pos_{pos}.csv"
    output: "results/rep{rep}/Shifted_Labeled_Images/hyb_{hyb}_ch_{ch}_pos_{pos}.png"
    conda: "../envs/bgsub.yaml"
    group: "fit_dots"
    resources: mem_mb = 3000
    script: "../scripts/shift_label_images.py"
