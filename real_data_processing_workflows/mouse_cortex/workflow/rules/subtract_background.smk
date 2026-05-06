rule subtract_background:
    input:
        "resources/imgs/pos{pos}/HybCycle_{hyb}_pos_{pos}_z_{z}.png", #"resources/imgs/Cy3B_Pos{pos}_png/pos{pos}_hyb_{hyb}.png",
        "resources/imgs/pos{pos}/initial_background_MMStack_Pos{pos}_z_{z}.png", #"resources/imgs/Cy3B_Pos{pos}_png/init_bg_pos{pos}.png",
        "results/alignment/offsets_pos_{pos}.csv",
        "results/Shifted_Labeled_Images/hyb_{hyb}_pos_{pos}_z_{z}.png"
    params: r_med_filt = config['r_med_filt'], r_ball = config['rb_radius'],bead_thresh=config['bead_thresh'], bead_dilation_radius=config['bead_dilation_radius']
    output: "results/ims_bg_sub/HybCycle_{hyb}_pos_{pos}_z_{z}.png"
    benchmark: "benchmarks/ims_bg_sub/HybCycle_{hyb}_pos_{pos}_z_{z}.txt"
    wildcard_constraints: hyb = "\d+", z= "\d+", pos = "\d+"
    group: "fit_dots"
    resources: mem_mb = 5000
    conda: "../envs/bgsub.yaml"
    script: "../scripts/subtract_bg.py"

rule shift_label_images:
    input: "results/filtered_labeled_imgs/lbld_img_pos{pos}_z_{z}.png", #"results/labeled_stacks/labeled_stack_pos{pos}.tif",
           "results/alignment/offsets_pos_{pos}.csv"
    output: "results/Shifted_Labeled_Images/hyb_{hyb}_pos_{pos}_z_{z}.png"
    benchmark: "benchmarks/Shifted_Labeled_Images/hyb_{hyb}_pos_{pos}_z_{z}.txt"
    conda: "../envs/bgsub.yaml"
    wildcard_constraints: hyb = "\d+", z= "\d+", pos = "\d+"
    group: "fit_dots"
    #resources: mem_mb = 3000
    script: "../scripts/shift_label_images.py"


rule filter_labeled_regions:
    input: "resources/labeled_stacks/segmentation_MMStack_Pos{pos}.tif"
    output: expand("results/filtered_labeled_imgs/lbld_img_pos{{pos}}_z_{z}.png", z=list(range(11)))
    benchmark: "benchmarks/filtered_labeled_imgs/lbld_img_pos{pos}.txt"
    script: "../scripts/filter_labeled_stacks.py"
