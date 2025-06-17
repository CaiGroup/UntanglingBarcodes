rule fit_fiducials_ref:
    input: "results/cell_masked_bead_stacks/chromatic_aberration/stack_ch{ch}_pos{pos}.png"
    params:
        sigma_xy_lb = config["fit_fiducials"]["sigma_lb"],
        sigma_xy_ub = config["fit_fiducials"]["sigma_ub"],
        sigma_z_lb = config["fit_fiducials"]["sigma_z_lb"],
        sigma_z_ub = config["fit_fiducials"]["sigma_z_ub"],
        final_loss_improvement = config["fit_fiducials"]["final_loss_improvement"],
        min_weight = config["fit_fiducials"]["min_weight"],
        max_iters = config["fit_fiducials"]["max_iters"],
        max_cd_iters = config["max_cd_iters"],
        min_allowed_separation = config["min_allowed_separation"],
        tile_main_width = config["tile_main_width"],
        tile_overlap = config["tile_overlap"],
        tile_depth = config["tile_depth"],
        tile_depth_overhang = config["tile_depth_overhang"]
    output:
        "results/fm_fits/chromatic_aberration_ch{ch}_pos{pos}.csv"#,
    resources: mem_mb=5000
    group: "fit_dots"
    script: "../scripts/tile_fit_3d_fiducials.jl"

rule fit_fiducials:
    input: "results/cell_masked_bead_stacks/HybCycle_{hyb}/masked_stack_ch{ch}_pos{pos}.png"
    params:
        sigma_xy_lb = config["fit_fiducials"]["sigma_lb"],
        sigma_xy_ub = config["fit_fiducials"]["sigma_ub"],
        sigma_z_lb = config["fit_fiducials"]["sigma_z_lb"],
        sigma_z_ub = config["fit_fiducials"]["sigma_z_ub"],
        final_loss_improvement = config["fit_fiducials"]["final_loss_improvement"],
        min_weight = config["fit_fiducials"]["min_weight"],
        max_iters = config["fit_fiducials"]["max_iters"],
        max_cd_iters = config["max_cd_iters"],
        min_allowed_separation = config["min_allowed_separation"],
        tile_main_width = config["tile_main_width"],
        tile_overlap = config["tile_overlap"],
        tile_depth = config["tile_depth"],
        tile_depth_overhang = config["tile_depth_overhang"]
    output:
        "results/fm_fits/HybCycle_{hyb}_ch{ch}_pos{pos}.csv"#,
    resources: mem_mb=5000
    group: "fit_dots"
    script: "../scripts/tile_fit_3d_fiducials.jl"


rule aggregate_fm_fits:
    #input: expand("fit_dots/{hyb}.csv", hyb = hybs)
    input: expand("results/fm_fits/HybCycle_{hyb}_ch{ch}_pos{{pos}}.csv", hyb = all_hybs, ch = config["channels"])
    output: "results/fm_fits/all_fms_pos_{pos}.csv"
    conda: "../envs/pdnp.yaml"
    resources: mem_mb=500
    script: "../scripts/aggregate_dots.py"

rule ref_cvt_to_regular_tiff:
    input: 
        "resources/RS_seqfish_half_pool/chromatic_aberration/chromatic_aberration_MMStack_Pos{pos}.ome.tif",
        "results/filtered_labeled_imgs/lbld_img_pos{pos}.png" #"results/labeled_stacks/labeled_stack_pos{pos}.tif"
    output: "results/cell_masked_bead_stacks/chromatic_aberration/stack_ch{ch}_pos{pos}.png"
    script: "../scripts/cvt_to_regular_tiff.py"

rule cvt_to_regular_tiff:
    input: 
        "resources/RS_seqfish_half_pool/HybCycle_{hyb}/HybCycle_{hyb}_MMStack_Pos{pos}.ome.tif",
        "results/filtered_labeled_imgs/lbld_img_pos{pos}.png" #"results/labeled_stacks/labeled_stack_pos{pos}.tif"
    output: "results/cell_masked_bead_stacks/HybCycle_{hyb}/masked_stack_ch{ch}_pos{pos}.png"
    resources: mem_mb=2000, disk_mb=5000
    script: "../scripts/get_masked_bead_stacks.py"
