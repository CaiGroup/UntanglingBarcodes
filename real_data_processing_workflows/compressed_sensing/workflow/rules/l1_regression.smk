import skimage
import numpy as np

rule find_cell_bounding_boxes:
    input: "resources/replicates/rep{rep}/Labeled_Images/MMStack_Pos{pos}.png" #labeled_img.tif"
    group: "find_dots"
    output: "results/cand_dots/cell_bounding_boxes_rep{rep}_pos{pos}.csv"
    script: "../scripts/get_cell_bounding_boxes.jl"


for rep in config['positions']:
    for pos in config["positions"][rep]:
        #rule split_cells
        lbld_img = skimage.io.imread("resources/replicates/" +rep +f"/Labeled_Images/MMStack_Pos{pos}.png") #.format(pos=pos))
        np.unique(lbld_img)
        cells = [cellID for cellID in np.unique(lbld_img) if cellID != 0]


        rule: #rule split_cell_imgss:
            input: 
                f"results/{rep}/ims_bg_sub/HybCycle_{{hyb}}_ch_{{ch}}_pos_{pos}.png", #.format(pos=pos, rep=rep),
                f"results/cand_dots/cell_bounding_boxes_{rep}_pos{pos}.csv",
                f"resources/replicates/{rep}/Labeled_Images/MMStack_Pos{pos}.png",
                f"results/{rep}/alignment/offsets_ch_{{ch}}_pos_{pos}.csv"
            group: "find_dots"
            output: expand(f"results/cell_imgs/{rep}_" +"HybCycle_{{hyb}}_ch_{{ch}}" + f"_pos_{pos}" + "_cell_{cell}.tif", cell=cells)
            wildcard_constraints: cell="\d+"
            script: "../scripts/split_cells.py"
  
        rule: #agg best curve cells
            input: expand(f"results/plots/est_fdr_vs_sens_best_sum_{rep}_pos_{pos}_"+"ch_{{ch}}_cell_{cell}.csv", cell=cells)
            output: f"results/plots/est_fdr_vs_sens_best_sum_{rep}_pos_{pos}_ch_"+"{ch}.csv"
            wildcard_constraints: cell="\d+"
            resources: mem_mb = 4000
            script: "../scripts/aggregate_sum_stats.py"

        #rule: # aggregate_cand_dots:
        #    input: expand("results/" + f"{rep}/cand_dots/cand_dots_{rep}" + "_HybCycle_{hyb}_ch_{ch}_" + f"pos_{pos}" + "_cell_{{cell}}.parquet", hyb=hybs, ch =config["channels"]) #expand("results/cand_dots/cand_dots_cell_{{cell}}_r{r}_pc{pc}.parquet", r = config["blocks"], pc = config["pseudocolors"])
        #    output: f"results/{rep}/cand_dots/cand_dots_pos_{pos}" + "_cell_{cell}.parquet"
        #    wildcard_constraints: cell="\d+"
        #    script: "../scripts/aggregate_cand_dots.jl"

        for ch in config["channels"]:
             for mpcd in config["mpcd"][ch]:
                rule:
                    input: expand(f"results/{rep}/lassoed_cpaths/lasso_path_summary_" + f"pos_{pos}_ch_{ch}" + "_cell_{cell}_lf{lf}_sr{sr}" + f"_mpcd{mpcd}.csv", cell=cells, lf=config["lfs"][ch], sr=config["search_radii"][ch])
                    output: f"results/summaries_agged/summaries_{rep}_pos_{pos}_ch_{ch}.csv"
                    wildcard_constraints: cell="\d+"
                    script: "../scripts/aggregate_sum_stats.py"
    
    #rule:
    #    input: expand("results/summaries_agged/" + f"summaries_{rep}" + "_pos_{pos}_ch_{ch}.csv", pos=config["positions"][rep], ch=config["channels"])
    #    output: f"results/summaries_agged/summaries_{rep}_all_pos_ch.csv"
    #    script: "../scripts/aggregate_sum_stats.py"

    rule: 
        input: expand(f"results/plots/est_fdr_vs_sens_best_sum_{rep}_"+"pos_{pos}_ch_{{ch}}.csv", pos=config["positions"][rep]) #, ch=config["channels"])
        output: f"results/plots/est_fdr_vs_sens_best_sum_{rep}_all_pos_"+"ch_{ch}.csv"
        resources: mem_mb = 4000
        script: "../scripts/aggregate_sum_stats.py"    

#rule aggregate_reps:
#    input: expand("results/summaries_agged/summaries_{rep}_all_pos_ch.csv", rep=config['positions'])
#    output: "results/summaries_agged/summaries_all_reps_all_pos_ch.csv"
#    script: "../scripts/aggregate_sum_stats.py"

rule aggregate_reps_best_curve_ch:
    input: expand("results/plots/est_fdr_vs_sens_best_sum_rep{rep}_all_pos_ch_{{ch}}.csv", rep=config['reps']) #, ch=config["channels"])
    output: "results/plots/est_fdr_vs_sens_best_sum_all_reps_all_pos_ch_{ch}.csv"
    resources: mem_mb = 4000
    script: "../scripts/aggregate_sum_stats.py"


rule aggregate_reps_best_curve:
    input: expand("results/plots/est_fdr_vs_sens_best_sum_all_reps_all_pos_ch_{ch}.csv", ch=config["channels"])
    output: "results/plots/est_fdr_vs_sens_best_sum_all_reps_all_pos_ch.csv"
    resources: mem_mb = 4000
    script: "../scripts/aggregate_sum_stats.py"

rule plot_best_curves_all_reps_ch:
    input: 
        bc488 = "results/plots/est_fdr_vs_sens_best_sum_all_reps_all_pos_ch_488.csv",
        bc561 = "results/plots/est_fdr_vs_sens_best_sum_all_reps_all_pos_ch_561.csv",
        bc643 = "results/plots/est_fdr_vs_sens_best_sum_all_reps_all_pos_ch_643.csv",
        cb488 = "resources/codebooks/channel_488.csv",
        cb561 = "resources/codebooks/channel_561.csv",
        cb643 = "resources/codebooks/channel_643.csv",
        smFISH_ref = "resources/validation/smFISH_results.csv"
    output: "results/plots/best_eFDR_sensitivity_overall.png", "results/plots/best_eFDR_sensitivity_overall.csv", "results/plots/overall_eFDR_0.5_sensitivity_fit.png"
    resources: mem_mb = 4000
    script: "../scripts/plot_best_eFDR_sensitivity_overall.py"

rule find_cand_dots:
    input: 
        "results/cell_imgs/rep{rep}_HybCycle_{hyb}_ch_{ch}_pos_{pos}_cell_{cell}.tif", #, "resources/labeled_img.tif", "results/cand_dots/cell_bounding_boxes.csv"
        env="workflow/envs/julia_environment",
        julia_env_installed="results/julia_environment_installed.txt"
    group: "find_dots"
    output: "results/rep{rep}/cand_dots/cand_dots_HybCycle_{hyb}_ch_{ch}_pos_{pos}_cell_{cell}_mpcd{mpcd}.parquet"
    params: sigma = config["psf_sigma"], min_peak=config["min_peak"], psf_cutoff = config["psf_cutoff"]
    benchmark: "benchmarks/rep{rep}_cand_dots_HybCycle_{hyb}_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.txt"
    resources: mem_mb = 4000
    script: "../scripts/find_cand_dots_sparse_mod.jl"

rule aggregate_cand_dots:
    input: expand("results/rep{{rep}}/cand_dots/cand_dots_HybCycle_{hyb}_ch_{{ch}}_pos_{{pos}}_cell_{{cell}}_mpcd{{mpcd}}.parquet", hyb=hybs) #, ch =config["channels"]) #expand("results/cand_dots/cand_dots_cell_{{cell}}_r{r}_pc{pc}.parquet", r = config["blocks"], pc = config["pseudocolors"])
    group: "find_dots"
    output: "results/rep{rep}/cand_dots/cand_dots_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.parquet"
    wildcard_constraints: cell="\d+"
    conda: "../envs/align.yaml"
    script: "../scripts/aggregate_cand_dots.jl"

rule apply_alignment_cand_dots:
    input: "results/rep{rep}/alignment/offsets_ch_{ch}_pos_{pos}.csv",
           "results/rep{rep}/cand_dots/cand_dots_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.parquet"
    group: "find_dots"
    output: "results/rep{rep}/cand_dots/aligned_cand_dots_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.parquet"
    wildcard_constraints: cell="\d+"
    script: "../scripts/apply_alignments_cand_dots.py"

rule find_cand_cpaths:
    input:
        "results/full_codebooks/codebook_ch_{ch}.csv",
        "resources/codebooks/H.txt",
        "results/rep{rep}/cand_dots/aligned_cand_dots_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.parquet",
        env="workflow/envs/julia_environment",
        julia_env_installed="results/julia_environment_installed.txt"
    output: "results/rep{rep}/cand_cpaths/cand_cpaths_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.csv"
    params: lat_search_radius = max(config["search_radii"]), z_search_radius = config["z_search_radius"]
    benchmark: "benchmarks/rep{rep}_cand_cpaths_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.txt"
    resources: mem_mb = 2000
    wildcard_constraints: cell="\d+"
    script: "../scripts/find_cand_cpaths.jl"

rule build_lasso_model:
    input:
        cand_cpaths="results/rep{rep}/cand_cpaths/cand_cpaths_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.csv",
        cand_dots_unregistered="results/rep{rep}/cand_dots/cand_dots_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.parquet",
        cand_dots_aligned="results/rep{rep}/cand_dots/aligned_cand_dots_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.parquet",
        imgs = expand("results/cell_imgs/rep{{rep}}_HybCycle_{hyb}_ch_{{ch}}_pos_{{pos}}_cell_{{cell}}.tif", hyb=hybs), #, ch=config["channels"])
        env="workflow/envs/julia_environment",
        julia_env_installed="results/julia_environment_installed.txt"
    output:
        "results/rep{rep}/lasso_mods/y_wnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.parquet",
        "results/rep{rep}/lasso_mods/A_wnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.parquet",
        "results/rep{rep}/lasso_mods/cpaths_wnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.csv",
        "results/rep{rep}/lasso_mods/y_nnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.parquet",
        "results/rep{rep}/lasso_mods/A_nnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.parquet",
        "results/rep{rep}/lasso_mods/cpaths_nnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.csv"
    params: psf_cutoff = config["psf_cutoff"], sigma = config["psf_sigma"], min_peak=config["min_peak"], search_radius="{sr}"
    benchmark: "benchmarks/rep{rep}_build_lasso_mod_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.txt"
    log: out = "logs/rep{rep}_lasso_mods_nnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.log",
         err = "logs/rep{rep}_lasso_mods_nnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.err"
    resources: mem_mb = 80000
    wildcard_constraints: cell="\d+"
    script: "../scripts/build_lasso_mod.jl"

rule run_lasso_mod:
    input:
        cpaths_wnc = "results/rep{rep}/lasso_mods/cpaths_wnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.csv",
        cpaths_nnc = "results/rep{rep}/lasso_mods/cpaths_nnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.csv",
        cdots = "results/rep{rep}/cand_dots/cand_dots_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.parquet",
        cdots_aligned = "results/rep{rep}/cand_dots/aligned_cand_dots_pos_{pos}_ch_{ch}_cell_{cell}_mpcd{mpcd}.parquet",
        y_wnc = "results/rep{rep}/lasso_mods/y_wnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.parquet",
        y_nnc = "results/rep{rep}/lasso_mods/y_nnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.parquet",
        A_wnc = "results/rep{rep}/lasso_mods/A_wnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.parquet",
        A_nnc = "results/rep{rep}/lasso_mods/A_nnc_pos_{pos}_ch_{ch}_cell_{cell}_sr{sr}_mpcd{mpcd}.parquet",
        env="workflow/envs/julia_environment",
        julia_env_installed="results/julia_environment_installed.txt"
    output:
        "results/rep{rep}/lassoed_cpaths/lassoed_cpaths_nnc_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}_mpcd{mpcd}.csv",
        "results/rep{rep}/lassoed_cpaths/lassoed_cpaths_wnc_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}_mpcd{mpcd}.csv",
        "results/rep{rep}/lassoed_cpaths/lasso_path_summary_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}_mpcd{mpcd}.csv",
        "results/rep{rep}/lassoed_cpaths/betas_nnc_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}_mpcd{mpcd}.png",
        "results/rep{rep}/lassoed_cpaths/betas_wnc_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}_mpcd{mpcd}.png",
    params: psf_cutoff = config["psf_cutoff"], sigma = config["psf_sigma"], min_peak=config["min_peak"]
    benchmark: "benchmarks/rep{rep}_lassoed_cpaths_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}_mpcd{mpcd}.txt"
    resources: mem_mb = 50000
    wildcard_constraints: cell="\d+"
    script: "../scripts/run_lasso_mod.jl"



"""
rule lasso_cpaths:
    input:
        "results/rep{rep}/cand_cpaths_pos_{pos}_ch_{ch}_cell_{cell}.csv",
        "results/rep{rep}/cand_dots/cand_dots_pos_{pos}_ch_{ch}_cell_{cell}.parquet",
        imgs = expand("results/cell_imgs/rep{{rep}}_HybCycle_{hyb}_ch_{{ch}}_pos_{{pos}}_cell_{{cell}}.tif", hyb=hybs) #, ch=config["channels"])
    output:
        "results/rep{rep}/lassoed_cpaths/lassoed_cpaths_nnc_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}.csv",
        "results/rep{rep}/lassoed_cpaths/lassoed_cpaths_wnc_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}.csv",
        "results/rep{rep}/lassoed_cpaths/lasso_path_summary_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}.csv",
        "results/rep{rep}/lassoed_cpaths/betas_nnc_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}.png",
        "results/rep{rep}/lassoed_cpaths/betas_wnc_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}.png",
    params: psf_cutoff = config["psf_cutoff"], sigma = config["psf_sigma"], min_peak=config["min_peak"], search_radius="{sr}"
    benchmark: "benchmarks/rep{rep}_lassoed_cpaths_pos_{pos}_ch_{ch}_cell_{cell}_lf{lf}_sr{sr}.txt"
    script: "../scripts/lasso_cpaths_penalize_pos_var.jl"

rule plot_lf_comparison_ch488:
    input:
        nnc_files = expand("results/rep{rep}/lassoed_cpaths/lassoed_cpaths_nnc_pos_{{pos}}_ch_{{ch}}_cell_{{cell}}_lf{lf}_sr{sr}.csv", lf=config["lfs"], sr=config["search_radii"]),
        nnc_beta_paths = expand("results/rep{rep}/lassoed_cpaths/betas_nnc_pos_{{pos}}_ch_{{ch}}_cell_{{cell}}_lf{lf}_sr{sr}.png", lf=config["lfs"], sr=config["search_radii"]),
        wnc_files = expand("results/rep{rep}/lassoed_cpaths/lassoed_cpaths_wnc_pos_{{pos}}_ch_{{ch}}_cell_{{cell}}_lf{lf}_sr{sr}.csv", lf=config["lfs"], sr=config["search_radii"]),
        wnc_beta_paths = expand("results/rep{rep}/lassoed_cpaths/betas_wnc_pos_{{pos}}_ch_{{ch}}_cell_{{cell}}_lf{lf}_sr{sr}.png", lf=config["lfs"], sr=config["search_radii"]),
    output:
        "results/plots/est_fdr_vs_ge_rep{rep}_pos_{pos}_ch_{ch}_cell_{cell}.png"
    params:
        lfs = config["lfs"],
        search_radii = config["search_radii"],
        beta_thresholds = config["beta_thresholds"]
    script:
        "../scripts/plot_cell_lf_sr_bt_curves.py"

"""


for ch in config["channels"]:
    """
    rule: # find_cand_cpaths:
        input:
            f"results/full_codebooks/codebook_ch_{ch}.csv",
            "resources/codebooks/H.txt",
            f"results/rep{{rep}}/cand_dots/aligned_cand_dots_pos_{{pos}}_ch_{ch}_cell_{{cell}}_mpcd{{mpcd}}.parquet",
            env="workflow/envs/julia_environment",
            julia_env_installed="results/julia_environment_installed.txt"
        output: f"results/rep{{rep}}/cand_cpaths/cand_cpaths_pos_{{pos}}_ch_{ch}_cell_{{cell}}_mpcd{{mpcd}}.csv"
        params: lat_search_radius = max(config["search_radii"][ch]), z_search_radius = config["z_search_radius"]
        benchmark: f"benchmarks/rep{{rep}}_cand_cpaths_pos_{{pos}}_ch_{ch}_cell_{{cell}}_mpcd{{mpcd}}.txt"
        resources: mem_mb = 2000
        script: "../scripts/find_cand_cpaths.jl"
        """


    for mpcd in config["mpcd"][ch]:
        rule: # plot_lf_comparison_chx:
            input:
                nnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_nnc_pos_{{pos}}_" + f"ch_{ch}"+"_cell_{{cell}}_lf{lf}_sr{sr}" + f"_mpcd{mpcd}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
                nnc_beta_paths = expand("results/rep{{rep}}/lassoed_cpaths/betas_nnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}" + f"_mpcd{mpcd}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
                wnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}" + f"_mpcd{mpcd}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
                wnc_beta_paths = expand("results/rep{{rep}}/lassoed_cpaths/betas_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}" + f"_mpcd{mpcd}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
            output:
                f"results/plots/est_fdr_vs_ge_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}_mpcd{mpcd}.png"
            params:
                lfs = config["lfs"][ch],
                search_radii = config["search_radii"][ch],
                beta_thresholds = config["beta_thresholds"][ch]
            script:
                "../scripts/plot_cell_lf_sr_bt_curves.py"

    rule: # plot_best_curves_chx:
        input:
            nnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_nnc_pos_{{pos}}_"+f"ch_{ch}"+"_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            nnc_betas = expand("results/rep{{rep}}/lassoed_cpaths/betas_nnc_pos_{{pos}}_"+f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            wnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            wnc_betas = expand("results/rep{{rep}}/lassoed_cpaths/betas_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            paths_summaries = expand("results/rep{{rep}}/lassoed_cpaths/lasso_path_summary_pos_{{pos}}_"+f"ch_{ch}"+"_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            codebook = f"results/full_codebooks/codebook_ch_{ch}.csv",
        output:
            f"results/plots/best_curve_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.png",
            f"results/plots/best_curve_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.csv"
        params:
            lfs = config["lfs"][ch],
            search_radii = config["search_radii"][ch],
            beta_thresholds = config["beta_thresholds"][ch],
            mpcds=config["mpcd"][ch]
        resources: mem_mb = 35000
        script: "../scripts/plot_cell_ch_best_curve.py"

    rule: # get_best_curves_chx_mp:
        input:
            nnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_nnc_pos_{{pos}}_"+f"ch_{ch}"+"_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{{mpcd}}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
            nnc_betas = expand("results/rep{{rep}}/lassoed_cpaths/betas_nnc_pos_{{pos}}_"+f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{{mpcd}}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
            wnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{{mpcd}}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
            wnc_betas = expand("results/rep{{rep}}/lassoed_cpaths/betas_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{{mpcd}}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
            codebook = f"results/full_codebooks/codebook_ch_{ch}.csv",
        output:
            f"results/plots/best_curve_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}_mpcd{{mpcd}}.csv"
        params:
            lfs = config["lfs"][ch],
            search_radii = config["search_radii"][ch],
            beta_thresholds = config["beta_thresholds"][ch]
        script:
            "../scripts/get_cell_ch_mp_best_curve.py"

    rule: # plot best curve chx mp
        input:
            expand("results/plots/best_curve_rep{{rep}}_pos_{{pos}}_" +f"ch_{ch}" + "_cell_{{cell}}_mpcd{mpcd}.csv", mpcd=config["mpcd"][ch]), 
        output:
            f"results/plots/best_curve_plot_ch_mp_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.png",
        script:
            "../scripts/plot_cell_ch_mp_best_curve.py"

    rule: # get cell best curve results
        input:
            nnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_nnc_pos_{{pos}}_"+f"ch_{ch}"+"_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{{mpcd}}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
            nnc_betas = expand("results/rep{{rep}}/lassoed_cpaths/betas_nnc_pos_{{pos}}_"+f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{{mpcd}}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
            wnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{{mpcd}}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
            wnc_betas = expand("results/rep{{rep}}/lassoed_cpaths/betas_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{{mpcd}}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch]),
            best_ge_eFDR_curve = f"results/plots/best_curve_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.csv",
            paths_summaries = expand("results/rep{{rep}}/lassoed_cpaths/lasso_path_summary_pos_{{pos}}_"+f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{{mpcd}}.csv", lf=config["lfs"], sr=config["search_radii"]),
            codebook = f"results/full_codebooks/codebook_ch_{ch}.csv"
        output:
            f"results/plots/best_curve_plot_ch_mp_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.csv"
        script:
            "../scripts/get_cell_ch_mp_best_curve_all_mp.py" 

    rule: # plot_lf_sr_comparison:
        input:
            nnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_nnc_pos_{{pos}}_"+f"ch_{ch}"+"_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            nnc_beta_paths = expand("results/rep{{rep}}/lassoed_cpaths/betas_nnc_pos_{{pos}}_"+f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            wnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            wnc_beta_paths = expand("results/rep{{rep}}/lassoed_cpaths/betas_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
        output:
            f"results/plots/est_fdr_vs_ge_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.png"
        params:
            lfs = config["lfs"][ch],
            search_radii = config["search_radii"][ch],
            beta_thresholds = config["beta_thresholds"][ch]
        script:
            "../scripts/plot_cell_lf_sr_bt_curves.py"

    rule: # plot_lf_sr_comparison:
        input:
            nnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_nnc_pos_{{pos}}_"+f"ch_{ch}"+"_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            nnc_beta_paths = expand("results/rep{{rep}}/lassoed_cpaths/betas_nnc_pos_{{pos}}_"+f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            wnc_files = expand("results/rep{{rep}}/lassoed_cpaths/lassoed_cpaths_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            wnc_beta_paths = expand("results/rep{{rep}}/lassoed_cpaths/betas_wnc_pos_{{pos}}_" + f"ch_{ch}" + "_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.png", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            paths_summaries = expand("results/rep{{rep}}/lassoed_cpaths/lasso_path_summary_pos_{{pos}}_"+f"ch_{ch}"+"_cell_{{cell}}_lf{lf}_sr{sr}_mpcd{mpcd}.csv", lf=config["lfs"][ch], sr=config["search_radii"][ch], mpcd=config["mpcd"][ch]),
            codebook = f"results/full_codebooks/codebook_ch_{ch}.csv",
            smFISH_ref = "resources/validation/smFISH_results.csv"
        output:
            f"results/plots/est_fdr_vs_sens_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.png",
            f"results/plots/est_fdr_vs_sens_best_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.png",
            f"results/plots/est_fdr_vs_sens_best_sum_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.csv",
        params:
            lfs = config["lfs"][ch],
            search_radii = config["search_radii"][ch],
            beta_thresholds = config["beta_thresholds"][ch]
        resources: mem_mb = 35000
        script: "../scripts/plot_cell_sens_lf_sr_bt_curves.py"

    """
    rule: # plot single channel sensitivity-estimated FDR curves
        input:
            best_curve_pnts = f"results/plots/best_curve_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.csv",#, ch= config["channels"])
            codebook = f"results/full_codebooks/codebook_ch_{ch}.csv",
            smFISH_reference = "resources/validation/smFISH_results.csv"
        output:
            f"results/plots/est_fdr_vs_sens_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.png"
        script:
            "../scripts/plot_single_channel_sensitivity_fdr_curves.py"
    """

    rule: # plot channel sensitivity-estimated FDR curves
        input:
            expand("results/plots/est_fdr_vs_ge_rep{{rep}}_pos_{{pos}}_ch_{ch}_cell_{{cell}}.png", ch= config["channels"])
        output:
            f"results/plots/est_fdr_vs_ge_rep{{rep}}_pos_{{pos}}_all_ch_cell_{{cell}}.png"
        script:
            "../scripts/plot_cell_all_ch_sensitivity_fdr_curves.py"

for rep in config['positions']:
    #aggregate_fov_counts
    rule:
        input: expand("results/"+rep+"/decoded/counts_{{filtered}}_ch_{ch}_pos_{pos}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv", pos=config['positions'][rep], ch=config['channels'])
        output: "results/"+rep+"/cell_barcode_counts/cell_barcode_counts_{filtered}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
        resources: mem_mb=1000
        conda: "../envs/pdnp.yaml"
        script: "../scripts/aggregate_sum_stats.py"


    #rule aggregate_sum_stats:
    rule:
        input: expand("results/"+rep+"/pos_ch_summaries/decode_summary_stats_{{filtered}}_ch_{ch}_pos_{pos}.csv",ch=config['channels'],pos=config['positions'][rep],sf=config['sfs'])
        output: "results/"+rep+"/decode_summary_stats_{filtered}.csv"
        resources: mem_mb = 1000
        conda: "../envs/pdnp.yaml"
        script: "../scripts/aggregate_sum_stats.py"
