
hybs_to_stack = list(range(config['initial_hyb'],config['final_hyb']))

rule stack_aligned_images:
    input:
        "results/rep{rep}/alignment/offsets_ch_{ch}_pos_{pos}.csv",
        'results/rep{rep}/ch_pngs/beads_initial_ch_{ch}_pos_{pos}.png',
        expand('results/rep{{rep}}/ch_pngs/HybCycle_{hyb}_ch_{{ch}}_pos_{{pos}}.png', hyb=hybs_to_stack)
        #expand('ch_pngs/HybCycle_{hyb}_ch_{{ch}}_pos_{{pos}}.png', hyb=list(range(1,81)))
    output: "results/rep{rep}/plots/aligned_stacks/aligned_stack_ch_{ch}_pos_{pos}.tif"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb = 20000
    conda: "../envs/get_ch_pngs.yaml"
    script: "../scripts/stack_aligned_images.py"

rule stack_aligned_bg_sub_images_dapi:
    input:
        "results/rep{rep}/alignment/offsets_ch_{ch}_pos_{pos}.csv",
        'results/rep{rep}/ch_pngs/beads_initial_ch_{ch}_pos_{pos}.png',
        expand("results/rep{{rep}}/ims_bg_sub/HybCycle_{hyb}_ch_{{ch}}_pos_{{pos}}.png",hyb=hybs)
    params: nhybs = config['final_hyb']-config['initial_hyb']+1,
            first_hyb =config['initial_hyb'],
            zslice = config['z_slice']
    output: "results/rep{rep}/plots/aligned_stacks/aligned_stack_bgsub_ch_{ch}_pos_{pos}.tif"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb = 20000
    conda: "../envs/get_ch_pngs.yaml"
    script: "../scripts/stack_aligned_images.py"

rule clean_up_2019_results:
    input:
        "resources/seqFISH+_NIH3T3_point_locations/all_gene_Names.mat",
        "resources/seqFISH+_NIH3T3_point_locations/RNA_locations_run_1.mat",
        "resources/seqFISH+_NIH3T3_point_locations/RNA_locations_run_2.mat"
    output:
        "results/seqFISH+_NIH3T3_point_locations_2019/RNA_locations_run_1.csv",
        "results/seqFISH+_NIH3T3_point_locations_2019/RNA_locations_run_2.csv"
    resources: mem_mb=10000
    script: "../scripts/read_published_results_reformat.py"

rule make_correlation_plots:
    input:
        "results/rep{rep}/decoded/barcodes{filtered}ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/seqFISH+_NIH3T3_point_locations_2019/RNA_locations_run_{rep}.csv",
        "resources/validation_files/smFISH_results.csv",
        "results/full_codebooks/codebook_ch_{ch}.csv",
        "results/rep{rep}/Labeled_Images/MMStack_Pos{pos}.png",
        "results/rep{rep}/alignment/aligned_dots_ch_{ch}_pos_{pos}.csv"
    resources: mem_mb = 2000
    output:
        "results/rep{rep}/plots/corr_plots/corr_plot{filtered}ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png",
        "results/rep{rep}/plots/corr_plots/pubcomp{filtered}ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/rep{rep}/plots/corr_plots/smFISHcomp{filtered}ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"#,
        #"results/rep{rep}/plots/corr_plots/corr_plot_filtered_ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    group: 'plot_pos_ch'
    conda: "../envs/make_corr_plots.yaml"
    resources: mem_mb=100
    script: "../scripts/make_corr_plots.py"

rule make_interactive_overlay:
    input:
        #"results/rep{rep}/decoded/points_ch_{ch}_pos_{pos}_"+min_params+".csv",
        "results/rep{rep}/decoded/points_no_neg_ctrl_ch_{ch}_pos_{pos}_"+plot_params+".csv",
        "results/rep{rep}/plots/aligned_stacks/aligned_stack_ch_{ch}_pos_{pos}.tif",
        "results/seqFISH+_NIH3T3_point_locations_2019/RNA_locations_run_{rep}.csv",
        "resources/codebooks/channel_{ch}.csv",#"validation_files/10k_codebook_ch_{ch}.csv",
        "results/full_codebooks/codebook_ch_{ch}.csv"#"codebooks/E2019_cb_all_control_ch_{ch}.txt",
    params: roi_width=config['overlay_roi_width']
    output: "results/rep{rep}/plots/interactive_plots/overlay_{filtered}_ch_{ch}_pos_{pos}_"+plot_params+".html" #"results/rep{rep}/plots/interactive_plots/overlay_{filtered}_ch_{ch}_pos_{pos}_"+min_params+".html"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb = 10000
    conda: "../envs/make_interactive_overlay.yaml"
    script: "../scripts/make_interactive_overlay_auto_roi.py"

rule make_interactive_overlay_bgsub_auto_roi:
    input:
        #"results/rep{rep}/decoded/points_ch_{ch}_pos_{pos}_"+min_params+".csv",
        "results/rep{rep}/decoded/points_no_neg_ctrl_ch_{ch}_pos_{pos}_"+plot_params+".csv",
        "results/rep{rep}/plots/aligned_stacks/aligned_stack_bgsub_ch_{ch}_pos_{pos}.tif",
        "results/seqFISH+_NIH3T3_point_locations_2019/RNA_locations_run_{rep}.csv",
        "resources/codebooks/channel_{ch}.csv",#"validation_files/10k_codebook_ch_{ch}.csv",
        "results/full_codebooks/codebook_ch_{ch}.csv"#"codebooks/E2019_cb_all_control_ch_{ch}.txt",
    params: roi_width=200
    output: "results/rep{rep}/plots/interactive_plots/bgsub_overlay_{filtered}_ch_{ch}_pos_{pos}_"+plot_params+".html" #"results/rep{rep}/plots/interactive_plots/bgsub_overlay_{filtered}_ch_{ch}_pos_{pos}_"+min_params+".html"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb = 10000
    conda: "../envs/make_interactive_overlay.yaml"
    script: "../scripts/make_interactive_overlay_auto_roi.py"

rule make_on_off_target_plots:
    input: "results/rep{rep}/pos_ch_summaries/decode_summary_stats_{filtered}_ch_{ch}_pos_{pos}.csv"
    output: "results/rep{rep}/plots/on_off_plots/on_off_{filtered}_ch_{ch}_pos_{pos}.png"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb = 10
    conda: "../envs/plot_on_off_params.yaml"
    script: "../scripts/plot_on_off_params.py"

rule make_overall_on_off_target_plot:
    input: "results/rep{rep}/decode_summary_stats_{filtered}.csv"
    output: "results/rep{rep}/plots/on_off_plots/on_off_{filtered}_overall.png"
    wildcard_constraints: rep="\d+"
    resources: mem_mb = 10
    conda: "../envs/plot_on_off_params.yaml"
    script: "../scripts/plot_on_off_params_overall.py"

rule plot_rmse_loocv_alignment_errors:
    input: expand("results/rep{{rep}}/alignment/loov_errors_ch_{ch}_pos_{{pos}}.csv", ch=config['channels'])
    output: "results/rep{rep}/plots/alignment_loocv_errors/pos_{pos}.png"
    wildcard_constraints: rep="\d+", pos="\d+"
    conda: "../envs/make_corr_plots.yaml"
    resources: mem_mb=10
    script: "../scripts/plot_hyb_rmse_loocv_errors.py"

rule plot_radial_density_function:
    input:
        expand('results/full_codebooks/codebook_ch_{ch}.csv', ch=config['channels']),
        expand_rep_all("results/rep{rep}/decoded/points_no_neg_ctrl_ch_{{ch}}_pos_{{pos}}_"+plot_params+".csv", config["reps"])
        #expand_rep_all("results/rep{rep}/decoded/points_ch_{{ch}}_pos_{{pos}}_"+min_params+".csv", config["reps"])
        #expand("results/rep{{rep}}/decoded/points_ch_{ch}_pos_{pos}_"+min_params+".csv", ch=config['channels'],pos=config['positions'])
    params:
        nchannels=len(config['channels']),
        rmax=config['radial_density_func_rmax'],
        delta_r=config['radial_density_func_delta_r']
    output:
        'results/rep{rep}/plots/dot_stats/radial_density_functions.csv',
        'results/rep{rep}/plots/dot_stats/radial_density_functions.png'
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    conda: "../envs/radial_density.yaml"
    resources: mem_mb=5000
    script: '../scripts/get_radial_density_functions.py'

rule plot_dot_stats:
    input:
        'results/full_codebooks/codebook_ch_{ch}.csv',
        "results/rep{rep}/decoded/points_no_neg_ctrl_ch_{ch}_pos_{pos}_"+plot_params+".csv"
        #"results/rep{rep}/decoded/points_ch_{ch}_pos_{pos}_"+min_params+".csv"
    output:
        'results/rep{rep}/plots/dot_stats/hyb_mean_intensities_pos_{pos}_ch_{ch}.png',
        'results/rep{rep}/plots/dot_stats/hyb_counts_pos_{pos}_ch_{ch}.png',
        #'results/plots/dot_stats/ecdfs/hyb_ecdf_pos_{pos}_ch_{ch}.png'
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    conda: "../envs/plot_on_off_params.yaml"
    group: 'plot_pos_ch'
    resources: mem_mb=500
    script: '../scripts/plot_hyb_dot_stats.py'

rule plot_overall_corrs:
    input:
        "results/seqFISH+_NIH3T3_point_locations_2019/RNA_locations_run_{rep}.csv",
        "results/rep{rep}/cell_barcode_counts/cell_barcode_counts_{filtered}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
    output:
        expand("results/rep{{rep}}/plots/corr_plots/corr_plot_{{filtered}}_ch_{ch}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.png", ch=config['channels']),
        "results/rep{rep}/plots/corr_plots/corr_plot_overall_{filtered}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png",
        #"results/rep{rep}/plots/corr_plots/rnaseq_comp_overall_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png"
    wildcard_constraints: rep="\d+"
    conda: "../envs/make_corr_plots.yaml"
    resources: mem_mb=2000
    script: "../scripts/plot_overall_corrs.py"

rule plot_all_rep_corrs:
    input:
        expand("results/seqFISH+_NIH3T3_point_locations_2019/RNA_locations_run_{rep}.csv", rep=config['reps']),
        expand("results/rep{rep}/cell_barcode_counts/cell_barcode_counts_{{filtered}}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv", rep=config['reps']),
        "resources/validation_files/smFISH_results.csv"
    output:
        expand("results/plots/corr_plots/corr_plot_{{filtered}}_ch_{ch}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.png", ch=config['channels']),
        "results/plots/corr_plots/corr_plot_overall_{filtered}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png",
        "results/plots/corr_plots/corr_plot_overall_smfish_{filtered}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png",
        #"results/rep{rep}/plots/corr_plots/rnaseq_comp_overall_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png"
    wildcard_constraints: rep="\d+"
    conda: "../envs/make_corr_plots.yaml"
    resources: mem_mb=3000
    script: "../scripts/plot_all_rep_overall_corrs.py"

for rep in config['positions']:
    #rule plot_nc_freq_hist:
    rule:
        input:
            expand("results/"+rep+"/decoded/barcodes_unfiltered_ch_{ch}_pos_{pos}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv",ch=config['channels'], pos=config['positions'][rep])
        output:
            "results/"+rep+"/plots/nc_freq_hists/nc_freq_hist_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png",
            "results/"+rep+"/plots/nc_freq_hists/nc_freq_hist_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
            "results/"+rep+"/plots/nc_freq_hists/nc_cell_freq_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
            "results/"+rep+"/plots/nc_freq_hists/nc_freq_hist_filt_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png",
            #"results/rep{rep}/plots/nc_freq_hists/nc_freq_hist_filt_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
            #"results/rep{rep}/plots/nc_freq_hists/nc_cell_freq_filt_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
        conda: "../envs/make_corr_plots.yaml"
        group: 'plot_pos_ch'
        resources: mem_mb=1000
        script: '../scripts/plt_nc_freq_hist.py'

rule test_gene_cell_cnts_poisson:
    input: "results/rep{rep}/plots/nc_freq_hists/nc_cell_freq_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
    output:
        "results/rep{rep}/plots/cell_pdt_hist/nc_codeword-cell_count_hist_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png",
        "results/rep{rep}/plots/cell_pdt_hist/cell_pdt_hist_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png",
        "results/rep{rep}/plots/cell_pdt_hist/cell_bc_mu_var_lf{lf}_wf{wf}_sf{sf}_dr{dr}.png"
    conda: "../envs/make_corr_plots.yaml"
    group: 'plot_pos_ch'
    resources: mem_mb=1000
    script: "../scripts/test_nc_cws_Poisson.py"
