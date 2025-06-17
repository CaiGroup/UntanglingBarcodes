#min_params = "lf{lf}_zf_{zf}_wf{wf}_sf{sf}_dr{dr}".format(lf=min(config['lfs']),zf=min(config['zfs']),wf=min(config['wfs']),sf=min(config['sfs']),dr=max(config['drops']))
min_params = "lf{lf}_zf_{zf}_wf{wf}_sf{sf}_dr{dr}".format(lf=3.0,zf=min(config['zfs']),wf=min(config['wfs']),sf=min(config['sfs']),dr=max(config['drops']))

"""
rule make_interactive_overlay_hyb_max_proj:
    input:
        "results/decoded/points_ch_{ch}_pos_{pos}_" + min_params + ".csv",
        "results/alignment/offsets_pos_{pos}.csv",
        expand("results/ims_bg_sub/HybCycle_{hyb}_ch{{ch}}_pos_{{pos}}.tif""results/ch_stacks/HybCycle_{hyb}_ch_{{ch}}_pos{{pos}}.tif", hyb=hybs)
    params: roi_width = 50
    output: "results/plots/hyb_max_proj_overlays/hyb_overlay_ch_{ch}_pos_{pos}.html"
    script: "../scripts/make_interactive_overlay_max_proj_auto_roi.py"
"""

rule plot_fit_overlay_bgsub_auto_roi:
    input:
        #"results/fit_dots/HybCycle_{hyb}_ch_{ch}_pos{pos}_mw_{mw}.csv",
        "results/fit_dots/HybCycle_{hyb}_ch_{ch}_pos{pos}.csv",
        "results/ims_bg_sub/HybCycle_{hyb}_ch{ch}_pos_{pos}.tif"
    params: roi_width=100
    output: "results/plots/fit_overlays/bgsub_overlay_{hyb}_ch_{ch}_pos_{pos}.html"
    wildcard_constraints: ch="\d+", pos="\d+"
    resources: mem_mb = 10000
    conda: "../envs/make_interactive_overlay.yaml"
    script: "../scripts/make_interactive_overlay_auto_roi.py"


rule make_on_off_target_plots:
    input: "results/decode_summary_stats_pos_{pos}.csv"
    output: "results/plots/on_off_plots/on_off_pos_{pos}.png"
    wildcard_constraints: ch="\d+", pos="\d+"
    resources: mem_mb = 1000
    conda: "../envs/plot_on_off_params.yaml"
    script: "../scripts/plot_on_off_params.py"

rule plot_cross_channel_registration_loocv:
    input: "results/alignment/ref_cross_channel_loocv_pos{pos}.csv"
    output: "results/plots/registrtion_loocv/on_off_pos_{pos}.png"
    resources: mem_mb = 100
    conda: "../envs/plot_on_off_params.yaml"
    script: "../scripts/plot_cross_channel_registration_loocv.py"
