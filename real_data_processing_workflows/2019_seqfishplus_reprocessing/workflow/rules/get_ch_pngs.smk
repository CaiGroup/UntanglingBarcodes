
rule get_ch_pngs_643:
    input: "resources/replicates/rep{rep}/HybCycle_{hyb}/MMStack_Pos{pos}.ome.tif"
    output: "results/rep{rep}/ch_pngs/HybCycle_{hyb}_ch_643_pos_{pos}.png"
    params: ch=643, z=config['z_slice']
    resources: mem_mb = 5000
    group: "fit_beads"
    conda: '../envs/get_ch_pngs.yaml'
    script: "../scripts/get_ch_pngs.py"

rule get_ch_pngs_561:
    input: "resources/replicates/rep{rep}/HybCycle_{hyb}/MMStack_Pos{pos}.ome.tif"
    output: "results/rep{rep}/ch_pngs/HybCycle_{hyb}_ch_561_pos_{pos}.png"
    params: ch=561, z=config['z_slice']
    resources: mem_mb = 5000
    group: "fit_beads"
    conda: '../envs/get_ch_pngs.yaml'
    script: "../scripts/get_ch_pngs.py"

#separate rule to fix switched hybs
rule get_ch_pngs_correct_488_rep1:
    input: expand("resources/replicates/rep1/HybCycle_{hyb}/MMStack_Pos{{pos}}.ome.tif", hyb = hybs)
    output: expand("results/rep1/ch_pngs/HybCycle_{hyb}_ch_488_pos_{{pos}}.png", hyb = hybs)
    params: ch=488, z=config['z_slice']
    resources: mem_mb = 5000
    conda: '../envs/get_ch_pngs.yaml'
    script: "../scripts/get_ch_pngs_488.py"

rule get_ch_pngs_correct_488_rep2:
    input: expand("resources/replicates/rep2/HybCycle_{hyb}/MMStack_Pos{{pos}}.ome.tif", hyb = hybs)
    output: expand("results/rep2/ch_pngs/HybCycle_{hyb}_ch_488_pos_{{pos}}.png", hyb = hybs)
    params: ch=488, z=config['z_slice']
    resources: mem_mb = 5000
    conda: '../envs/get_ch_pngs.yaml'
    script: "../scripts/get_ch_pngs_488.py"

rule get_bead_pngs:
    input: "resources/replicates/rep{rep}/beads_{bead}/MMStack_Pos{pos}.ome.tif"
    output: "results/rep{rep}/ch_pngs/beads_{bead}_ch_{ch}_pos_{pos}.png"
    params: z=config['z_slice']
    resources: mem_mb = 5000
    group: "fit_beads"
    conda: '../envs/get_ch_pngs.yaml'
    script:"../scripts/get_ch_pngs.py"
    #notebook: "../notebooks/get_ch_pngs.py.ipynb" #"scripts/to_dlm.py"
