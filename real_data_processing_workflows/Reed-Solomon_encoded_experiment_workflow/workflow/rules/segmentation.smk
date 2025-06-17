"""
rule run_cellpose:
    input:
        "resources/RS_seqfish_half_pool/HybCycle_52/HybCycle_52_MMStack_Pos{pos}.ome.tif"
    output:
        "results/labeled_stacks/labeled_stack_pos{pos}.tif"
    conda: "../envs/segmentation.yaml"
    script: "../scripts/run_segmentation.py"
"""

rule get_segmentation_training_imgs:
    input: "resources/RS_seqfish_half_pool/HybCycle_35/MMStack_Pos{pos}.ome.tif"
    output:
        "results/seg_training/pos{pos}_polyT_dapi_z0.tiff",
        "results/seg_training/pos{pos}_polyT_dapi_z8.tiff",
        "results/seg_training/pos{pos}_polyT_dapi_z16.tiff"
    script: "../scripts/get_segmentation_training_imgs.py"


rule get_segmentation_stacks:
    input: "resources/RS_seqfish_half_pool/HybCycle_52/HybCycle_52_MMStack_Pos{pos}.ome.tif"
    output:
        "results/seg_stacks/pos{pos}_polyT_dapi.tiff"
    script: "../scripts/get_segmentation_stacks.py"

rule filter_labeled_regions:
    input: "resources/labeled_images/pos{pos}_polyT_dapi_cp_masks.png"
    output: "results/filtered_labeled_imgs/lbld_img_pos{pos}.png"
    script: "../scripts/filter_labeled_images.py"