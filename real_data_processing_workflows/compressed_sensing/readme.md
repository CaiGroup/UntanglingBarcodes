# Compressed sensing example workflow

This workflow uses the preprocessed images from the spots-first 2019_seqfishplus_reprocessing workflow as its inputs. To run this workflow, get the preprocessed images from the FigShare version of the 2019_seqfishplus_reprocessing workflow, and add them to results folder of this workflow so that the preprocessed images have paths relative to the main folder of

```
results/{rep}/ims_bg_sub/HybCycle_{hyb}_ch_{ch}_pos_{pos}.png
```

and add the labeled images to the resources folder so that they have paths of 

```
resources/replicates/{rep}/Labeled_Images/MMStack_Pos{pos}.png
```


Once the preprocessed and labeled images are in place, the compressed sensing snakemake workflow can be run using the same commands as the spots-first workflows.