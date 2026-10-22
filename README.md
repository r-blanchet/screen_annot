# Screen annotation tool

Annotate a vcf file from SCREEN Encode database.

Usage on Glicid:
```
module load apptainer
SCREEN_DATADIR="/LAB-DATA/GLiCID/projects/ITX/DOUBLE_DISCORDANCE/SCREEN/data"
apptainer run raphb44/screen_annot:latest .screen_annot \
  --invcf <Your vcf> \
  --ccres $SCREEN_DATADIR/GRCh38-cCREs.bed.gz \
  --gene_link $SCREEN_DATADIR/V4-hg38.Gene-Links.3D-Chromatin.bed.gz,$SCREEN_DATADIR/V4-hg38.Gene-Links.CRISPR.bed.gz,$SCREEN_DATADIR/V4-hg38.Gene-Links.eQTLs.bed.gz \
  --biodata CTCF:$SCREEN_DATADIR/GRCh38.CTCF-zscore.rDHS-V3.bed.gz,DNASE:$SCREEN_DATADIR/GRCh38.DNase-zscore.rDHS-V3.bed.gz,H3K27ac:$SCREEN_DATADIR/GRCh38.H3K27ac-zscore.rDHS-V3.bed.gz,H3k4me3:$SCREEN_DATADIR/GRCh38.H3K4me3-zscore.rDHS-V3.bed.gz
```

This will add in the following fields in the VCF:
- INFO/ScreenElements: The cCREs from screen overlapping the variant and their nature (distal enhancer, promoter, ....)
- INFO/ScreenGenes: The genes associated with the overlapping cCREs and the method used to identify this association
- INFO/ScreenBiodata: The biosamples for which the zscore is >2 (top 95%) for each of the provided biodata source.  
