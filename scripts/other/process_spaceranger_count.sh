#!/bin/bash

### Sample: NB1

# Generate  sample root names
parallel -j 1 echo {1}_{2} ::: NB_LU_01_Pre NB_LU_02_Post ::: 1 2  > raw/root_ids.csv    

# run spaceranger for multiple samples in parallel
cat raw/root_ids.csv | parallel \
spaceranger count --id={}\
                   --transcriptome=downloads/genomes/spaceranger_refdata/refdata-gex-GRCh38-2020-A\
                   --probe-set=genomes/spaceranger_refdata/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv\
                   --fastqs=datasets/{}_fastqs\
                   --image=datasets/{}_image.tif\
                   --unknown-slide=visium-1\
                   --localcores=16\
                   --localmem=30  
                 
### Sample: NB2
# Generate  sample root names
parallel -j 1 echo {1}{2} ::: NBLU02Pre NBLU02Post ::: 1 2  > raw/root_ids2.csv                 

# Run spaceranger count
cat raw/root_ids2.csv | parallel spaceranger count --id={} --transcriptome= downloads/genomes/spaceranger_refdata/refdata-gex-GRCh38-2020-A --probe-set=downloads/genomes/spaceranger_refdata/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv --fastqs=raw/P30603/NBLU02/{}_fastqs --cytaimage=raw/P30603/NBLU02/{}_image.tif --unknown-slide=visium-2 --localcores=16 --localmem=30  



