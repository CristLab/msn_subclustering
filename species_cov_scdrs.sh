#!/bin/bash

if [ -f /etc/profile.d/modules.sh ]; then
	source /etc/profile.d/modules.sh
fi

module load python/3.9.1
export PATH="/home/crist/.local/bin/:$PATH"

scdrs compute-score \
    --h5ad-file /home/crist/nac_rn7/scDRS/human_msn.h5ad\
    --h5ad-species human\
    --gs-file /home/crist/scDRS/sud_zscores.gs\
    --gs-species human\
    --out-folder /home/crist/scDRS/msn_manuscript/human/\
	--cov-file /home/crist/nac_rn7/scDRS/human_msn.cov\
    --flag-filter-data True\
    --flag-raw-count True\
    --n-ctrl 1000\
    --flag-return-ctrl-raw-score False\
    --flag-return-ctrl-norm-score True

scdrs perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/human_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/human/AUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/human/\
    --group-analysis cellType \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True
	
scdrs perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/human_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/human/AUDITC.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/human/\
    --group-analysis cellType \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True

scdrs perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/human_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/human/OUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/human/\
    --group-analysis cellType \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True

scdrs perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/human_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/human/TUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/human/\
    --group-analysis cellType \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True
	
scdrs compute-score \
    --h5ad-file /home/crist/nac_rn7/scDRS/mouse_msn.h5ad\
    --h5ad-species mouse\
    --gs-file /home/crist/scDRS/sud_zscores.gs\
    --gs-species human\
    --out-folder /home/crist/scDRS/msn_manuscript/mouse/\
	--cov-file /home/crist/nac_rn7/scDRS/mouse_msn.cov\
    --flag-filter-data True\
    --flag-raw-count True\
    --n-ctrl 1000\
    --flag-return-ctrl-raw-score False\
    --flag-return-ctrl-norm-score True

scdrs perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/mouse_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/mouse/AUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/mouse/\
    --group-analysis orig.ident \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True
	
scdrs perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/mouse_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/mouse/AUDITC.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/mouse/\
    --group-analysis orig.ident \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True

scdrs perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/mouse_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/mouse/OUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/mouse/\
    --group-analysis orig.ident \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True

scdrs perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/mouse_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/mouse/TUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/mouse/\
    --group-analysis orig.ident \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True


scdrs_rat compute-score \
    --h5ad-file /home/crist/nac_rn7/scDRS/rat_msn.h5ad\
    --h5ad-species rat\
    --gs-file /home/crist/scDRS/sud_zscores.gs\
    --gs-species human\
    --out-folder /home/crist/scDRS/msn_manuscript/rat/\
	--cov-file /home/crist/nac_rn7/scDRS/rat_msn.cov\
    --flag-filter-data True\
    --flag-raw-count True\
    --n-ctrl 1000\
    --flag-return-ctrl-raw-score False\
    --flag-return-ctrl-norm-score True

scdrs_rat perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/rat_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/rat/AUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/rat/\
    --group-analysis seurat_clusters \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True
	
scdrs_rat perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/rat_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/rat/AUDITC.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/rat/\
    --group-analysis seurat_clusters \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True

scdrs_rat perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/rat_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/rat/OUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/rat/\
    --group-analysis seurat_clusters \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True

scdrs_rat perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/rat_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/rat/TUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/rat/\
    --group-analysis seurat_clusters \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True


scdrs_rat perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/rat_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/rat/AUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/rat/\
    --group-analysis subclusters \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True
	
scdrs_rat perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/rat_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/rat/AUDITC.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/rat/\
    --group-analysis subclusters \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True

scdrs_rat perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/rat_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/rat/OUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/rat/\
    --group-analysis subclusters \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True

scdrs_rat perform-downstream \
    --h5ad-file /home/crist/nac_rn7/scDRS/rat_msn.h5ad\
    --score-file /home/crist/scDRS/msn_manuscript/rat/TUD.full_score.gz\
    --out-folder /home/crist/scDRS/msn_manuscript/rat/\
    --group-analysis subclusters \
	--gene-analysis \
    --flag-filter-data True\
    --flag-raw-count True