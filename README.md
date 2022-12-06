# The impact of violating the independence assumption in meta-analysis on biomarker discovery

In this study, we (1) review and compare different meta-analyses to estimate variations across studies along with biomarker discoveries using preclinical pharmacogenomics data, and (2) evaluate the performance of conventional meta-analysis where the dependence of the effects was ignored via simulation studies and pharmacogenomics data (breast and pan-cancer). 

## Data

We used transcriptomic (RNA-Sequencing and gene expression microarray) and drug response data from pharmacogenomic cancer cell line sensitivity screenings, including 

- Cancer Cell Line Encyclopedia (CCLE: Broad-Novartis)
- Genomics of Drug Sensitivity in Cancer (GDSC: Wellcome Trust Sanger Institute)
- Genentech Cell Line Screening Initiative (gCSI)
- Cancer Therapeutics Response Portal (CTRP: Broad Institute)
- Oregon Health and Science University breast cancer screen (GRAY)
- University Health Network Breast Cancer Screen (UHNBreast) 

Molecular information was obtained from the PharmacoGx R package along with details on data processing. Cell line drug response data, in the form of area above the curve (AAC) recomputed information, was also obtained from the PharmacoGx R package.  
