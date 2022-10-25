# The impact of violating the independence assumption in meta-analysis on biomarker discovery

With rapid advancements in high-throughput sequencing technologies, massive amounts of ‘-omics’ data are now available in almost every biomedical field. Due to variance in biological models and analytic methods, findings from clinical and biological studies are often not generalizable when tested in independent cohorts. Meta-analysis, a set of statistical tools to integrate independent studies addressing similar research questions, has been proposed to improve the accuracy and robustness of new biological insights. However, it is common practice among biomarker discovery studies using preclinical pharmacogenomic data to borrow molecular profiles of cancer cell lines from one study to another, creating dependence across studies. The impact of violating the independence assumption in meta-analyses is largely unknown. In this repo, we (1) review and compare different meta-analyses to estimate variations across studies along with biomarker discoveries using preclinical pharmacogenomics data, and (2) evaluate the performance of conventional meta-analysis where the dependence of the effects was ignored via simulation studies and pharmacogenomics data (breast and pan-cancer). 

# Data

We used transcriptomic (RNA-Sequencing and gene expression microarray) and drug response data from pharmacogenomic cancer cell line sensitivity screenings, including the Cancer Cell Line Encyclopedia (CCLE: Broad-Novartis), the Genomics of Drug Sensitivity in Cancer (GDSC: Wellcome Trust Sanger Institute), the Genentech Cell Line Screening Initiative (gCSI), the Cancer Therapeutics Response Portal (CTRP: Broad Institute), Oregon Health and Science University breast cancer screen (GRAY), and University Health Network Breast Cancer Screen (UHNBreast). Molecular information was obtained from the \textbf{\textit{PharmacoGx} R package} along with details on data processing. Cell line drug response data, in the form of area above the curve (AAC) recomputed information, was also obtained from the \textbf{\textit{PharmacoGx} R package} 

