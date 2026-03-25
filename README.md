# NatCom_2026_WM_Multiomics
This repository contains the R scripts used to generate the figures, supplementary figures, and supplementary tables found in **Hunter et al The Evolution and Subtypes of Waldenstrom Macroglobulinemia: Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients. 2026. Nat Comm.**

Code is organized by figure panel. Supplementary figures and tables are embedded with the corresponding figure panel from the main text. All elements included in each script are listed in a commented header section for convenience. These scripts make extensive use of the **R** package ‘bvt’ (*https://github.com/ZachHunter/bvt*). This is a set of tools used to streamline data visualization, particularly from Bioconductor data structures. All code is publicly available and can be installed by running `devtools::install_github(“ZachHunter/bvt”)`. 

## Data Availability

- Data from this manuscript is available through the European Genome-Phenome Archive and accessible at *Link Pending*
- Healthy donor immune cell transcripts per million data was obtained from the Human Protein Atlas and is available for download at  (*https://www.proteinatlas.org/humanproteome/single+cell/immune+cell/data#immune_cells_hpa* ).<SUP>1</SUP>
- Gene expression profiling data used in EScore and subtype gene signature validation is described in Trojani et al.<SUP>2</SUP> (GSE171739; *https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171739*) and accessed using the geoQuery package in **R**.
- The single cell RNASeq data described in Hao et al<SUP>3</SUP> is available through the European Genome-Phenome Archive and accessible at *https://ega-archive.org/dacs/EGAC50000000863*.

## References

1.	Uhlen, M. et al. A genome-wide transcriptomic analysis of protein-coding genes in human blood cells. Science (1979). 366, (2019).
2.	Trojani, A. et al. Identification of a candidate gene set signature for the risk of progression in igm mgus to smoldering/symptomatic waldenström macroglobulinemia (Wm) by a comparative transcriptome analysis of b cells and plasma cells. Cancers (Basel). 13, (2021).
3.	Sun, H. et al. Evolution of tumor subclones and T-cell dynamics underlie variable ibrutinib responses in Waldenström macroglobulinemia. Blood https://doi.org/10.1182/BLOOD.2025032268 (2026) doi:10.1182/BLOOD.2025032268.
 
