# NatCom_2026_WM_Multiomics
This repository contains the R scripts used to generate the figures, supplementary figures, and supplementary tables found in **Hunter et al The Evolution and Subtypes of Waldenstrom Macroglobulinemia: Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients. 2026. Nat Comm.**

Code is organized by figure panel. Supplementary figures and tables are embedded with the corresponding figure panel from the main text. All elements included in each script are listed in a commented header section for convenience. These scripts make extensive use of the **R** package ‘bvt’ (*https://github.com/ZachHunter/bvt*). This is a set of tools used to streamline data visualization, particularly from Bioconductor data structures. All code is publicly available and can be installed by running `devtools::install_github(“ZachHunter/bvt”)`.

## Guide to the files

Each file generates the associated main figure panel in addition to related supplementary content list as subheadings below.
- **NC_2026_Figure1.R:** *Identification of WM Subtypes*
  - Supplementary Table 1: *List of heatmap genes used to indentify WM subtypes*
  - Supplementary Table 2: *Clinical associations with Early WM and subtyped WM (BCL/PCL)*
  - Supplementary Table 3: *Clinical associations with BCL and PCL subtyped WM*
- **NC_2026_Figure2.R:** *Identification of the WM Evolutionary Score*
  - Supplementary Figure 1: *Exploring EScore and WM Subtype*
  - Supplementary Table 4: *Signature genes for EScore and WM Subtype from NMF analysis*
  - Supplementary Table 5: *Clinical associations with Early and Late EScore*
- **NC_2026_Figure3.R:** *Clinical associations and external validation of Subtype and EScore*
- **NC_2026_Figure4.R:** *Subtype signature gene expression and validation*
  - Supplementary Figure 2: *WM subtype gene expression*
  - BCL versus PCL differential gene expression analysis
- **NC_2026_Figure5.R:** *Characterization of the BCL and PCL WM subtypes*
  - Supplementary Figure 3: *Confirmation of WM Subtype findings*
  - Supplementary Table 6: *Clinical associations with EBCL and BCL subtyped WM*
  - Supplementary Table 7: *Clinical associations with EPCL and PCL subtyped WM*
  - Supplementary Table 8: *First therapy received by WM subtype*
  - Supplementary Table 9: *Type of first therapy received by WM subtype*
  - EBCL versus BCL differential gene expression analysis
  - EPCL versus PCL differential gene expression analysis
- **NC_2026_Figure6.R:** *Characterization of the WM EScore*
  - Supplementary Figure 4: *Analysis of the WM EScore*
  - Early versus Late EScore differential gene expression analysis
  - EScore GAM differential gene expression analysis
- **NC_2026_Figure7.R:** *Non-B lineage gene expression in early EScore WM*
- **NC_2026_Figure8.R:** *Validation of non-B lineage expression in early EScore WM*
- **NC_2026_Figure9.R:** *EScore validation studies*
  - Supplementary Figure 5: *Additional supportive data for EScore validation*
- **NC_2026_Figure10.R:** *Clinical significance of the WM EScore*

## Data Availability

- Data from this manuscript is available through the European Genome-Phenome Archive and accessible at *Link Pending*
- Healthy donor immune cell transcripts per million data was obtained from the Human Protein Atlas and is available for download at  (*https://www.proteinatlas.org/humanproteome/single+cell/immune+cell/data#immune_cells_hpa* ).<SUP>1</SUP>
- Gene expression profiling data used in EScore and subtype gene signature validation is described in Trojani et al.<SUP>2</SUP> (GSE171739; *https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171739*) and accessed using the geoQuery package in **R**.
- The single cell RNASeq data described in Hao et al<SUP>3</SUP> is available through the European Genome-Phenome Archive and accessible at *https://ega-archive.org/dacs/EGAC50000000863*.

## References

1.	Uhlen, M. et al. A genome-wide transcriptomic analysis of protein-coding genes in human blood cells. Science (1979). 366, (2019).
2.	Trojani, A. et al. Identification of a candidate gene set signature for the risk of progression in igm mgus to smoldering/symptomatic waldenström macroglobulinemia (Wm) by a comparative transcriptome analysis of b cells and plasma cells. Cancers (Basel). 13, (2021).
3.	Sun, H. et al. Evolution of tumor subclones and T-cell dynamics underlie variable ibrutinib responses in Waldenström macroglobulinemia. Blood https://doi.org/10.1182/BLOOD.2025032268 (2026) doi:10.1182/BLOOD.2025032268.
 
