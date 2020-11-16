# Re-identification of patient subgroups in Uveal Melanoma
---
This repository contains source code and original/preprocessed datasets for the paper "Re-identification of patient subgroups in Uveal Melanoma". 
#### 1. Introduction
---
Uveal melanoma (UM) is a comparatively rare cancer but requires significant consideration since patients with developing metastatic diagnosis survive only for about 6-12 months. Fortunately, increasingly large multi-omics databases allow us to further understand cancer initiation and development. Moreover, previous studies have observed that the association between copy number aberrations (CNA) and methylation (MET) has affected these processes. From that, we decided to explore the effect of this association on a case study of UM. Also, the current subtypes of UM display its weak association with biological phenotypes and its lack of therapy suggestions. Therefore, the re-identification of molecular subtype is essential for UM patients. In this study, we recruited three omics profiles, including CNA, MET, and mRNA expression, in a UM cohort from The Cancer Genome Atlas (TCGA). Firstly, we identified two sets of genes, CNAexp and METexp, whose CNA and MET significantly correlated with their corresponding expression levels, respectively. Then, single and integrative analyses of the three data types were performed using a tool named PINSPlus. As a result, we discovered two novel integrative subgroups (i.e., IntSub1 and IntSub2), which could be a useful alternative classification for UM patients in the future. To further explore molecular events behind each subgroup, we identified subgroup-specific genes computationally. Accordingly, IntSub1-specific genes were highly associated with cellular cation homeostasis, which responds effectively to chemotherapy, using ion channel inhibitors drugs. On the other side, the most highest expressed genes among IntSub2-specific genes were mostly enriched with immune-related processes. In addition, we detected that different integrative subgroups showed different age-related risks and survival rates. These discoveries can influence the frequency of metastatic surveillance, support medical practitioner to choose appropriate strategies.

**NOTE:** All statistical analyses were performed using R statistical software.

#### 2. Analysis Pipeline
---
![Figure1](https://imgur.com/M0YryNG.png)
**Figure 1. Analysis pipeline.** Firstly, CNA and MET datasets with the corresponding mRNA data were inputted to the function [geneCor](https://github.com/huynguyen250896/geneCor) to identify a list of CNAcor and METcor genes, respectively. Then, we detected prognostic subgroups for individual CNAcor and METcor datasets, and an integration of CNAcor + METcor + mRNA data through single and integrated analyses using PINPlus, respectively. Finally, we examined the relationships between our classification system and existing subtypes of UM patients before characterizing our two integrated subgroups (i.e., IntSub1 and IntSub2) using the R package [GeneCluster](https://github.com/huynguyen250896/GeneCluster). Abbreviation: UM, Uveal melanoma

Here, we used [Uveal melanoma (UM) data](https://www.cell.com/cancer-cell/fulltext/S1535-6108(17)30295-7#secsectitle0110) from The Cancer Genome Atlas (TCGA).
#### 3. Implementation 
---
**UPDATING...**

#### 4. Contact
---
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) (huynguyen96.dnu@gmail.com) or [Duc-Hau Le](https://github.com/hauldhut) (hauldhut@gmail.com) for any questions about the paper, datasets, code and results.


