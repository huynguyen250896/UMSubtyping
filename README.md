# Re-identification of patient subgroups in Uveal Melanoma
---
This repository contains source code and original/preprocessed datasets for the paper "[Re-identification of patient subgroups in Uveal Melanoma](https://www.frontiersin.org/articles/10.3389/fonc.2021.731548/full)". 
#### 1. Introduction
---
Uveal melanoma (UM) is a comparatively rare cancer but requires serious consideration since patients with developing metastatic UM survive only for about 6-12 months. Fortunately, increasingly large multi-omics databases allow us to further understand cancer initiation and development. Moreover, previous studies have observed that associations between copy number aberrations (CNA) or methylation (MET) versus mRNA expression (mRNA) have affected these processes. From that, we decide to explore the effect of these associations on a case study of UM. Also, the current subtypes of UM display its weak association with biological phenotypes and its lack of therapy suggestions. Therefore, the re-identification of molecular subtypes is a pressing need. In this study, we recruit three omics profiles, including CNA, MET, and mRNA, in a UM cohort from The Cancer Genome Atlas (TCGA). Firstly, we identify two sets of genes, CNAexp and METexp, whose CNA and MET significantly correlated with their corresponding mRNA, respectively. Then, single and integrative analyses of the three data types are performed using the PINSPlus tool. As a result, we discover two novel integrative subgroups, IntSub1 and IntSub2, which could be a useful alternative classification for UM patients in the future. To further explore molecular events behind each subgroup, we identify their subgroup-specific genes computationally. Accordingly, the highest expressed genes among IntSub1-specific genes are mostly enriched with immune-related processes. On the other hand, IntSub2-specific genes are highly associated with cellular cation homeostasis, which responds effectively to chemotherapy using ion channel inhibitor drugs. In addition, we detect that the two integrative subgroups show different age-related risks and survival rates. These discoveries can influence the frequency of metastatic surveillance and support medical practitioners to choose an appropriate treatment regime.

**NOTE:** All statistical analyses were performed using R statistical software.

#### 2. Analysis Pipeline
---
![Figure1](https://imgur.com/aRe4ca7.png)
**Figure 1. Analysis pipeline.** Firstly, we inputted CNA and MET datasets with their corresponding mRNA data to the function [geneCor](https://github.com/huynguyen250896/geneCor) to identify a list of CNAexp and METexp genes, respectively. Then, PINPlus was used to extract different patient subgroups for individual CNAexp and METexp datasets and integration of CNAexp + METexp + mRNA data through single and integrated analyses, respectively. Finally, we discovered subtype-specific genes within each identified integrated subgroups, IntSub1 and IntSub2, using the R package [GeneCluster](https://github.com/huynguyen250896/GeneCluster). Abbreviation: UM, Uveal melanoma.

Here, we used [Uveal melanoma (UM) data](https://www.cell.com/cancer-cell/fulltext/S1535-6108(17)30295-7#secsectitle0110) from The Cancer Genome Atlas (TCGA).

#### 3.Citation
If you found our work interesting, please cite it as follows
```sh
Reference Type: Journal Article
Author: Thi Hai Yen Nguyen
Tin Nguyen
Nguyen, Quang-Huy
Le, Duc-Hau
Year: 2021
Title: Re-identification of patient subgroups in Uveal Melanoma
Journal: Frontiers in Oncology
Date: 2021/10/20
DOI: 10.3389/fonc.2021.731548
```

Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) (huynguyen96.dnu@gmail.com) or [Duc-Hau Le](https://github.com/hauldhut) (hauldhut@gmail.com) for any questions about the paper, datasets, code and results.


