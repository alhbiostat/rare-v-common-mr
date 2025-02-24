# rare-v-common-mr

Comparing the performance of rare and common genetic instruments in an MR framework.

## Background

Rare variant associations with molecular and complex traits are now testable in large-scale population cohorts of whole-exome or whole-genome sequenced individuals. These varaiants typically have differing functional properties to common variants, having been suppressed to low population frequency by negative selective pressures acting against deleterious biological effects. In the context of mendelian randomisation (MR), it may be expected that genetic variants that have larger effect sizes, and tangibly interfere with molecular pathways, have greater utility in instrumenting molecular exposures. These vartiants may, for example, provide more accurate estimates of the causal effect of extreme modulation of a molecular exposure on a outcome (as in therapeutic inhibiton of a drug target).

This projects serves to investigate and compare the performance of rare genetic instruments for molecular traits with their common counterparts. It also explores differing methods of utilising rare genetic instruments for molecular traits, for example by generating gene-based aggregate scores across differing genomic regions (coding, regulatory, protein domain specific) or using gene imparement scores.  

**Two components of the work involve:**
1. Conduting MR of **protein level exposures --> complex traits** using:  
   - coding cis-pQTLs and trans-pQTLs (common, rare and burden and gene impairment instruments)  
   - proximal regulatory pQTLs (e.g promoter region; common, rare and burden and DeepRVAT gene impairment instruments)  
   - distal regulatory pQTLs (e.g. in enhancer regions; common, rare and burden and DeepRVAT gene impairment instruments)
2. Conducting MR of **complex trait --> complex trait**:  
   - comparing single variant (rare and common) and gene-based (burden masks and DeepRVAT gene impairment scores) instruments

DeepRVAT scores are unique in that they allow instrumentation of gene "disruptions" (inferred from sequence) in a way that does not necessitate a known statistical associations with a measured molecular trait (e.g. gene or protein expression). As many features of encoded translated/transcibed transcripts are unmeasured (e.g. protein activity, biding affinity, isoform ratio) this enables instrumentation of genes in a way that does not assume a patricular molecular mode of causality.

## Resources

**Common-variant pQTLs: [Plasma proteomic associations with genetics and health in the UK Bionank (Sun et al.)](https://www.nature.com/articles/s41586-023-06592-6#Sec49)**

* 2,923 plasma proteins measured in 54,219 UK Biobank participants
    * Reported associations with age, sex, BMI,smoking, medication use, liver function, renal function and 20 prevalent diseases
* pQTL discovery: 34,557 EUR participants, 16.1 million imputed variants (from genotypes)
    * 14,287 primary associations (across 3,760 genes) at strict multiple testing (P<1.7X10^-11)
        * 29,420 independent signals with SuSiE regression for fine mapping
    * 66.9% of protiens showing at least one cis associatiion (within 1MB)
    * 92% of cis associaitons within 100kb of the gene start site

**UKB-PPP portal of pQTL summary statistics: https://metabolomips.org/ukbbpgwas/**

**Rare-variant pQTLs: [Rare variant associations with plasma protein levels in the UK Biobank (Dhindsa et al.)](https://www.nature.com/articles/s41586-023-06547-x)**

* 2,923 plasma proteins measured in 49,736 UK Biobank participants
* pQTL discovery: 617,073 variants with MAF down to 0.006% in EUR
    * 5,433 rare genotype-protein associations
    * 1,962 gene-protein associations (aggregate test)

Variant classes (rare = MAF <0.1%):
* *cis-CDS*: **2,227** rare variants (in or nearby given protein)
* *cis-position trans-CDS* (within 1MB of gene encoding altered protein)
* *trans-CDS* (greater than 1MB away)

The absolute effect of both rare **cis** and **trans** pQTLs was significantly larger than common variants

**PheWAS portal for pQTL summary statistics: https://azphewas.com/** \
**pQTL browser: https://astrazeneca-cgr-publications.github.io/pqtl-browser/**

**Rare-variant trait associations: [Exome sequencing analysis of 454,787 UK Bionank participants (Backman et al.)](https://www.nature.com/articles/s41586-021-04103-z)**

**Genebass: [Systematic single-variant and gene-based association testing in thousands of phenotypes in 394,841 UK Biobank exomes (Karczewski et al.)](https://app.genebass.org/)**

## Required datasets

Rare-variant and aggregate pQTL effect sizes: **SNP --> PROTEIN**

1. Significant cis- and trans- CDS pQTLs for olink proteins (with variant classification: PTV, missense, synonymous etc.) - Dhindsa
2. Significant aggregate gene pQTLs (various burden masks) - Dhindsa

Rare-variant and burden effect of trait: **SNP --> TRAIT**

3. Variant effects on UKB trait - Backman
4. Burdan mask effects on UKB trait - Dhindsa 

How to pick low hanging fruit for example MR?: **PROTEIN --> TRAIT**

4. PheWAS protein associations with UKB traits (either plasma protein levels, or burdan masks generarted from variants in olink protein encoding genes)

## So, what now...

We wish to compare MR causal effect estimates for molecular (protein) exposures on health traits using rare and common instruments. Ideally, this requires there to be a detectable causal effect we can compare the magnitude of across instrument sets...

The low hanging fruit to start off with are those UKB olink proteins that show an association with measured traits. This has be tested using the measured protein concentrations themselves (as done for a small number of traits in Sun et al.) and using gene based tests of amalgamated PTVs in olink genes (essentially, leveraging genetically predicted protein levels). This information is available in the PheWAS portal above, so may be the best place to start...

*PTV = protein truncating variant*

