Code used for QTL mapping

document what I did for QTL mapping

```
# first of all, make the gentic linkage map using Lepmap3
using this cal_distance.sh script

# then make the recombination fraction plot
using this make_final_rf_plot.R script

# then add phenotypes to the genotype csv file
using this Jenny_pheno_qtl.R

# finally with the genotype csv file ready, run Rqtl
using this phenoallqtl_plus_em.R script

# summarize the qtl mapping resutls
using this define_qtl_loci.R script
```
