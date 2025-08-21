# UBD
Incorporating uncertainty in cell type proportion estimates from bulk samples to infer cell-type-specific profiles

## Tutorial
### To get the cell-type-specific profiles from bulk profiles given cell type proportion estimates (from BayesPrism), please run 
```
Rscript UBD_refineW.R [file for bulk expression/DNA methylation] [file for BayesPrism output] [path for outputs] \
ref_file=[reference file from signature matrix]
```
- **file for bulk expression/DNA methylation:** This is the bulk expression/DNA methylation matrix. Each row represents an individual and each column represents one gene or one CpG site.
```
            gene1     gene2     gene3    gene4    gene5
sample1  3.272629 1.9255926 1.4394065 2.229085 2.316632
sample2  2.754884 3.0336091 2.1822057 2.512674 3.015411
sample3  3.691906 1.9835263 1.1408063 2.869823 2.934138
sample4  2.978322 2.2354733 0.9303458 2.425820 2.613886
sample5  2.689748 0.9205953 2.4768426 2.640798 2.695839
sample6  2.808944 1.1204922 2.3739149 2.719597 3.417008
...
```
- **file for BayesPrism output:** This is the direct output from BayesPrism, which contains information on sample-level cell type proportion estimates and uncertainty (Chu, T., Wang, Z., Pe’er, D. et al. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. *Nat Cancer* 3, 505–517 (2022). https://doi.org/10.1038/s43018-022-00356-3)
- **path for outputs:** For example, `/Mypath/result`.
- **Optional arguments:** `scale_MH` specifies the proposal scaling in Metropolis-Hastings (MH). It adjusts the size (scale) of the proposal distribution's steps and also controls the acceptance rate. The default value is `scale_MH=200`. `ref_file` specifies xxx. We strongly recommend to provide xxx. If no reference information is available, we will estimate it within the UBD algorithm. In `ref_file`, each row corresponds to one gene/CpG site and each column corresponds to one cell type.
```
            CT1      CT2      CT3
gene1 3.4560016 2.672303 3.360823
gene2 0.9288779 1.620721 1.575300
gene3 1.6935820 1.447028 1.284312
gene4 2.7038122 2.566223 2.446321
gene5 3.7489763 3.098162 3.192967
gene6 1.6620560 1.144611 1.520568
...
```

Using the example inputs, a typical run is like
```
Rscript UBD_refineW.R bulkM.RData bp.res.RData ./result \
ref_file=alpha_hat.RData 
```

