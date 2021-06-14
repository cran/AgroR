# AgroR 1.2.0

## Major changes

* Added the `emerge` dataset 

* `TBARPLOT.reverse` function was implemented

* `FAT2DIC.ad` and `FAT2DBC.ad` function was implemented

* The `polynomial` and `polynomial2` functions have been improved

* When a significant double interaction is found when one of the factors is quantitative, it now performs multiple comparisons of the other factor within each level of the quantitative factor.

* The `PSUBSUB` function bugs have been fixed. In addition, the function returns breakdown of the analysis of variance, assumptions, and additional information. 

* `dunnett` function was implemented 

* `cor_ic` function was implemented 

* `bar_graph` function was implemented 

* `seg_graph` function was implemented 

* The graphical output of the `conjdic` and `conjdbc` functions has been changed. It is now possible to use the outputs of these functions in `sk_graph`, `radargraph` and `barplot_positive`.

* Now `FAT2DIC`, `FAT2DBC`, `PSUBDIC`, `PSUBDBC`, `FAT3DIC` and `FAT3DBC` can be used in the `sk_graph`, `radargraph` and `barplot_positive` commands, in cases of significant isolated effect.

# AgroR 1.1.0

## Major changes

* More graphical arguments were added to the commands, such as `angle.label`.

* The layout of the variance analysis output in the `conjdic` and `conjdbc` commands has been changed.

* The variation coefficient (CV %) was added in the `DICT`, `DBCT` and `DQLT` functions. 

* Analysis using a generalized linear model of the binomial or poisson family was added for a completely randomized design (`DIC.glm`) in randomized blocks (`DBC.glm`) with a factor of interest. 

* Additional information was inserted in FAT2DIC, FAT2DBC, FAT3DIC and FAT3DBC. 

* `plot_cor` function has been added. 

* `plot_jitter` function has been added.

* Graphics themes default has been changed from `theme_bw` to `theme_classic`. 

* The graphical output of standardized graphics residues for ggplot2 has been changed.

* The `aristolochia` dataset has been added 

* The `aacp` function has been added, which calculates the area under the progress curve 

## Bug fixes

* The bar width error in the polynomial and polynomial2 commands has been fixed and can now be defined by the `width.bar` argument.

* The mean and general median of the DIC, DBC and DQL commands have been corrected. The function returned the measurements of the transformed data. 


