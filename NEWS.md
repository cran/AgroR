# AgroR 1.3.6

* fixed the bug when there are more than 10 treatments for `dunnett` function

* The `summarize_conj` function has been implemented. So far, it is limited to summarizing the main analysis of variance framework, in addition to returning the QMres ratio and the means test when the assumptions for jointly analyzing the data are met.

* the `fat2_table` function was implemented to facilitate the generation of dynamic documents for tests in a factorial scheme with two factors or split plot

* The `grid_onefactor` function was implemented to group graphs of DIC, DBC and DQL outputs. It is possible to join up to six graphs in one figure.

* The `correlation` function has been improved

# AgroR 1.3.5

* Fixed test_two function bug for paired data.

* Implemented the option to merge non-parametric and parametric tests in `barplot_positive`

# AgroR 1.3.4

* Fixed the following error: The Scott-Knott test failed when data are transformed for factorial and split-plot schemes.

* As suggested by users, the `names.factor` argument was implemented to change the name of factors in the output of functions that encompass more than one factor

* Fixed function summarise_anova in case of analysis using kruskal-wallis test

# AgroR 1.3.3

* Fixed residual graph bug for `FAT2DIC`, `FAT2DBC`, `PSUBDIC`, `PSUBDBC`, `FAT2DIC.ad`, `FAT2DBC.ad`, `FAT3DIC`, `FAT3DBC`, `FAT3DIC.ad`, `FAT3DBC .ad` when data is unsorted or ordered by block.

* `ylim` argument have been implemented in `dic`, `dbc` and `dql` function.

* fixed bug when `transf` is set to "angular"

* Sum of a `constant` is now available in all experimental projects

# AgroR 1.3.2

* `dic.ad` and `dbc.ad` function bug fixed

* Added to `quantfat2desd` and `jointcluster` functions

* The `conjdic` and `conjdbc` functions have been improved
 
# AgroR 1.3.1

* fix a bug in `barplot_positive`

* Fixed DIC, DBC and DQL function bug when point="mean_se" and factor is qualitative

* Now point in the `FAT2DIC`, `FAT2DBC`, `FAT2DIC.ad`, `FAT2DBC.ad`, `FAT3DIC`, `FAT3DBC`, `FAT3DIC.ad`, `FAT3DBC.ad`, `PSUBDIC` and `PSUBDBC` functions use the point argument ("mean_sd" or "mean_se") to define the information that will be used in the error bars in qualitative factors.

* FAT2DBC.ad function sum of squares and degree of freedom bug fixed

* The `conjfat2dbc` function has been implemented

* The `dic.ad` and `dbc.ad` function has been implemented

* The shiny app has been translated into English

# AgroR 1.3.0

* Fixed `sk_graph` function offset

* Added `ibarplot_double` function

* Added `bargraph2` function

* `plot_jitter` has been improved

* fixed significance level bug for Dunnett's test (`dunnett` function)

# AgroR 1.2.9

 * Fix a bug standard residuals graph in functions `FAT2DIC`, `FAT2DBC`, `FAT2DIC.ad`, `FAT2DBC.ad`, `FAT3DIC`, `FAT3DBC`, `FAT3DIC.ad`, `FAT3DBC.ad`, `PSUBDIC` and `PSUBDBC`
 
 * fix a bug in double interaction significant in `FAT3DIC`,`FAT3DBC`,`FAT3DIC.ad` and `FAT3DIC.ad`

# AgroR 1.2.8

## Minor changes

-   The sorting bug in the case of factorial designs including functions `FAT2DIC`, `FAT2DBC`, `FAT2DIC.ad`, `FAT2DBC.ad`, `FAT3DIC`, `FAT3DBC`, `FAT3DIC.ad`, `FAT3DBC.ad`, `PSUBDIC` and `PSUBDBC` has been fixed. Now the final worksheet, when meaningful interaction doesn't go wrong when the user supplies input parameters in a disorderly fashion.

-   `axissize` argument has been added to the `sketch` function

-   `summarise_dunnett` has been added

# AgroR 1.2.7

-   Add `labelsize` argument in `FAT2DIC`, `FAT2DBC`, `PSUBDIC`, `PSUBDBC`, `FAT3DIC`, `FAT3DBC`, `FAT2DIC.ad`, `FAT2DBC.ad`, `FAT3DIC.ad` and `FAT3DBC.ad`

-   Add `xlab.factor` argument in `FAT2DIC`, `FAT2DBC`, `PSUBDIC`, `PSUBDBC`, `FAT3DIC`, `FAT3DBC`, `FAT2DIC.ad`, `FAT2DBC.ad`, `FAT3DIC.ad` and `FAT3DBC.ad`

-   The output of the `desc3fat` function has been improved

-   The sketch of a strip-plot has been added to the `sketch` function

-   Analysis of variance of the strip-plot scheme has been added in `STRIPLOT`function

-   It is now possible to add the degree of the polynomial for each type of interaction in the case of quantitative factor.

# AgroR 1.2.6

-   The `bargraph_onefactor` and `bargraph_twofactor` functions have been implemented

-   The bug of sorting to `FAT3DBC` in case of significant interaction has been fixed

# AgroR 1.2.5

-   The `FAT2DIC.art` and `FAT2DBC.art` functions have been discontinued

-   Dependencies `ARTool`, `reshape2`, `Hmisc`, `stringr` and `ScottKnott` have been removed

-   Bug fixing `conjdic` and `conjdbc` functions for linux

# AgroR 1.2.4

-   The `.welcome` startup function has been removed

-   The `dunn` function has been added.

-   LSD testing on `FAT2DBC` was added

-   Dependence on the `grid` and `gridExtra` packages have been removed. Functions with more than one graph are now joined by the `cowplot` package.

-   Added `add.letters` argument in `DICT`, `DBCT` and `DQLT` function.

# AgroR 1.2.3

-   Added the `soybean` dataset

-   Added the `bean` dataset

-   Added the `corn` dataset

-   Added the `covercrops` dataset

-   Added the `orchard` dataset

-   Added the `pepper` dataset

-   The summarise_anova function was implemented, which performs a summary of the outputs of the `DIC`, `DBC` and `DQL` functions, when the factor is qualitative. In the case of `FAT2DIC`, `FAT2DBC`, `PSUBDIC` and `PSUBDBC` a summary of the analysis of variance frame is returned.

-   Fixed bug of `conjdic` and `conjdbc` functions in case of analysis of separate experiments.

-   Now for `DIC`, `DBC`, `DQL`, `FAT2DIC`, `FAT2DBC`, `FAT3DIC` and `FAT3DBC` functions, the standardized residual graph is returned in list form, before the graphs.

-   The standardized residual graph for the `PSUBDIC` and `PSUBDBC` functions was implemented

# AgroR 1.2.2

-   Scott-Knott's test has been improved. Now composite names can be used in factor vectors.

-   Fixed the `FAT2DIC.art` function bug

-   `PCA_function` function was implemented

# AgroR 1.2.1

## Minor changes

-   Correction of `CV` for `FAT2DIC` and `FAT2DBC`. Before, I was extracting the QM from the interaction by making the QM from the residue.

-   The `transf` function suggestion issue has been fixed.

-   The `linesize` argument was implemented in `polynomial`, `polynomial2` and `polynomial2_color`

# AgroR 1.2.0

## Major changes

-   Added the `emerge` dataset

-   `TBARPLOT.reverse` function was implemented

-   `FAT2DIC.ad` and `FAT2DBC.ad` function was implemented

-   The `polynomial` and `polynomial2` functions have been improved

-   When a significant double interaction is found when one of the factors is quantitative, it now performs multiple comparisons of the other factor within each level of the quantitative factor.

-   The `PSUBSUB` function bugs have been fixed. In addition, the function returns breakdown of the analysis of variance, assumptions, and additional information.

-   `dunnett` function was implemented

-   `cor_ic` function was implemented

-   `bar_graph` function was implemented

-   `seg_graph` function was implemented

-   The graphical output of the `conjdic` and `conjdbc` functions has been changed. It is now possible to use the outputs of these functions in `sk_graph`, `radargraph` and `barplot_positive`.

-   Now `FAT2DIC`, `FAT2DBC`, `PSUBDIC`, `PSUBDBC`, `FAT3DIC` and `FAT3DBC` can be used in the `sk_graph`, `radargraph` and `barplot_positive` commands, in cases of significant isolated effect.

# AgroR 1.1.0

## Major changes

-   More graphical arguments were added to the commands, such as `angle.label`.

-   The layout of the variance analysis output in the `conjdic` and `conjdbc` commands has been changed.

-   The variation coefficient (CV %) was added in the `DICT`, `DBCT` and `DQLT` functions.

-   Analysis using a generalized linear model of the binomial or poisson family was added for a completely randomized design (`DIC.glm`) in randomized blocks (`DBC.glm`) with a factor of interest.

-   Additional information was inserted in FAT2DIC, FAT2DBC, FAT3DIC and FAT3DBC.

-   `plot_cor` function has been added.

-   `plot_jitter` function has been added.

-   Graphics themes default has been changed from `theme_bw` to `theme_classic`.

-   The graphical output of standardized graphics residues for ggplot2 has been changed.

-   The `aristolochia` dataset has been added

-   The `aacp` function has been added, which calculates the area under the progress curve

## Bug fixes

-   The bar width error in the polynomial and polynomial2 commands has been fixed and can now be defined by the `width.bar` argument.

-   The mean and general median of the DIC, DBC and DQL commands have been corrected. The function returned the measurements of the transformed data.
