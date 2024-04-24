
# [AgroR](https://agronomiar.github.io/AgroR_package/index.html) <img src='logo.png' align="right" height="139" />

<!-- badges: start -->

![Downloads](http://cranlogs.r-pkg.org/badges/AgroR)
[![CRAN status](https://www.r-pkg.org/badges/version-ago/AgroR)](https://CRAN.R-project.org/package=AgroR)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable-1)
![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/AgroR)

<!-- badges: end -->

Package: AgroR

Type: Package

Title: Experimental Statistics and Graphics for Agricultural Sciences 

Version: 1.3.6

Date: 2023-26-12

Authors: 
 
 - Gabriel Danilo Shimizu
 - Rodrigo Yudi Palhaci Marubayashi
 - Leandro Simoes Azeredo Goncalves

Maintainer: Gabriel Danilo Shimizu <shimizu@uel.br>

Description: Performs the analysis of completely randomized experimental designs (CRD), randomized blocks (RBD) and Latin square (LSD), experiments in double and triple factorial scheme (in CRD and RBD), experiments in subdivided plot scheme (in CRD and RBD), subdivided and joint analysis of experiments in CRD and RBD, linear regression analysis, test for two samples. The package performs analysis of variance, ANOVA assumptions and multiple comparison test of means or regression, according to Pimentel-Gomes (2009, ISBN: 978-85-7133-055-9), nonparametric test (Conover, 1999, ISBN: 0471160687), test for two samples, joint analysis of experiments according to Ferreira (2018, ISBN: 978-85-7269-566-4) and generalized linear model (glm) for binomial and Poisson family in CRD and RBD (Carvalho, FJ (2019), <doi: 10.14393/ufu.te.2019.1244>). It can also be used to obtain descriptive measures and graphics, in addition to correlations and creative graphics used in agricultural sciences (Agronomy, Zootechnics, Food Science and related areas).

Encoding: UTF-8

RoxygenNote: 7.1.1

Imports: ggplot2, nortest, lme4, crayon, lmtest, emmeans, multcomp, ggrepel, MASS, cowplot, multcompView, RColorBrewer, drc, dunn.test, gtools

Suggests: DT, knitr, rmarkdown, roxygen2

Depends: R (>= 3.6.0)

License: GPL (>= 2)

# Installation

<br><br>

``` r
# Install release version from CRAN
install.packages("AgroR")

# Install development version from GitHub
devtools::install_github("https://github.com/AgronomiaR/AgroR.git")
```

<br><br>

# References

## Data set

 - `aristolochia`: Germination of seeds of _Aristolochia_ sp. as a function of temperature.
 - `bean`: Bean data
 - `cloro`: Sodium dichloroisocyanurate in soybean
 - `corn`: Corn data
 - `covercrops`: Covercrops data
 - `emerg`: Emergence of passion fruit seeds over time .
 - `enxofre`: Sulfur data
 - `laranja`: Orange plants under different rootstocks
 - `mirtilo`: Cutting blueberry data
 - `orchard`: Orchard data
 - `passiflora`: Substrate data in the production of passion fruit seedlings
 - `pepper`: Pepper data
 - `phao`: Osmocote in *Phalaenopsis* sp.
 - `pomegranate`: Pomegranate data
 - `porco`: Pig development and production
 - `sensorial`: Sensorial data
 - `simulate1`: Simulated data dict
 - `simulate2`: Simulated data dbct
 - `simulate3`: Simulated data dqlt
 - `soybean`: Soybean data
 - `tomate`: Tomato data
 - `weather`: Weather data

## Descritive analysis

 - `desc`: Descriptive analysis
 - `desc2fat`: Descriptive analysis (Two factors)
 - `desc3fat`: Descriptive analysis (Three factors)
 - `dispvar`: Boxplot with standardized data
 - `tabledesc`: Table descritive analysis

## Analysis function

*t test to compare means with a reference value*

 - `tonetest`: t test to compare means with a reference value
 
*Analysis for testing of two independent or dependent samples by parametric or non-parametric method*

 - `test_two`: Test for two samples

*Analysis of simple experiments*

 - `DIC`: Completely randomized design
 - `DBC`: Randomized block design
 - `DQL`: Latin square design

*Analysis of simple experiments with an additional treatment for quantitative factor*

 - `dic.ad`: Completely randomized design with an additional treatment for quantitative factor
 - `dbc.ad`: Randomized block design with an additional treatment for quantitative factor

*Analysis of simple experiments in DIC and DBC by generalized linear model (Binomial or Poisson)*

 - `DIC.glm`: Completely randomized design by glm
 - `DBC.glm`: Randomized block design by glm

*Analysis of experiments in DIC, DBC or DQL with multiple assessments over time or disregarding the effect of another factor*

 - `DICT`: Completely randomized design evaluated over time
 - `DBCT`: Randomized block design evaluated over time
 - `DQLT`: Latin square design evaluated over time
 
*Analysis of groups of experiments in DIC and DBC*

 - `conjdbc`: Joint analysis of experiments in randomized block design
 - `conjdic`: Joint analysis of experiments in completely randomized design
 
*Analysis of groups of experiments in FAT2DBC*

 - `conjfat2dbc`: Joint analysis of experiments in randomized block design in double factorial

*Analysis of experiments in double factorial design in DIC and DBC*

 - `FAT2DIC`: DIC experiments in double factorial
 - `FAT2DBC`: DBC experiments in double factorial
 
*Analysis of double factorial design experiments in DIC or DBC with an additional treatment*

 - `FAT2DIC.ad`: DIC experiment in double factorial design with an additional treatment
 - `FAT2DBC.ad`: DBC experiment in double factorial design with an additional treatment
 
*Analysis of DIC or DBC experiments in a factorial scheme with three factors*

 - `FAT3DIC`: DIC experiments in triple factorial
 - `FAT3DBC`: DBC experiments in triple factorial

*Analysis of triple factorial design experiments in DIC or DBC with an additional treatment*

 - `FAT3DIC.ad`: DIC experiment in double factorial design with an additional treatment
 - `FAT3DBC.ad`: DBC experiment in double factorial design with an additional treatment
 
*Split-plot scheme in DIC or DBC*

 - `PSUBDBC`: DBC experiments in split-plot
 - `PSUBDIC`: DIC experiments in split-plot

*Splitsplitplot parcels scheme in DBC*

 - `PSUBSUBDBC`: DBC experiments in split-split-plot

*Plot subdivided into randomized blocks with a subplot in a double factorial scheme*

 - `PSUBFAT2DBC`: Plot subdivided into randomized blocks with a subplot in a double factorial scheme

*Dunnett's Test for Comparison of Control vs. Treatments*
 
 - `dunnett`: Dunnett test

*Dunn's non-parametric test*

- `dunn`: Dunn test

*Logistic regression 3 or 4 parameters*

 - `logistic`: Logistic regression
 
*Polynomial Regression to the Third Degree*

 - `polynomial`: Linear regression graph
 - `polynomial2`: Linear regression graph in double factorial
 - `polynomial2_color`: Linear regression graph in double factorial with color graph

*Principal component analysis*

 - `PCA_function`: Principal components analysis
 
## Graphs

 - `barfacet`: Bar graph for one factor with facets
 - `bargraph_onefactor`: Group DIC, DBC and DQL functions column charts
 - `bargraph_twofactor`: Group FAT2DIC, FAT2DBC, PSUBDIC or PSUBDBC functions column charts
 - `barplot_positive`: Positive barplot
 - `bar_dunnett`: Barplot for Dunnett test
 - `bar_graph`: Bar graph for one factor
 - `bar_graph2`: Bar graph for one factor model 2
 - `corgraph`: Correlogram
 - `cor_ic`: Plot Pearson correlation with interval of confidence
 - `ibarplot.double`: Invert letters for two factor chart
 - `line_plot`: Line chart
 - `plot_cor`: Plot correlation
 - `plot_interaction`: Interaction plot
 - `plot_jitter`: Column, box or segment chart with observations
 - `plot_TH`: Climate chart of temperature and humidity
 - `plot_TH1`: Climate chart of temperature and humidity (Model 2)
 - `radargraph`: Circular column chart
 - `seg_graph`: Segment graph for one factor
 - `seg_graph2`: Point graph for one factor model 2
 - `sk_graph`: Scott-Knott graphics
 - `spider_graph`: Spider graph for sensorial analysis
 - `TBARPLOT.reverse`: Reverse graph of DICT, DBCT and DQL output when geom="bar"
 - `plot_tonetest`: Graphic for t test to compare means with a reference value

## Utils

 - `sketch`: Experimental sketch
 - `aacp`: Area under the curve
 - `transf`: Data transformation (Box-Cox, 1964)
 - `summarise_anova`: Summary of analysis of variance and test of means
 - `summarise_dunnett`: Dunnett's Test Summary
 - `confinterval`: Interval of confidence for groups
