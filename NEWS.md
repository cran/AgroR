# AgroR 1.1.0

## Major changes

* More graphical arguments were added to the commands, such as `angle.label`.

* The layout of the variance analysis output in the conjdic and conjdbc commands has been changed.

* The variation coefficient (CV %) was added in the DICT, DBCT and DQLT functions. 

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


