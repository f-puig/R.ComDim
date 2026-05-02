# ComDim

## What's new (April 2026 update)

This package is an updated version of the original R.ComDim package published in 2021 (available at https://github.com/f-puig/R.ComDim). The most important changes introduced in this update are:

* **`MultiBlock()` is now the sole constructor**: `BuildMultiBlock()` has been deprecated and removed. The `MultiBlock()` constructor has been extended with new parameters (`ignore.names`, `ignore.size`) to better handle unbalanced designs, and the S4 accessor API has been overhauled — the old `getBlockNames()` / `setBlockNames()` / `getSampleNames()` / ... style functions are replaced by S4 generics with assignment forms (`blockNames()`, `blockNames<-()`, `sampleNames()`, `sampleNames<-()`, `variableNames()`, `variableNames<-()`) and standard `nrow()` / `ncol()`.
* **`ComDim_OPLS()`**: A new OPLS variant of the ComDim algorithm has been added, complementing the existing `ComDim_PLS()`.

---

## Introduction: the ComDim method
ComDim (also known as CCSWA) is an unsupervised multi-block method that aims to simultaneously consider multiple data tables to find the latent components that are common to all the tables as well as those that are specific to each data table, along with the contribution of each of the tables to each of these components. ComDim determines a common space describing the dispersion of the samples in all the blocks, each block having a specific weight (__salience__) associated with each dimension in this common space. Significant differences in the saliences for a given dimension reflect the fact that the dimension contains different amounts of information coming from each block. In addition to the saliences, __Local loadings__ for each analyzed block and two different sets of scores are obtained. The first set corresponds to the __Local scores__ for each analyzed block while the second set is composed of the __Global scores__, common to all the blocks.

## Why should I use ComDim?
* To analyze __different types of data__ (ex. multi-omics) and understand how the information is distributed across them.
* To extract the common profiles of __related variables__ (ex. metabolites detected in the same pathway).
* To deal with __unbalanced multi-block__ datasets (ex. different number of sample replicates in the blocks). However, ComDim can also deal with __balanced multi-block__ datasets.
* Within the data from the same analytical platform, to evaluate __inter-sample variability__ and __batch effects__ related to the analytical platform.
* To investigate __cross-platform variability__, which is useful to detect errors in the sample preparation.

## Functions
To successfully extract all the potential of the ComDim method, several functions coded in R are proposed. Some of them are listed below:
* __MultiBlock()__: To initialize a MultiBlock object with the first data-block(s). 
* __ExpandMultiBlock()__: To split a single data block into multiple blocks according to a metadata grouping (e.g. metabolic pathways, gene ontology terms), allowing variables to appear in more than one block simultaneously.
* __NormalizeMultiBlock()__: To normalize (some or all) the data-blocks of the MultiBlock object.
* __ProcessMultiBlock()__: To apply customized data transformation to (some or all) the data-blocks of the MultiBlock object.
* __ComDim_PCA()__: This function applies the ComDim-PCA algorithm on a MultiBlock object.

For more information on the usage of these functions, please consult the __tutorial__ from the __docs__ folder.

## Install and load R.ComDim package
```r
  if (!require("devtools")) install.packages("devtools")
  library("devtools")
  install_github("f-puig/R.ComDim")
  
  # Load R.ComDim
  library("R.ComDim")
```

## References
* Puig-Castellví, F.; Jouan-Rimbaud Bouveresse, D.; Mazéas, L.; Chapleur, O.; Rutledge, D. N. Rearrangement of incomplete multi-omics datasets combined with ComDim for evaluating replicate cross-platform variability and batch influence. 	Chemom. Intell. Lab. Syst. 2021, 18 (104422). https://doi.org/10.1016/j.chemolab.2021.104422
* Qannari, E. M.; Courcoux, P.; Vigneau, E. Common Components and Specific Weights Analysis Performed on Preference Data. Food Qual. Prefer. 2001, 12 (5–7), 365–368. https://doi.org/10.1016/S0950-3293(01)00026-X.
* Mazerolles, G.; Hanafi, M.; Dufour, E.; Bertrand, D.; Qannari, E. M. Common Components and Specific Weights Analysis: A Chemometric Method for Dealing with Complexity of Food Products. Chemom. Intell. Lab. Syst. 2006, 81 (1), 41–49. https://doi.org/10.1016/J.CHEMOLAB.2005.09.004.
* Claeys-Bruno, M.; Béal, A.; Rutledge, D. N.; Sergent, M. Use of the Common Components and Specific Weights Analysis to Interpret Supersaturated Designs. Chemom. Intell. Lab. Syst. 2016, 152, 97–106. https://doi.org/10.1016/j.chemolab.2016.01.014.
