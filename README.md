# timeCalibrated_hec
Negative human- wild elephant interaction, especially in the form of crop depredation, impacted communities' livelihood and deteriorate conservation succes. Distribution of this interaction across landscape remain limited. Here, we attempted to predict distribution of conflict using occurrence data from news report and time-calibrated Maximum Enterpy method (MaxEnt). Two better identify ultimate drivers of distribution changes overtime, we model two main drivers saperately: resource availability and direct human pressure. Two-dimension conflict matrix was then applied to classify different probability of conflict into *High*, *Likely*, *Low*, and *Rare* conflict category.  

![](/0-graphical_abstract.png)

This project is divided into 3 parts
1. data preparation
  + random-sampling of occurrence dataset
  + raster pre-process (partially after-GEE)
2. modeling
  + extract predictors value across multiple seasons
  + evaluate occurence data: wilcox-test, multicollinearity
  + parameter optimization (EMNeval) and MaxEnt (dismo) modeling
3. analysis
  + conflict classification and hotspot
  + trend
