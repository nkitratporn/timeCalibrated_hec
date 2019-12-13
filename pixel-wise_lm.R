###########################################################
#       Fit Linear Regression and Export as Rasters       #
###########################################################

#--- Read in rasters from folder ####
result.folder <- c('data/sdm_output/occ-timeDependent/prediction/25-results/')

conflict.wet <- stack(paste0(result.folder, list.files(result.folder, pattern='conflictStack_wet.tif$')))
conflict.dry <- stack(paste0(result.folder, list.files(result.folder, pattern='conflictStack_dry.tif$')))

scenario.phy <- stack(paste0(result.folder, list.files(result.folder, pattern='pred_phyStack.tif$')))
scenario.human <- stack(paste0(result.folder, list.files(result.folder, pattern='pred_humanStack.tif$')))

#--- Change raster name ####
names(conflict.wet) <- c('c_wet_2009','c_wet_2010','c_wet_2011','c_wet_2012','c_wet_2013','c_wet_2014','c_wet_2015','c_wet_2016','c_wet_2017','c_wet_2018')
names(conflict.dry) <- c('c_dry_2009','c_dry_2010','c_dry_2011','c_dry_2012','c_dry_2013','c_dry_2014','c_dry_2015','c_dry_2016','c_dry_2017','c_dry_2018')

names(scenario.phy) <- c('p2009_dry','p2009_wet','p2010_dry','p2010_wet','p2011_dry','p2011_wet','p2012_dry','p2012_wet','p2013_dry','p2013_wet',
                         'p2014_dry','p2014_wet','p2015_dry','p2015_wet','p2016_dry','p2016_wet','p2017_dry','p2017_wet','p2018_dry','p2018_wet')
names(scenario.human) <- c('h2009_dry','h2009_wet','h2010_dry','h2010_wet','h2011_dry','h2011_wet','h2012_dry','h2012_wet','h2013_dry','h2013_wet',
                           'h2014_dry','h2014_wet','h2015_dry','h2015_wet','h2016_dry','h2016_wet','h2017_dry','h2017_wet','h2018_dry','h2018_wet')

#--- Use keywords to filter stack ####

scenario.dry.phy <- dropLayer(scenario.phy, grep(names(scenario.phy),pattern = '*wet')) 
scenario.wet.phy <- dropLayer(scenario.phy, grep(names(scenario.phy),pattern = '*dry'))
scenario.dry.human <- dropLayer(scenario.human, grep(names(scenario.human),pattern = '*wet')) 
scenario.wet.human <- dropLayer(scenario.human, grep(names(scenario.human),pattern = '*dry'))

#--- Perform pixel-wise linear regression and retrieve i) slope and ii) intercept values ####
time <- 1:nlayers(scenario.wet.phy) # set time period

fun.slope = function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[2] }} #fit lm per pixel, only fit time is not NA
fun.intercept = function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[1] }}


wet.phy.slope <- calc(scenario.wet.phy, fun.slope); wet.phy.intercept <- calc(scenario.wet.phy,fun.intercept)
wet.human.slope <- calc(scenario.wet.human, fun.slope); wet.human.intercept <- calc(scenario.wet.human,fun.intercept)
dry.phy.slope <- calc(scenario.dry.phy, fun.slope); dry.phy.intercept <- calc(scenario.dry.phy,fun.intercept)
dry.human.slope <- calc(scenario.dry.human, fun.slope); dry.human.intercept <- calc(scenario.dry.human,fun.intercept)


writeRaster(stack(wet.phy.slope,wet.phy.intercept),'data/sdm_output/occ-timeDependent/prediction/25-results/trend_phy_wet.tif',overwrite=T)
writeRaster(stack(wet.human.slope,wet.human.intercept),'data/sdm_output/occ-timeDependent/prediction/25-results/trend_human_wet.tif',overwrite=T)
writeRaster(stack(dry.phy.slope,dry.phy.intercept),'data/sdm_output/occ-timeDependent/prediction/25-results/trend_phy_dry.tif',overwrite=T)
writeRaster(stack(dry.human.slope,dry.human.intercept),'data/sdm_output/occ-timeDependent/prediction/25-results/trend_human_dry.tif',overwrite=T)