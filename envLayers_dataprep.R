# library
# spaital related
library("raster")
library("rgeos")
library("rgdal")
library("sp")
#########################
# set working folder if necessary
setwd('D:/Academic/Research/0_maximus-materials/analysis/')

###### EVI ######
#--- load monthly 2008-2018 and remove 2008 data
evi.monthly <- stack('data/sdm_data/data_prep/monthly/EVI_month.tif')
names(evi.monthly) <- c(rep(1:12,11)) # rename layer by month
#--- load 2014-2018 median by month data
evi.median10yrs <- stack('data/sdm_data/data_prep/monthly/EVI_month_median.tif')
names(evi.median10yrs) <- c('evi1','evi2','evi3','evi4','evi5','evi6','evi7','evi8','evi9','evi10','evi11','evi12') # rename

#--- subset stack into yearly 
evi08 <- subset(evi.monthly,1:12)
evi09 <- subset(evi.monthly,13:24)
evi10 <- subset(evi.monthly,25:36)
evi11 <- subset(evi.monthly,37:48)
evi12 <- subset(evi.monthly,49:60)
evi13 <- subset(evi.monthly,61:72)
evi14 <- subset(evi.monthly,73:84)
evi15 <- subset(evi.monthly,85:96)
evi16 <- subset(evi.monthly,97:108)
evi17 <- subset(evi.monthly,109:120)
evi18 <- subset(evi.monthly,121:132)

#--- loop through by month and replace missing value usnig median over 10 yrs
for(i in 1:12){
  evi08[[i]] <- cover(evi08[[i]],evi.median10yrs[[i]])
  evi09[[i]] <- cover(evi09[[i]],evi.median10yrs[[i]])
  evi10[[i]] <- cover(evi10[[i]],evi.median10yrs[[i]])
  evi11[[i]] <- cover(evi11[[i]],evi.median10yrs[[i]])
  evi12[[i]] <- cover(evi12[[i]],evi.median10yrs[[i]])
  evi13[[i]] <- cover(evi13[[i]],evi.median10yrs[[i]])
  evi14[[i]] <- cover(evi14[[i]],evi.median10yrs[[i]])
  evi15[[i]] <- cover(evi15[[i]],evi.median10yrs[[i]])
  evi16[[i]] <- cover(evi16[[i]],evi.median10yrs[[i]])
  evi17[[i]] <- cover(evi17[[i]],evi.median10yrs[[i]])
  evi18[[i]] <- cover(evi18[[i]],evi.median10yrs[[i]])
}

#-- create only dry season 2009 to wet season 2018
evi.cover <- stack(evi09,evi10,evi11,evi12,evi13,evi14,evi15,evi16,evi17,evi18)
names(evi.cover) <- c(rep(1:12,11))
evi.cover <- dropLayer(evi.cover,1:10)
evi.cover <- dropLayer(evi.cover,121:122)
evi.monthly 
writeRaster(evi.cover,'data/sdm_data/data_prep/season_by_year/evi_filled_gap1108-1018.tif')
evi.mean <- stack(evi.monthly[[1]])
rm(evi.08,evi.09,evi.10,evi.11,evi.12,evi.13,evi.14,evi.15,evi.16,evi.17,evi.18)

#--- separate by season every 6 month and calculate standard deviation value
evi.sd <- stack(evi.monthly[[1]])
for(i in 1:20){
  e <- 6*i; s <- e-5
  evi.sd[[i]] <- calc(evi.cover[[s:e]],na.rm=T,sd)
  cat('step:',i,'\n')
}

names(evi.sd) <- c('evi_dry_sd1','evi_wet_sd2','evi_dry_sd3','evi_wet_sd4','evi_dry_sd5','evi_wet_sd6','evi_dry_sd7','evi_wet_sd8','evi_dry_sd9','evi_wet_sd10',
                   'evi_dry_sd11','evi_wet_sd12','evi_dry_sd13','evi_wet_sd14','evi_dry_sd15','evi_wet_sd16','evi_dry_sd17','evi_wet_sd18','evi_dry_sd19','evi_wet_sd20')

#--- separate by season every 6 month and calculate mean value
for(i in 1:20){
  e <- 6*i; s <- e-5
  evi.mean[[i]] <- mean(evi.cover[[s:e]],na.rm=T)
  cat('step:',i,'\n')
}

names(evi.mean) <- c('evi_dry1','evi_wet2','evi_dry3','evi_wet4','evi_dry5','evi_wet6','evi_dry7','evi_wet8','evi_dry9','evi_wet10',
                     'evi_dry11','evi_wet12','evi_dry13','evi_wet14','evi_dry15','evi_wet16','evi_dry17','evi_wet18','evi_dry19','evi_wet20')

#--- calculate heterogeneity using glcm texture for 500m, 1000m, 1500m
library(glcm)
evi.het3 <- stack(evi.mean[[1]]);evi.het5 <- stack(evi.mean[[1]]);evi.het7 <- stack(evi.mean[[1]])

for(i in seq(nlayers(evi.mean))){
    evi.het3[[i]] <- glcm(evi.mean[[i]],shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),window=c(3,3),statistics = c('homogeneity'))
    evi.het5[[i]] <- glcm(evi.mean[[i]],shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),window=c(5,5),statistics = c('homogeneity'))
    evi.het7[[i]] <- glcm(evi.mean[[i]],shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),window=c(7,7),statistics = c('homogeneity'))
}

names(evi.het3) <- c('evi_dry_3het1','evi_wet_3het2','evi_dry_3het3','evi_wet_3het4','evi_dry_3het5','evi_wet_3het6','evi_dry_3het7','evi_wet_3het8','evi_dry_3het9','evi_wet_3het10',
                     'evi_dry_3het11','evi_wet_3het12','evi_dry_3het13','evi_wet_3het14','evi_dry_3het15','evi_wet_3het16','evi_dry_3het17','evi_wet_3het18','evi_dry_3het19','evi_wet_3het20')
names(evi.het5) <- c('evi_dry_5het1','evi_wet_5het2','evi_dry_5het3','evi_wet_5het4','evi_dry_5het5','evi_wet_5het6','evi_dry_5het7','evi_wet_5het8','evi_dry_5het9','evi_wet_5het10',
                     'evi_dry_5het11','evi_wet_5het12','evi_dry_5het13','evi_wet_5het14','evi_dry_5het15','evi_wet_5het16','evi_dry_5het17','evi_wet_5het18','evi_dry_5het19','evi_wet_5het20')
names(evi.het7) <- c('evi_dry_7het1','evi_wet_7het2','evi_dry_7het3','evi_wet_7het4','evi_dry_7het5','evi_wet_7het6','evi_dry_7het7','evi_wet_7het8','evi_dry_7het9','evi_wet_7het10',
                     'evi_dry_7het11','evi_wet_7het12','evi_dry_7het13','evi_wet_7het14','evi_dry_7het15','evi_wet_7het16','evi_dry_7het17','evi_wet_7het18','evi_dry_7het19','evi_wet_7het20')

#--- EVI Slope
evi.slope <- stack('data/sdm_data/data_prep/season_by_year/evi_slope_dry2008-wet2018.tif')
names(evi.slope) <- c('evi_dry_slope1','evi_wet_slope2','evi_dry_slope3','evi_wet_slope4','evi_dry_slope5','evi_wet_slope6','evi_dry_slope7','evi_wet_slope8','evi_dry_slope9','evi_wet_slope10',
                      'evi_dry_slope11','evi_wet_slope12','evi_dry_slope13','evi_wet_slope14','evi_dry_slope15','evi_wet_slope16','evi_dry_slope17','evi_wet_slope18','evi_dry_slope19','evi_wet_slope20')

###### KBDI #####
kbdi.mean <- stack('data/sdm_data/data_prep/season_by_year/kbdi_dry2008-wet2018.tif')
names(kbdi.mean) <- c('kbdi_dry1','kbdi_wet2','kbdi_dry3','kbdi_wet4','kbdi_dry5','kbdi_wet6','kbdi_dry7','kbdi_wet8','kbdi_dry9','kbdi_wet10',
                      'kbdi_dry11','kbdi_wet12','kbdi_dry13','kbdi_wet14','kbdi_dry15','kbdi_wet16','kbdi_dry17','kbdi_wet18','kbdi_dry19','kbdi_wet20')

###### NTL #####
ntl <- stack('data/sdm_data/data_prep/season_by_year/ntl_stack2009-2019.tif')
lit <- ntl
lit[lit < 20] <- 0
lit[lit != 0] <- 1
names(lit) <- c('litup2007','litup2008','litup2009','litup2010','litup2011','litup2012','litup.2013','litup.2014',
                'litup.2015','litup.2016','litup.2017','litup.2018','litup.2019')

### create distance for yearly
disLitArea <- lit[[3:12]]
for (i in 1:nlayers(disLitArea)) {
  disLitArea[[i]] <- gridDistance(lit[[i]], origin = 1)
}

disLitArea <- stack(stack(replicate(2,disLitArea[[1]])),stack(replicate(2,disLitArea[[2]])),stack(replicate(2,disLitArea[[3]])),stack(replicate(2,disLitArea[[4]])),
                    stack(replicate(2,disLitArea[[5]])),stack(replicate(2,disLitArea[[6]])),stack(replicate(2,disLitArea[[7]])),stack(replicate(2,disLitArea[[8]])),
                    stack(replicate(2,disLitArea[[9]])),stack(replicate(2,disLitArea[[10]])))
names(disLitArea) <- c('disLitArea_dry1','disLitArea_wet2','disLitArea_dry3','disLitArea_wet4','disLitArea_dry5','disLitArea_wet6','disLitArea_dry7','disLitArea_wet8','disLitArea_dry9','disLitArea_wet10',
                       'disLitArea_dry11','disLitArea_wet12','disLitArea_dry13','disLitArea_wet14','disLitArea_dry15','disLitArea_wet16','disLitArea_dry17','disLitArea_wet18','disLitArea_dry19','disLitArea_wet20')
###### Water #####
water <- stack('data/sdm_data/data_prep/season_by_year/disWaterStack_3month.tif')
water <- stack(stack(replicate(2,water[[1]])),stack(replicate(2,water[[2]])),stack(replicate(2,water[[3]])),stack(replicate(2,water[[4]])),
               stack(replicate(2,water[[5]])),stack(replicate(2,water[[6]])),stack(replicate(2,water[[7]])),stack(replicate(2,water[[8]])),
               stack(replicate(2,water[[9]])),stack(replicate(2,water[[10]])))
names(water) <- c('water_dry1','water_wet2','water_dry3','water_wet4','water_dry5','water_wet6','water_dry7','water_wet8','water_dry9','water_wet10',
                  'water_dry11','water_wet12','water_dry13','water_wet14','water_dry15','water_wet16','water_dry17','water_wet18','water_dry19','water_wet20')

###### human density ######
humanDen <- stack('data/sdm_data/data_prep/season_by_year/popDen_2009-2017.tif')
humanDen <- stack(humanDen,humanDen[[nlayers(humanDen)]])
humanDen <- stack(stack(replicate(2,humanDen[[1]])),stack(replicate(2,humanDen[[2]])),stack(replicate(2,humanDen[[3]])),stack(replicate(2,humanDen[[4]])),
                  stack(replicate(2,humanDen[[5]])),stack(replicate(2,humanDen[[6]])),stack(replicate(2,humanDen[[7]])),stack(replicate(2,humanDen[[8]])),
                  stack(replicate(2,humanDen[[9]])),stack(replicate(2,humanDen[[10]])))
names(humanDen) <- c('humanDen_dry1','humanDen_wet2','humanDen_dry3','humanDen_wet4','humanDen_dry5','humanDen_wet6','humanDen_dry7','humanDen_wet8','humanDen_dry9','humanDen_wet10',
                     'humanDen_dry11','humanDen_wet12','humanDen_dry13','humanDen_wet14','humanDen_dry15','humanDen_wet16','humanDen_dry17','humanDen_wet18','humanDen_dry19','humanDen_wet20')


###### Protected Natural Habitat #####
pa <- stack('data/sdm_data/data_prep/variables_oct2019/resources_related/disRefuge.tif')
pa <- stack(replicate(20,pa))
names(pa) <- c('pa_dry1','pa_wet2','pa_dry3','pa_wet4','pa_dry5','pa_wet6','pa_dry7','pa_wet8','pa_dry9','pa_wet10',
               'pa_dry11','pa_wet12','pa_dry13','pa_wet14','pa_dry15','pa_wet16','pa_dry17','pa_wet18','pa_dry19','pa_wet20')

###### Main Road #####
road <- stack('data/sdm_data/data_prep/variables_oct2019/human_pressure/disMainRd.tif')
road <- stack(replicate(20,road))
names(road) <- c('road_dry1','road_wet2','road_dry3','road_wet4','road_dry5','road_wet6','road_dry7','road_wet8','road_dry9','road_wet10',
                 'road_dry11','road_wet12','road_dry13','road_wet14','road_dry15','road_wet16','road_dry17','road_wet18','road_dry19','road_wet20')

###### land cover ######
#--- readin data for each year 2009-2018
forest <- stack('data/sdm_data/data_prep/season_by_year/forestStack.tif')
woody <- stack('data/sdm_data/data_prep/season_by_year/woodSavStack.tif')
savanna <- stack('data/sdm_data/data_prep/season_by_year/savStack.tif')
shrub <- stack('data/sdm_data/data_prep/season_by_year/shrubStack.tif')
grass <- stack('data/sdm_data/data_prep/season_by_year/grassStack.tif')
crop <- stack('data/sdm_data/data_prep/season_by_year/cropStack.tif')
unused <- stack('data/sdm_data/data_prep/season_by_year/noUsedStack.tif')

#--- replicate 2 layers per year (one for each season)
forest <- stack(stack(replicate(2,forest[[1]])),stack(replicate(2,forest[[2]])),stack(replicate(2,forest[[3]])),stack(replicate(2,forest[[4]])),
                stack(replicate(2,forest[[5]])),stack(replicate(2,forest[[6]])),stack(replicate(2,forest[[7]])),stack(replicate(2,forest[[8]])),
                stack(replicate(2,forest[[9]])),stack(replicate(2,forest[[10]])))
woody <- stack(stack(replicate(2,woody[[1]])),stack(replicate(2,woody[[2]])),stack(replicate(2,woody[[3]])),stack(replicate(2,woody[[4]])),
               stack(replicate(2,woody[[5]])),stack(replicate(2,woody[[6]])),stack(replicate(2,woody[[7]])),stack(replicate(2,woody[[8]])),
               stack(replicate(2,woody[[9]])),stack(replicate(2,woody[[10]])))
savanna <- stack(stack(replicate(2,savanna[[1]])),stack(replicate(2,savanna[[2]])),stack(replicate(2,savanna[[3]])),stack(replicate(2,savanna[[4]])),
                 stack(replicate(2,savanna[[5]])),stack(replicate(2,savanna[[6]])),stack(replicate(2,savanna[[7]])),
                 stack(replicate(2,savanna[[8]])),stack(replicate(2,savanna[[9]])),stack(replicate(2,savanna[[10]])))
shrub <- stack(stack(replicate(2,shrub[[1]])),stack(replicate(2,shrub[[2]])),stack(replicate(2,shrub[[3]])),stack(replicate(2,shrub[[4]])),
               stack(replicate(2,shrub[[5]])),stack(replicate(2,shrub[[6]])),stack(replicate(2,shrub[[7]])),
               stack(replicate(2,shrub[[8]])),stack(replicate(2,shrub[[9]])),stack(replicate(2,shrub[[10]])))
grass <- stack(stack(replicate(2,grass[[1]])),stack(replicate(2,grass[[2]])),stack(replicate(2,grass[[3]])),stack(replicate(2,grass[[4]])),
               stack(replicate(2,grass[[5]])),stack(replicate(2,grass[[6]])),stack(replicate(2,grass[[7]])),
               stack(replicate(2,grass[[8]])),stack(replicate(2,grass[[9]])),stack(replicate(2,grass[[10]])))
crop <- stack(stack(replicate(2,crop[[1]])),stack(replicate(2,crop[[2]])),stack(replicate(2,crop[[3]])),stack(replicate(2,crop[[4]])),
              stack(replicate(2,crop[[5]])),stack(replicate(2,crop[[6]])),stack(replicate(2,crop[[7]])),
              stack(replicate(2,crop[[8]])),stack(replicate(2,crop[[9]])),stack(replicate(2,crop[[10]])))
unused <- stack(stack(replicate(2,unused[[1]])),stack(replicate(2,unused[[2]])),stack(replicate(2,unused[[3]])),stack(replicate(2,unused[[4]])),
                stack(replicate(2,unused[[5]])),stack(replicate(2,unused[[6]])),stack(replicate(2,unused[[7]])),
                stack(replicate(2,unused[[8]])),stack(replicate(2,unused[[9]])),stack(replicate(2,unused[[10]])))
names(forest) <- c('forest_dry1','forest_wet2','forest_dry3','forest_wet4','forest_dry5','forest_wet6','forest_dry7','forest_wet8','forest_dry9','forest_wet10',
                   'forest_dry11','forest_wet12','forest_dry13','forest_wet14','forest_dry15','forest_wet16','forest_dry17','forest_wet18','forest_dry19','forest_wet20')
names(woody) <- c('woody_dry1','woody_wet2','woody_dry3','woody_wet4','woody_dry5','woody_wet6','woody_dry7','woody_wet8','woody_dry9','woody_wet10',
                  'woody_dry11','woody_wet12','woody_dry13','woody_wet14','woody_dry15','woody_wet16','woody_dry17','woody_wet18','woody_dry19','woody_wet20')
names(savanna) <- c('savanna_dry1','savanna_wet2','savanna_dry3','savanna_wet4','savanna_dry5','savanna_wet6','savanna_dry7','savanna_wet8','savanna_dry9','savanna_wet10',
                    'savanna_dry11','savanna_wet12','savanna_dry13','savanna_wet14','savanna_dry15','savanna_wet16','savanna_dry17','savanna_wet18','savanna_dry19','savanna_wet20')
names(shrub) <- c('shrub_dry1','shrub_wet2','shrub_dry3','shrub_wet4','shrub_dry5','shrub_wet6','shrub_dry7','shrub_wet8','shrub_dry9','shrub_wet10',
                  'shrub_dry11','shrub_wet12','shrub_dry13','shrub_wet14','shrub_dry15','shrub_wet16','shrub_dry17','shrub_wet18','shrub_dry19','shrub_wet20')
names(grass) <- c('grass_dry1','grass_wet2','grass_dry3','grass_wet4','grass_dry5','grass_wet6','grass_dry7','grass_wet8','grass_dry9','grass_wet10',
                  'grass_dry11','grass_wet12','grass_dry13','grass_wet14','grass_dry15','grass_wet16','grass_dry17','grass_wet18','grass_dry19','grass_wet20')
names(crop) <- c('crop_dry1','crop_wet2','crop_dry3','crop_wet4','crop_dry5','crop_wet6','crop_dry7','crop_wet8','crop_dry9','crop_wet10',
                 'crop_dry11','crop_wet12','crop_dry13','crop_wet14','crop_dry15','crop_wet16','crop_dry17','crop_wet18','crop_dry19','crop_wet20')
names(unused) <- c('unused_dry1','unused_wet2','unused_dry3','unused_wet4','unused_dry5','unused_wet6','unused_dry7','unused_wet8','unused_dry9','unused_wet10',
                   'unused_dry11','unused_wet12','unused_dry13','unused_wet14','unused_dry15','unused_wet16','unused_dry17','unused_wet18','unused_dry19','unused_wet20')
#--- % cover of forest
library(raster)
lc <- stack('data/sdm_data/data_prep/season_by_year/lcStack.tif')
lc.forest <- lc
lc.forest[lc.forest!=1] <- 0

lc.savanna <- lc
lc.savanna[lc.savanna != 4] <- 0

lc.forest <- mask(lc.forest,r.studyArea)
lc.savanna <- mask(lc.savanna,r.studyArea)

perc.forest <- stack()
for (i in 1:nlayers(lc.forest)) {
  radius <-6000
  
  fw.i <- ifelse(focalWeight(lc.forest[[i]], radius,'circle')>0,1,0)
  
  focalimg <- focal(lc.forest[[i]]==1,fw.i,na.rm=T)
  focalPortion <- focalimg/sum(fw.i)*100
  perc.forest <- stack(perc.forest,focalPortion)
  rm(focalimg,focalPortion)
}
perc.forest <- stack(stack(replicate(2,perc.forest[[1]])),stack(replicate(2,perc.forest[[2]])),stack(replicate(2,perc.forest[[3]])),
                     stack(replicate(2,perc.forest[[4]])),stack(replicate(2,perc.forest[[5]])),stack(replicate(2,perc.forest[[6]])),
                     stack(replicate(2,perc.forest[[7]])),stack(replicate(2,perc.forest[[8]])),stack(replicate(2,perc.forest[[9]])),
                     stack(replicate(2,perc.forest[[10]])))
names(perc.forest) <- c('perc.forest_dry1','perc.forest_wet2','perc.forest_dry3','perc.forest_wet4','perc.forest_dry5','perc.forest_wet6',
                        'perc.forest_dry7','perc.forest_wet8','perc.forest_dry9','perc.forest_wet10','perc.forest_dry11','perc.forest_wet12',
                        'perc.forest_dry13','perc.forest_wet14','perc.forest_dry15','perc.forest_wet16','perc.forest_dry17','perc.forest_wet18',
                        'perc.forest_dry19','perc.forest_wet20')


##### TRI #####
tri <- raster('data/sdm_data/data_prep/season_by_year/TRI.tif')
tri <- stack(replicate(20,tri))
names(tri) <- c('tri_dry1','tri_wet2','tri_dry3','tri_wet4','tri_dry5','tri_wet6','tri_dry7','tri_wet8','tri_dry9','tri_wet10',
                'tri_dry11','tri_wet12','tri_dry13','tri_wet14','tri_dry15','tri_wet16','tri_dry17','tri_wet18','tri_dry19','tri_wet20')


##### phenology #####
install.packages("greenbrown", repos="http://R-Forge.R-project.org")
library(greenbrown)

##### SAVE BY LAYER TO EACH FOLDER
dir.create('data/sdm_data/data_time-dependent/predictors')
var.folder <- c('data/sdm_data/data_time-dependent/predictors/')
for(i in 1:nlayers(water)){
  dir.create(paste0(var.folder,i))
  cat('start:',i,'/',nlayers(water),'... ')
  folder.name <- paste0(var.folder,i,'/')
  season <- ifelse(i%%2==1,'dry','wet')
    writeRaster(perc.forest[[i]],paste0(folder.name,'percForest_',season,i,'.tif'),overwrite=T)
    writeRaster(evi.mean[[i]],paste0(folder.name,'EVI_',season,i,'.tif'),overwrite=T)
    writeRaster(evi.slope[[i]],paste0(folder.name,'EVI_',season,'_slope',i,'.tif'),overwrite=T)
    writeRaster(evi.het3[[i]],paste0(folder.name,'EVI_',season,'_het3',i,'.tif'),overwrite=T)
    writeRaster(evi.het5[[i]],paste0(folder.name,'EVI_',season,'_het5',i,'.tif'),overwrite=T)
    writeRaster(evi.het7[[i]],paste0(folder.name,'EVI_',season,'_het7',i,'.tif'),overwrite=T)
    writeRaster(kbdi.mean[[i]],paste0(folder.name,'KBDI_',season,i,'.tif'),overwrite=T)
    writeRaster(forest[[i]],paste0(folder.name,'forest_',season,i,'.tif'),overwrite=T)
    writeRaster(woody[[i]],paste0(folder.name,'woody_',season,i,'.tif'),overwrite=T)
    writeRaster(savanna[[i]],paste0(folder.name,'sav_',season,i,'.tif'),overwrite=T)
    writeRaster(shrub[[i]],paste0(folder.name,'shrub_',season,i,'.tif'),overwrite=T)
    writeRaster(grass[[i]],paste0(folder.name,'grass_',season,i,'.tif'),overwrite=T)
    writeRaster(crop[[i]],paste0(folder.name,'crop_',season,i,'.tif'),overwrite=T)
    writeRaster(water[[i]],paste0(folder.name,'water_',season,i,'.tif'),overwrite=T)
    writeRaster(tri[[i]],paste0(folder.name,'TRI_',season,i,'.tif'),overwrite=T)
    writeRaster(evi.sd[[i]],paste0(folder.name,'EVI_',season,'_sd',i,'.tif'),overwrite=T)
    writeRaster(unused[[i]],paste0(folder.name,'unused_',season,i,'.tif'),overwrite=T)
    writeRaster(humanDen[[i]],paste0(folder.name,'humanDen_',season,i,'.tif'),overwrite=T)
    writeRaster(road[[i]],paste0(folder.name,'road_',season,i,'.tif'),overwrite=T)
    writeRaster(pa[[i]],paste0(folder.name,'pa_',season,i,'.tif'),overwrite=T)
    writeRaster(disLitArea[[i]],paste0(folder.name,'litArea_',season,i,'.tif'),overwrite=T)
  cat('saved \n')
}
