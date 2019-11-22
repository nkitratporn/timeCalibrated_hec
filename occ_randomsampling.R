# library
# spaital related
library("raster")
library("rgeos")
library("rgdal")
library("sp")

#########################
# set working folder if necessary
setwd('D:/Academic/Research/0_maximus-materials/analysis/')

dir.create('data/sdm_data/studyarea/occ_time-dependent/all_data/random-sampling')
folder <- paste0('data/sdm_data/studyarea/occ_time-dependent/all_data/random-sampling/')
set.seed(123)
for(k in 2014:2018){
  villageBuffer <- readOGR(paste0('data/sdm_data/studyarea/occ_time-dependent/all_data/village-occ-all-',k,'.shp'))
  village.id <- villageBuffer$NRDCODE
  village.occWet <- as.integer(villageBuffer$wet)
  village.occDry <- as.integer(villageBuffer$dry)

  #adjusting missing value
  if(k==2017){village.occDry[78] <- 4} #2017
  if(k==2018){village.occWet[78] <- 4; village.occWet[73] <- 5;} #2018
  set.seed(30)
  for(j in 1:5){
    cat('strat:',j,' ...')
    occ.wet <- SpatialPoints(data.frame(x = 0, y = 0))[-1,]
    crs(occ.wet) <- crs(villageBuffer)
    for (i in 1:length(village.id)) {
      cat(i,'-')
      if(is.na(village.occWet[i])) next
      temp.sp <- spsample(villageBuffer[i,],n=village.occWet[i],'random',iter=50)
      temp.sp$vill_id <- rep(village.id[i],village.occWet[i])
      temp.sp$season <- rep('occ_wet',village.occWet[i])
      occ.wet <- union(occ.wet,temp.sp)
      cat(i,' ')
    }
    cat('end occ.wet ...')
    occ.dry <- SpatialPoints(data.frame(x = 0, y = 0))[-1,]
    crs(occ.dry) <- crs(villageBuffer)
    for (i in 1:length(village.id)) {
      if(is.na(village.occDry[i])) next
      temp.sp <- spsample(villageBuffer[i,],n=village.occDry[i],'random',iter=50)
      temp.sp$vill_id <- rep(village.id[i],village.occDry[i])
      temp.sp$season <- rep('occ_dry',village.occDry[i])
      occ.dry <- union(occ.dry,temp.sp)
    }
    cat('end occ.dry ...')
    
    name.occWet <- paste0(folder,'occ_wet_',k,'_',j,'.shp')
    name.occDry <- paste0(folder,'occ_dry_',k,'_',j,'.shp')

    #save
    shapefile(occ.wet,name.occWet)
    shapefile(occ.dry,name.occDry)

    cat('saved...\n')
  }
  
}


#########################
#check value
#sum(village.occDry,na.rm=T); sum(village.occWet,na.rm=T);
#df <- data.frame(no = 1:length(village.id))
#df$id <- village.id;
#df$occ_dry <- village.occDry; df$occ_wet <- village.occWet; #df$val_dry <- village.valDry; df$val_wet <- village.valWet;
#View(df)