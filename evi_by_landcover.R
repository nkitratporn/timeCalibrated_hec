# LOAD LIBRARY ###################
# spaital related
library("raster")
library("rgeos")
library("rgdal")
library("sp")

# modeling
library("dismo")     # contain various sdm analysis tools
library("ENMeval")   # model selection and evaluation tools
library('usdm')      # VIF
# visualization
library('ggplot2')   # graphing/plotting tools
library('reshape')   # data.frameame manipulation tools
library('rasterVis')
library('RColorBrewer')

library('tidyr')
library(dplyr)
library(stringr)
################## SET WORKING DIR ###################
setwd('D:/Academic/Research/0_maximus-materials/analysis/')


################## Read in EVI and LC ###################
studyArea <- readOGR('data/sdm_data/studyarea/studyArea_TH.shp')
r.raster <- raster()
env.ref <- raster('data/sdm_data/data_prep/season_average/KBDI_dry_anomaly.tif')
extent(r.raster) <- extent(env.ref)
res(r.raster) <- 500
r.studyArea <- rasterize(studyArea, r.raster)

#-- read in raster
evi.slope <- stack('data/sdm_data/data_prep/season_by_year/evi_slope_dry2008-wet2018.tif')
evi <- stack('data/sdm_data/data_prep/season_by_year/evi_dry2008-wet2018.tif')
evi.sd <- stack('data/sdm_data/data_prep/season_by_year/eviSD_dry2008-wet2018.tif')
names(evi.slope) <- c('evi_dry1','evi_wet2','evi_dry3','evi_wet4','evi_dry5','evi_wet6','evi_dry7','evi_wet8','evi_dry9','evi_wet10',
                'evi_dry11','evi_wet12','evi_dry13','evi_wet14','evi_dry15','evi_wet16','evi_dry17','evi_wet18','evi_dry19','evi_wet20')
names(evi) <- c('evi_dry1','evi_wet2','evi_dry3','evi_wet4','evi_dry5','evi_wet6','evi_dry7','evi_wet8','evi_dry9','evi_wet10',
                'evi_dry11','evi_wet12','evi_dry13','evi_wet14','evi_dry15','evi_wet16','evi_dry17','evi_wet18','evi_dry19','evi_wet20')
names(evi.sd) <- c('evi_dry1','evi_wet2','evi_dry3','evi_wet4','evi_dry5','evi_wet6','evi_dry7','evi_wet8','evi_dry9','evi_wet10',
                'evi_dry11','evi_wet12','evi_dry13','evi_wet14','evi_dry15','evi_wet16','evi_dry17','evi_wet18','evi_dry19','evi_wet20')
lc <- stack('data/sdm_data/data_prep/season_by_year/lcStack.tif')

########### extract values from related LC pixels ######

lc.mask <- mask(lc,r.studyArea)
lc.mask<- stack(stack(replicate(2,lc.mask[[1]])),stack(replicate(2,lc.mask[[2]])),stack(replicate(2,lc.mask[[3]])),
                stack(replicate(2,lc.mask[[4]])),stack(replicate(2,lc.mask[[5]])),stack(replicate(2,lc.mask[[6]])),
                stack(replicate(2,lc.mask[[7]])),stack(replicate(2,lc.mask[[8]])),stack(replicate(2,lc.mask[[9]])),
                stack(replicate(2,lc.mask[[10]])))


variable_within_landcover <- function(varstack,lcstack,class,name){
  lcstack[lcstack != class] <- NA
  varstack[is.na(lcstack)] <- NA
  var.df <- as.data.frame(varstack, xy=TRUE) %>%
    melt(id.vars=c('x','y'))
  var.df <- na.omit(var.df)
  var.df$class <- class
  var.df$season <- ifelse(grepl('wet',var.df$variable,ignore.case = T),'wet','dry')
  var.df$year <- str_extract(var.df$variable, "\\-*\\d+\\.*\\d*")
  var.df$var <- name
  var.df <- var.df[,!(names(var.df)) %in% 'variable']
  return(var.df)
}

#forest
evi.forest <- variable_within_landcover(evi,lc.mask,1,'evi')
evi.slope.forest <- variable_within_landcover(evi.slope,lc.mask,1,'evi_slope')
evi.sd.forest <- variable_within_landcover(evi.sd,lc.mask,1,'evi_sd')
forest_df <- rbind(evi.forest,evi.slope.forest,evi.sd.forest);
rm(evi.forest,evi.slope.forest,evi.sd.forest)

# shrub
evi.shrub <- variable_within_landcover(evi,lc.mask,2,'evi')
evi.slope.shrub <- variable_within_landcover(evi.slope,lc.mask,2,'evi_slope')
evi.sd.shrub <- variable_within_landcover(evi.sd,lc.mask,2,'evi_sd')
shrub_df <- rbind(evi.shrub,evi.slope.shrub,evi.sd.shrub);
rm(evi.shrub,evi.slope.shrub,evi.sd.shrub)

#woody
evi.woody <- variable_within_landcover(evi,lc.mask,3,'evi')
evi.slope.woody <- variable_within_landcover(evi.slope,lc.mask,3,'evi_slope')
evi.sd.woody <- variable_within_landcover(evi.sd,lc.mask,3,'evi_sd')
woody_df <- rbind(evi.woody,evi.slope.woody,evi.sd.woody);
rm(evi.woody,evi.slope.woody,evi.sd.woody)

# savanna
evi.savanna <- variable_within_landcover(evi,lc.mask,4,'evi')
evi.slope.savanna <- variable_within_landcover(evi.slope,lc.mask,4,'evi_slope')
evi.sd.savanna <- variable_within_landcover(evi.sd,lc.mask,4,'evi_sd')
savanna_df <- rbind(evi.savanna,evi.slope.savanna,evi.sd.savanna);
rm(evi.savanna,evi.slope.savanna,evi.sd.savanna)

# grassland
evi.grassland <- variable_within_landcover(evi,lc.mask,5,'evi')
evi.slope.grassland <- variable_within_landcover(evi.slope,lc.mask,5,'evi_slope')
evi.sd.grassland <- variable_within_landcover(evi.sd,lc.mask,5,'evi_sd')
grassland_df <- rbind(evi.grassland,evi.slope.grassland,evi.sd.grassland);
rm(evi.grassland,evi.slope.grassland,evi.sd.grassland)

# crop
evi.crop <- variable_within_landcover(evi,lc.mask,6,'evi')
evi.slope.crop <- variable_within_landcover(evi.slope,lc.mask,6,'evi_slope')
evi.sd.crop <- variable_within_landcover(evi.sd,lc.mask,6,'evi_sd')
crop_df <- rbind(evi.crop,evi.slope.crop,evi.sd.crop);
rm(evi.crop,evi.slope.crop,evi.sd.crop)

# urban
evi.urban <- variable_within_landcover(evi,lc.mask,7,'evi')
evi.slope.urban <- variable_within_landcover(evi.slope,lc.mask,7,'evi_slope')
evi.sd.urban <- variable_within_landcover(evi.sd,lc.mask,7,'evi_sd')
urban_df <- rbind(evi.urban,evi.slope.urban,evi.sd.urban);
rm(evi.urban,evi.slope.urban,evi.sd.urban)

#-- combine
evi.df <- rbind(forest_df,shrub_df,woody_df,savanna_df,grassland_df,crop_df,urban_df)
rm(forest_df,shrub_df,woody_df,savanna_df,grassland_df,crop_df,urban_df)
write.csv(evi.df,'data/sdm_output/variable_within_lc_reasponses.csv')

################## plot evi characteristic by LC ###################

#-- prepare data
summary(evi.df)
evi.df$year <- as.integer(evi.df$year); evi.df$season <- as.factor(evi.df$season); evi.df$var <- as.factor(evi.df$var);
head(evi.df,10)

evi.df_from2014 <- evi.df[evi.df$year >= 11, ]
summary(evi.df_from2014)

#library(devtools)
#install_github("clauswilke/ggridges")
library(ggridges)


evi.class <- ggplot(evi.df_from2014[evi.df_from2014$var=='evi',],
                  aes(x=value, y=factor(class),fill=paste(class,season))) +
  geom_density_ridges(alpha =0.7, color='black', scale=1.5,from=0,to=1) +
  scale_fill_cyclical(
    breaks = c('1 wet','1 dry'),
    labels = c('1 wet' = 'Wet', '1 dry'='Dry'),
    values = c('#ff8080','#8080ff'),
    name='',guide='legend'
  ) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  theme(legend.position = c(0.75,0.2)) +
  labs(x = 'EVI', y = '')

evi.slope.class <- ggplot(evi.df_from2014[evi.df_from2014$var=='evi_slope',],
                        aes(x=value, y=factor(class),fill=paste(class,season))) +
  geom_density_ridges(alpha =0.7, color='black', scale=1,from=-0.1,to=0.1) +
  scale_fill_cyclical(
    breaks = c('forest wet','forest dry'),
    labels = c('forest wet' = 'Wet', 'forest dry'='Dry'),
    values = c('#ff8080','#8080ff'),
    name='Season',guide='none') +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  labs(x = 'EVI Slope', y = '')

evi.sd.class <- ggplot(evi.df_from2014[evi.df_from2014$var=='evi_sd',],
                        aes(x=value, y=factor(class),fill=paste(class,season))) +
  geom_density_ridges(alpha =0.7, color='black', scale=1.5,from=0,to=0.15) +
  scale_fill_cyclical(
    breaks = c('forest wet','forest dry'),
    labels = c('forest wet' = 'Wet', 'forest dry'='Dry'),
    values = c('#ff8080','#8080ff'),
    name='Season',guide='none') +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  labs(x = 'EVI Standard Deviation', y = '')

library(gridExtra)
grid.arrange(evi.slope.class,evi.sd.class,evi.class, ncol=3)

############################################################
################## area by LC time-seires ###################

#--- load raster for land cover
lc <- stack('data/sdm_data/data_prep/season_by_year/lcStack.tif')
names(lc) <- c('lc2009','lc2010','lc2011','lc2012','lc2013','lc2014','lc2015','lc2016','lc2017','lc2018')
lc <- mask(lc,r.studyArea)
#--- load polygon of provinces within study area
province <- readOGR('data/sdm_data/studyarea/studyArea_provinces.shp')

#--- calc frequency of each factor/class for each raster one-by-one
lc.area <- data.frame()
lc.area <- freq(lc,merge=T)
lc.area$prov <- c('All')
for (i in seq_along(province)) {
  scenario.prov <- mask(lc,province[i,]) # raster by province
  lc.prov <- freq(scenario.prov,merge=T) # get freqency of class
  lc.prov$prov <- province@data[i,5] # add name of province
  lc.area <- rbind(lc.area,lc.prov) # add to data.frame
}
lc.area <- na.omit(lc.area,cols=1)
lc.area.melt <- melt(lc.area,id=c('prov','value')) # melt to basic data by key
lc.area.melt <- as.data.frame(unclass(lc.area.melt))
colnames(lc.area.melt) <- c('province','class','year','count')

#--- extract season, scenario, year from raster name as saved in type
lc.area.melt$area <- lc.area.melt$count*0.5*0.5
lc.area.melt$year <- substring(lc.area.melt$year,3)
lc.area.melt$year <- as.factor(lc.area.melt$year)

#--- save
write.csv(lc.area.melt,'data/sdm_output/occ-timeDependent/prediction/25-results/lc_area_all.csv',row.names = F,col.names = T)

## create graphs to visulize area #####
lc.area <- read.csv('data/sdm_output/occ-timeDependent/prediction/25-results/lc_area_all.csv')
lc.area.overall <- lc.area[(lc.area$province=='All'),];

class.labs <- c("Forest", "Shrub", "Woody savannas",'Savannas','Grasslands','Croplands','Built-up')
names(class.labs) <- c("1", "2", "3", "4",'5','6','7')

lc.area.barAll <- ggplot(lc.area.overall[lc.area.overall$class > 0 & lc.area.overall$class <=7 ,],
                         aes(y=area,x=factor(year), fill = factor(class,levels=c(1,2,3,4,5,6,7)))) +
  geom_bar(position='stack',stat='identity') +
  scale_fill_brewer(palette='Dark2',labels = c("Forest", "Shrub", "Woody savannas",'Savannas','Grasslands','Croplands','Built-up'))+
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = 'bottom', legend.direction = "horizontal",
        legend.title = element_blank(), legend.text = element_text(size = 15), strip.text = element_text(size=16),
        axis.text = element_text(size=14),axis.title = element_text(size=15)) +
  theme(panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())+
  labs(y=bquote('Area'~(km^2)))

lc.area.line <- ggplot(lc.area.overall[lc.area.overall$class > 0 & lc.area.overall$class <=7 ,],
                       aes(y=area,x=factor(year), group = factor(class,levels=c(1,2,3,4,5,6,7)))) +
  geom_line(aes(color=factor(class)),size=1.5) +
  geom_point(size=1.5) +
  scale_color_brewer(palette='Dark2',labels = c("Forest", "Shrub", "Woody savannas",'Savannas','Grasslands','Croplands','Built-up'))+
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = 'bottom', legend.direction = "horizontal",
        legend.title = element_blank(), legend.text = element_text(size = 14), strip.text = element_text(size=16),
        axis.text = element_text(size=12),axis.title = element_text(size=14)) +
  theme(strip.background =element_rect(fill='#EDEDED', color='#EDEDED'),
        #panel.grid.major = element_line(color = "#EDEDED"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())+
  guides(colour = guide_legend(reverse=T))+
  labs(y=bquote('Area'~(km^2)))

