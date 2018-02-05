####################################################################################
#######    object: INTERSECT LOGGING ROADS WITH AA POINTS       ####################
#######    Update : 2017/11/27                                  ####################
#######    contact: remi.dannunzio@fao.org                      ####################
####################################################################################
library(foreign)
library(plyr)
library(rgeos)
library(rgdal)
library(raster)
library(ggplot2)

options(stringsAsFactors = F)

#################### Define AOI as country boundaries 
workdir    <- "/media/dannunzio/OSDisk/Users/dannunzio/Documents/countries/guinea_ecuatorial/"

setwd(workdir)

################################################################################
#################### BIOKO DATASET
################################################################################

#################### READ POINT FILE AND MAP
pts <- read.csv("estudio_dd/analisis_all_points_20171127.csv")

#################### SPATIALIZE POINTS
pt_df_geo <- SpatialPointsDataFrame(
  coords = pts[,c("location_x","location_y")],
  data   = data.frame(pts),
  proj4string=CRS("+init=epsg:4326")
)

plot(pt_df_geo)

#################### READ ROAD MAP AND BUFFER 50m
road <- readOGR("gis_data/GNQ_loggingroads/GNQ_loggingroads.shp")

plot(road)
proj4string(road) <- proj4string(pt_df_geo)

road_buff      <- gBuffer(road,byid=T,width=0.002)
road_buff@data <- road_buff@data[,c("oid","highway")]

#writeOGR(road_buff,"gis_data/GNQ_loggingroads/buffer_200m.shp","buffer_200m","ESRI Shapefile")

#################### INTERSECT WITH POINTS
pt_df_geo$type <- over(x = pt_df_geo,y=road_buff)[,"highway"]
df             <- pt_df_geo@data

ch <- df[df$causa_dd == "infraestructura",]

table(ch$type,useNA = "always")

write.csv(df,"estudio_dd/points_logging_road.csv",row.names = F)
