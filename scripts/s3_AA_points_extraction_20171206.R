####################################################################################
#######    object: ACCURACY ASSESSMENT FROM OLD POINTS NEW MAP  ####################
#######    Update : 2017/11/08                                  ####################
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
workdir    <- "/media/dannunzio/OSDisk/Users/dannunzio/Documents/countries/guinea_ecuatorial/aa_gnq_forest_change/aa_new_map_20171108/"

setwd(workdir)

################################################################################
#################### BIOKO DATASET
################################################################################

#################### READ POINT FILE AND MAP
rasname <- "../../gis_data/degradation_option2_por_agente_poligonos/uni_version_20171108/uni_map_dd_bioko_aea_20171206.tif"
ptsname          <- "input_gnq_aa_20171108/BIOKO_hotspot_excel.csv"
aa_rasname       <- "results_gnq_20171206/bioko_dd_20171206.tif"
aa_point_name    <- "results_gnq_20171206/pts_bioko_aa_ddmap_20171206.csv"
area_point_name  <- "results_gnq_20171206/area_bioko_20171206.csv"

pts <- read.csv(ptsname,sep=";")
ras <- raster(rasname)

#################### HARMONIZE COLUMN NAMES IN POINT FILE
pts$x <- as.numeric(gsub(pts$x,pattern = ",",replacement = "."))
pts$y <- as.numeric(gsub(pts$y,pattern = ",",replacement = "."))


#################### SPATIALIZE POINTS
pt_df_geo <- SpatialPointsDataFrame(
  coords = pts[,c("x","y")],
  data   = data.frame(pts),
  proj4string=CRS("+init=epsg:4326")
)

pt_df_aea <- spTransform(pt_df_geo,proj4string(ras))
plot(pt_df_aea)
#################### INTERSECT POINTS WITH MAP

pt_df_aea@data$DD_code <- extract(ras,pt_df_aea)

points <- pt_df_aea@data
points$unique_id <- row(points)[,1]
names(points)
summary(points)
head(points)
table(points$DD_code,points$ref_class)

out <- pt_df_aea[points$DD_code ==0,]
plot(out,add=T,col="red")

#################### GENERATE TO COLUMN TO MATCH MAP AND REFERENCE DATA
unicode <- read.table("results_gnq_20171206/reclass_uni_aa_codes.txt")
names(unicode) <- c("uni_code","aa_code")

df <- merge(points,unicode,by.x="DD_code",by.y="uni_code",all.x=T)
df <- arrange(df,unique_id)
table(df$aa_code,df$DD_code,useNA="always")
table(df$aa_code,df$ref_class,useNA="always")

reclass <- reclassify(ras,unicode)


#################### DOUBLE CHECK THAT CHANGING COLUMN LABELS AND CHANGING MAP GIVE SAME RESULTS
df$verif <- extract(reclass,pt_df_aea)
table(df$verif,df$aa_code,useNA="always")

df1 <- df[!is.na(df$aa_code),]

#################### EXPORT
write.csv(df1,aa_point_name,row.names = F)
writeRaster(reclass,aa_rasname,overwrite=T)

#### COMPUTE MAP AREAS
system(sprintf("oft-stat -i %s -o %s -um %s -nostd",
               aa_rasname,
               "results_gnq_20171206/tmp_stats_dd_20171206.txt",
               aa_rasname
))

areas <- read.table("results_gnq_20171206/tmp_stats_dd_20171206.txt")[,1:2]
names(areas) <- c("map_code","map_area")
pix <- res(raster(aa_rasname))[1]
areas$map_area <- areas$map_area*pix*pix/10000
write.csv(areas,area_point_name,row.names = F)

################################################################################
#################### CONTINENTE DATASET
################################################################################

#################### READ POINT FILE AND MAP
rasname <- "../../gis_data/degradation_option2_por_agente_poligonos/uni_version_20171108/uni_map_dd_continente_aea_20171206.tif"
ptsname          <- "input_gnq_aa_20171108/continente_hotspot_excel.csv"
aa_rasname       <- "results_gnq_20171206/continente_dd_20171206.tif"
aa_point_name    <- "results_gnq_20171206/pts_continente_aa_ddmap_20171206.csv"
area_point_name  <- "results_gnq_20171206/area_continente_20171206.csv"

pts <- read.csv(ptsname,sep=";")
ras <- raster(rasname)

#################### HARMONIZE COLUMN NAMES IN POINT FILE
pts$x <- as.numeric(gsub(pts$x,pattern = ",",replacement = "."))
pts$y <- as.numeric(gsub(pts$y,pattern = ",",replacement = "."))


#################### SPATIALIZE POINTS
pt_df_geo <- SpatialPointsDataFrame(
  coords = pts[,c("x","y")],
  data   = data.frame(pts),
  proj4string=CRS("+init=epsg:4326")
)

pt_df_aea <- spTransform(pt_df_geo,proj4string(ras))
plot(pt_df_aea)
#################### INTERSECT POINTS WITH MAP

pt_df_aea@data$DD_code <- extract(ras,pt_df_aea)

points <- pt_df_aea@data
points$unique_id <- row(points)[,1]
names(points)
summary(points)
head(points)
table(points$DD_code,points$ref_class)

out <- pt_df_aea[points$DD_code ==0,]
plot(out,add=T,col="red")

#################### GENERATE TO COLUMN TO MATCH MAP AND REFERENCE DATA
unicode <- read.table("results_gnq_20171206/reclass_uni_aa_codes.txt")
names(unicode) <- c("uni_code","aa_code")

df <- merge(points,unicode,by.x="DD_code",by.y="uni_code",all.x=T)
df <- arrange(df,unique_id)
table(df$aa_code,df$DD_code,useNA="always")
table(df$aa_code,df$ref_class,useNA="always")

reclass <- reclassify(ras,unicode)


#################### DOUBLE CHECK THAT CHANGING COLUMN LABELS AND CHANGING MAP GIVE SAME RESULTS
df$verif <- extract(reclass,pt_df_aea)
table(df$verif,df$aa_code,useNA="always")

df1 <- df[!is.na(df$aa_code),]

#################### EXPORT
write.csv(df1,aa_point_name,row.names = F)
writeRaster(reclass,aa_rasname,overwrite=T)

#### COMPUTE MAP AREAS
system(sprintf("oft-stat -i %s -o %s -um %s -nostd",
               aa_rasname,
               "results_gnq_20171206/tmp_stats_dd_20171206.txt",
               aa_rasname
))

areas <- read.table("results_gnq_20171206/tmp_stats_dd_20171206.txt")[,1:2]
names(areas) <- c("map_code","map_area")
pix <- res(raster(aa_rasname))[1]
areas$map_area <- areas$map_area*pix*pix/10000
write.csv(areas,area_point_name,row.names = F)


################################################################################
#################### ANNOBON DATASET
################################################################################

#################### READ POINT FILE AND MAP
rasname <- "../../gis_data/degradation_option2_por_agente_poligonos/uni_version_20171108/uni_map_dd_annobon_aea_20171206.tif"
aa_rasname       <- "results_gnq_20171206/annobon_dd_20171206.tif"
area_point_name  <- "results_gnq_20171206/area_annobon_20171206.csv"

ras <- raster(rasname)


reclass <- reclassify(ras,unicode)

#################### EXPORT
writeRaster(reclass,aa_rasname,overwrite=T)

#### COMPUTE MAP AREAS
system(sprintf("oft-stat -i %s -o %s -um %s -nostd",
               aa_rasname,
               "results_gnq_20171206/tmp_stats_dd_20171206.txt",
               aa_rasname
))

areas <- read.table("results_gnq_20171206/tmp_stats_dd_20171206.txt")[,1:2]
names(areas) <- c("map_code","map_area")
pix <- res(raster(aa_rasname))[1]
areas$map_area <- areas$map_area*pix*pix/10000
write.csv(areas,area_point_name,row.names = F)
