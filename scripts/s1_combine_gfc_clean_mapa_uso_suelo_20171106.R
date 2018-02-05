####################################################################################
####### Object:  Processing chain                
####### Author:  remi.dannunzio@fao.org                               
####### Update:  2017/06/05                                      
####################################################################################
library(foreign)
library(plyr)
library(rgeos)
library(rgdal)
library(raster)
library(ggplot2)

options(stringsAsFactors = F)


rootdir <- "/media/dannunzio/OSDisk/Users/dannunzio/Documents/countries/guinea_ecuatorial/gis_data/degradation_option2_por_agente_poligonos/"
setwd(rootdir)

##########################################################################################
#### Clip continente LANDSAT MOSAIC from GFC product
system(sprintf("gdal_translate -projwin 8.86708860759 2.52848101266 12.1329113924 0.610759493671 -co COMPRESS=LZW %s %s",
               paste0(rootdir,"../GFC_hansen_LSAT_mosaics_2014/Hansen_GFC-2015-v1.3_last_10N_010E.tif"),
               paste0(rootdir,"tmp_10_10_cont.tif")
))

system(sprintf("gdal_translate -projwin 8.86708860759 2.52848101266 12.1329113924 0.610759493671 -co COMPRESS=LZW %s %s",
               paste0(rootdir,"../GFC_hansen_LSAT_mosaics_2014/Hansen_GFC-2015-v1.3_last_10N_000E.tif"),
               paste0(rootdir,"tmp_10_00_cont.tif")
))

##########################################################################################
#### Merge bands together
system(sprintf("gdal_merge.py -v -n 0 -co COMPRESS=LZW -o %s %s",
               paste0(rootdir,"tmp_merge.tif"),
               paste0(rootdir,"tmp_10_*.tif")
))

#### Compress the resulting raster
system(sprintf("gdal_translate  -co COMPRESS=LZW %s %s",
               paste0(rootdir,"tmp_merge.tif"),
               paste0(rootdir,"lsat_mosaic_continente_gfc.tif")
))


##########################################################################################
#### Clip bioko
system(sprintf("gdal_translate -projwin 8.38502175633 3.80346509098 8.97886036392 3.19708059731 -co COMPRESS=LZW %s %s",
               paste0(rootdir,"../GFC_hansen_LSAT_mosaics_2014/Hansen_GFC-2015-v1.3_last_10N_000E.tif"),
               paste0(rootdir,"lsat_mosaic_bioko_gfc.tif")
))

##########################################################################################
##########################################################################################
##########################################################################################
#### BIOKO

##########################################################################################
#### Clip classification to correct extent
system(sprintf("gdal_translate -a_nodata none -projwin 8.38502175633 3.80346509098 8.97886036392 3.19708059731 -co COMPRESS=LZW %s %s",
               paste0(rootdir,"../Clasificaciones_originales_LatLong/b2_5c_x.tif"),
               paste0(rootdir,"bioko_classificacion.tif")
))

##########################################################################################
#### Set no_data to zero for classification
system(sprintf("gdal_calc.py -A %s --outfile=%s --calc=\"%s\"",
               paste0(rootdir,"bioko_classificacion.tif"),
               paste0(rootdir,"bioko_classificacion_nd.tif"),
               "(A==255)*0+(A<255)*A"
))

#### Clip GFC clean to correct extent
system(sprintf("gdal_translate -a_nodata none -co COMPRESS=LZW %s %s",
               paste0(rootdir,"../gfc_cleaning/cleaning_20170612/biokogfc2_ndd.tif"),
               paste0(rootdir,"bioko_gfc_clean.tif")
))

##########################################################################################
#### Set no_data to zero for GFC clean
system(sprintf("gdal_calc.py -A %s --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               paste0(rootdir,"bioko_gfc_clean.tif"),
               paste0(rootdir,"bioko_gfc_clean_nd.tif"),
               "(A==65535)*0+(A<65535)*A"
               ))

##########################################################################################
#### Clump results
system(sprintf("oft-clump -i %s -o %s -um %s",
               paste0(rootdir,"bioko_classificacion_nd.tif"),
               paste0(rootdir,"bioko_classificacion_clump.tif"),
               paste0(rootdir,"bioko_classificacion_nd.tif")
))

##########################################################################################
#### Compute clump stats
system(sprintf("oft-stat -i %s -o %s -um %s",
               paste0(rootdir,"bioko_classificacion_clump.tif"),
               paste0(rootdir,"stat_bioko_classificacion_clump.txt"),
               paste0(rootdir,"bioko_classificacion_clump.tif")
))

#########################################################################################
### Read and create a reclass code for clumps > 50 pixels
df <- read.table(paste0(rootdir,"stat_bioko_classificacion_clump.txt"))[,1:2]

names(df) <- c("clump_id","size")
summary(df)

df$new <- 1
df[df$size > 50,]$new <- 2

table(df$new)

write.table(df,paste0(rootdir,"bioko_clump_gt_50.txt"),row.names = F,col.names = F)

####### Reclassificar para la masquera
system(sprintf("(echo %s; echo 1; echo 1; echo 3; echo 0) | oft-reclass  -oi %s  -um %s %s",
               paste0(rootdir,"bioko_clump_gt_50.txt"),
               paste0(rootdir,"bioko_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"bioko_classificacion_clump.tif"),
               paste0(rootdir,"bioko_classificacion_clump.tif")

))

####### Segmentacion del imagen
system(sprintf("(echo 0; echo 0 ; echo 0)|oft-seg -region -ttest -automax %s %s",
               paste0(rootdir,"lsat_mosaic_bioko_gfc.tif"),
               paste0(rootdir,"seg_bioko_gfc.tif")
)
)

####### Reclip segments
system(sprintf("oft-clip.pl %s %s %s",
               paste0(rootdir,"bioko_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"seg_bioko_gfc.tif"),
               paste0(rootdir,"seg_bioko_gfc_clip.tif")
)
)

####### Polygonize
system(sprintf("gdal_polygonize.py -f \"ESRI Shapefile\" %s %s",
               paste0(rootdir,"seg_bioko_gfc_clip.tif"),
               paste0(rootdir,"seg_bioko_gfc_clip.shp")
)
)

####### Reclip GFC_perdidas cleaned
system(sprintf("oft-clip.pl %s %s %s",
               paste0(rootdir,"bioko_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"bioko_gfc_clean_nd.tif"),
               paste0(rootdir,"bioko_gfc_clean_nd_clip.tif")
)
)

####### Apply mask to identify segments from both products
system(sprintf("gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc=\"%s\"",
               paste0(rootdir,"seg_bioko_gfc_clip.tif"),
               paste0(rootdir,"bioko_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"bioko_classificacion_clump.tif"),
               paste0(rootdir,"tmp_bioko_clump_mixed.tif"),
               "(B>0)*((B==1)*C+(B==2)*A)"
))

#################### SIEVE RESULTS x2
system(sprintf("gdal_sieve.py -st %s %s %s",
               2,
               paste0(rootdir,"tmp_bioko_clump_mixed.tif"),
               paste0(rootdir,"tmp_bioko_clump_mixed_sieve02.tif")
))

#################### SIEVE RESULTS x4
system(sprintf("gdal_sieve.py -st %s %s %s",
               4,
               paste0(rootdir,"tmp_bioko_clump_mixed_sieve02.tif"),
               paste0(rootdir,"tmp_bioko_clump_mixed_sieve04.tif")
))

#################### SIEVE RESULTS x8
system(sprintf("gdal_sieve.py -st %s %s %s",
               8,
               paste0(rootdir,"tmp_bioko_clump_mixed_sieve04.tif"),
               paste0(rootdir,"tmp_bioko_clump_mixed_sieve08.tif")
))


#################### SIEVE RESULTS x12
system(sprintf("gdal_sieve.py -st %s %s %s",
               12,
               paste0(rootdir,"tmp_bioko_clump_mixed_sieve08.tif"),
               paste0(rootdir,"tmp_bioko_clump_mixed_sieve12.tif")
))

####### Reclump results
system(sprintf("oft-clump -i %s -o %s -um %s",
               paste0(rootdir,"tmp_bioko_clump_mixed_sieve12.tif"),
               paste0(rootdir,"bioko_clump_mix.tif"),
               paste0(rootdir,"tmp_bioko_clump_mixed_sieve12.tif")
))

####### Compute zonal stats for LOSSES
system(sprintf("oft-zonal_large_list.py -i %s -o %s -um %s",
               paste0(rootdir,"bioko_gfc_clean_nd_clip.tif"),
               paste0(rootdir,"zonal_perdidas_bioko.txt"),
               paste0(rootdir,"bioko_clump_mix.tif")
))

####### Compute zonal stats for CLASSIFICATION
system(sprintf("oft-zonal_large_list.py -i %s -o %s -um %s",
               paste0(rootdir,"bioko_classificacion_nd.tif"),
               paste0(rootdir,"zonal_classif_bioko.txt"),
               paste0(rootdir,"bioko_clump_mix.tif")
))


####### Read zonal stats and rename columns
df_perdidas <- read.table(paste0(rootdir,"zonal_perdidas_bioko.txt"))
df_classif  <- read.table(paste0(rootdir,"zonal_classif_bioko.txt"))

names(df_classif)  <- c("clump_id","total","nodata",paste0("class_",1:5))

codes_gfclimpio <- c(1,2,3,4,5,6,7,8,9,10,105,106,107,108,109,110,111,112,113,114,190,195,205,206,207,208,209,210,211,212,213,214,290,295,305,306,307,308,309,310,311,312,313,314,390,395,405,406,407,408,409,410,411,412,413,414,490,495,505,506,507,508,509,510,511,512,513,514,590)
df_perdida2 <- df_perdidas[,c(1,2,3,3+codes_gfclimpio)]
names(df_perdida2) <- c("clump_id","total","nodata",paste0("perd_",codes_gfclimpio))

codes_bosque  <- paste0("perd_",c(1,2,3,4,5,190,195,290,295,390,395,490,495,590))
codes_nobosq  <- paste0("perd_",c(6,7,8,9))
codes_agua    <- paste0("perd_",c(10))
codes_perd_bs <- paste0("perd_",c(105,106,107,108,109,110,111,112,113,114,305,306,307,308,309,310,311,312,313,314,405,406,407,408,409,410,411,412,413,414,505,506,507,508,509,510,511,512,513,514))
codes_perd_mg <- paste0("perd_",c(205,206,207,208,209,210,211,212,213,214))

dbf1 <- cbind(df_classif[,c(1,4:8)],df_perdida2)
head(dbf1)


dbf1$new_class <- max.col(dbf1[,paste0("class_",1:5)])

dbf1$new_perd  <- 1

dbf1[rowSums(dbf1[,codes_bosque]) <= 0.3*dbf1$total & rowSums(dbf1[,codes_perd_bs]) > 0,]$new_perd <- 2
dbf1[rowSums(dbf1[,codes_bosque]) >  0.3*dbf1$total & rowSums(dbf1[,codes_perd_bs]) > 0,]$new_perd <- 3

dbf1[rowSums(dbf1[,codes_bosque]) <= 0.3*dbf1$total & rowSums(dbf1[,codes_perd_mg]) > 0,]$new_perd <- 4
dbf1[rowSums(dbf1[,codes_bosque]) >  0.3*dbf1$total & rowSums(dbf1[,codes_perd_mg]) > 0,]$new_perd <- 5

dbf1[rowSums(dbf1[,codes_nobosq]) > 0.9*dbf1$total,]$new_perd <- 6
dbf1[dbf1[,codes_agua]            > 0.9*dbf1$total,]$new_perd <- 7

dbf1$combi <- dbf1$new_perd*10+dbf1$new_class

table(dbf1$combi)
sum(dbf1[dbf1$total < 12,"total"]) / sum(dbf1$total)*100

write.table(dbf1[,c("clump_id","total","combi")],paste0(rootdir,"bioko_reclass.txt"),row.names = F,col.names = F)

####### Reclassificar para la masquera
system(sprintf("(echo %s; echo 1; echo 1; echo 3; echo 0) | oft-reclass  -oi %s  -um %s %s",
               paste0(rootdir,"bioko_reclass.txt"),
               paste0(rootdir,"tmp_bioko_clump_mix_reclass.tif"),
               paste0(rootdir,"bioko_clump_mix.tif"),
               paste0(rootdir,"bioko_clump_mix.tif")
               
))

####### Reproject in Albers
system(sprintf("gdalwarp -t_srs %s -ot Byte -co COMPRESS=LZW %s %s",
               "../albers_equal_area_projection.txt",
               paste0(rootdir,"tmp_bioko_clump_mix_reclass.tif"),
               paste0(rootdir,"aea_bioko_map_dd_20171106.tif")
)
)

####### Comprimir el resultado
system(sprintf("gdal_translate -ot Byte -co COMPRESS=LZW %s %s",
               paste0(rootdir,"tmp_bioko_clump_mix_reclass.tif"),
               paste0(rootdir,"geo_bioko_map_dd_20171106.tif")
)
)

##########################################################################################
##########################################################################################
##########################################################################################
#### CONTINENTE

##########################################################################################
#### Clip classification to correct extent
system(sprintf("gdal_translate -a_nodata none -projwin 8.86708860759 2.52848101266 12.1329113924 0.610759493671 -co COMPRESS=LZW %s %s",
               paste0(rootdir,"../Clasificaciones_originales_LatLong/c2_5c_20170607.tif"),
               paste0(rootdir,"continente_classificacion.tif")
))

##########################################################################################
#### Set no_data to zero for classification
system(sprintf("gdal_calc.py -A %s --outfile=%s --calc=\"%s\"",
               paste0(rootdir,"continente_classificacion.tif"),
               paste0(rootdir,"continente_classificacion_nd.tif"),
               "(A==255)*0+(A<255)*A"
))

#### Clip GFC clean to correct extent
system(sprintf("gdal_translate -a_nodata none -co COMPRESS=LZW %s %s",
               paste0(rootdir,"../gfc_cleaning/cleaning_20170612/cont2gfc_nd_cmp.tif"),
               paste0(rootdir,"continente_gfc_clean.tif")
))

##########################################################################################
#### Set no_data to zero for GFC clean
system(sprintf("gdal_calc.py -A %s --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               paste0(rootdir,"continente_gfc_clean.tif"),
               paste0(rootdir,"continente_gfc_clean_nd.tif"),
               "(A==65535)*0+(A<65535)*A"
))

##########################################################################################
#### Clump results
system(sprintf("oft-clump -i %s -o %s -um %s",
               paste0(rootdir,"continente_classificacion_nd.tif"),
               paste0(rootdir,"continente_classificacion_clump.tif"),
               paste0(rootdir,"continente_classificacion_nd.tif")
))

##########################################################################################
#### Compute clump stats
system(sprintf("oft-stat -i %s -o %s -um %s",
               paste0(rootdir,"continente_classificacion_clump.tif"),
               paste0(rootdir,"stat_continente_classificacion_clump.txt"),
               paste0(rootdir,"continente_classificacion_clump.tif")
))

##########################################################################################
#### Read and create a reclass code for clumps > 50 pixels
df <- read.table(paste0(rootdir,"stat_continente_classificacion_clump.txt"))[,1:2]

names(df) <- c("clump_id","size")
summary(df)

df$new <- 1
df[df$size > 50,]$new <- 2

table(df$new)

write.table(df,paste0(rootdir,"continente_clump_gt_50.txt"),row.names = F,col.names = F)

####### Reclassificar para la masquera
system(sprintf("(echo %s; echo 1; echo 1; echo 3; echo 0) | oft-reclass  -oi %s  -um %s %s",
               paste0(rootdir,"continente_clump_gt_50.txt"),
               paste0(rootdir,"continente_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"continente_classificacion_clump.tif"),
               paste0(rootdir,"continente_classificacion_clump.tif")

))

####### Segmentacion del imagen
system(sprintf("(echo 0; echo 0 ; echo 0)|oft-seg -region -ttest -automax %s %s",
               paste0(rootdir,"lsat_mosaic_continente_gfc.tif"),
               paste0(rootdir,"seg_continente_gfc.tif")
)
)

####### Reclip segments
system(sprintf("oft-clip.pl %s %s %s",
               paste0(rootdir,"continente_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"seg_continente_gfc.tif"),
               paste0(rootdir,"seg_continente_gfc_clip.tif")
)
)

####### Reclip GFC_perdidas cleaned
system(sprintf("oft-clip.pl %s %s %s",
               paste0(rootdir,"continente_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"continente_gfc_clean_nd.tif"),
               paste0(rootdir,"continente_gfc_clean_nd_clip.tif")
)
)

####### Apply mask to identify segments from both products
system(sprintf("gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc=\"%s\"",
               paste0(rootdir,"seg_continente_gfc_clip.tif"),
               paste0(rootdir,"continente_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"continente_classificacion_clump.tif"),
               paste0(rootdir,"tmp_continente_clump_mixed.tif"),
               "(B>0)*((B==1)*C+(B==2)*A)"
))

#################### SIEVE RESULTS x2
system(sprintf("gdal_sieve.py -st %s %s %s",
               2,
               paste0(rootdir,"tmp_continente_clump_mixed.tif"),
               paste0(rootdir,"tmp_continente_clump_mixed_sieve02.tif")
))

#################### SIEVE RESULTS x4
system(sprintf("gdal_sieve.py -st %s %s %s",
               4,
               paste0(rootdir,"tmp_continente_clump_mixed_sieve02.tif"),
               paste0(rootdir,"tmp_continente_clump_mixed_sieve04.tif")
))

#################### SIEVE RESULTS x8
system(sprintf("gdal_sieve.py -st %s %s %s",
               8,
               paste0(rootdir,"tmp_continente_clump_mixed_sieve04.tif"),
               paste0(rootdir,"tmp_continente_clump_mixed_sieve08.tif")
))


#################### SIEVE RESULTS x12
system(sprintf("gdal_sieve.py -st %s %s %s",
               12,
               paste0(rootdir,"tmp_continente_clump_mixed_sieve08.tif"),
               paste0(rootdir,"tmp_continente_clump_mixed_sieve12.tif")
))

####### Reclump results
system(sprintf("oft-clump -i %s -o %s -um %s",
               paste0(rootdir,"tmp_continente_clump_mixed_sieve12.tif"),
               paste0(rootdir,"continente_clump_mix.tif"),
               paste0(rootdir,"tmp_continente_clump_mixed_sieve12.tif")
))

####### Compute zonal stats for LOSSES
system(sprintf("oft-zonal_large_list.py -i %s -o %s -um %s",
               paste0(rootdir,"continente_gfc_clean_nd_clip.tif"),
               paste0(rootdir,"zonal_perdidas_continente.txt"),
               paste0(rootdir,"continente_clump_mix.tif")
))

####### Compute zonal stats for CLASSIFICATION
system(sprintf("oft-zonal_large_list.py -i %s -o %s -um %s",
               paste0(rootdir,"continente_classificacion_nd.tif"),
               paste0(rootdir,"zonal_classif_continente.txt"),
               paste0(rootdir,"continente_clump_mix.tif")
))

####### Read zonal stats and rename columns
df_perdidas <- read.table(paste0(rootdir,"zonal_perdidas_continente.txt"))
df_classif  <- read.table(paste0(rootdir,"zonal_classif_continente.txt"))

names(df_classif)  <- c("clump_id","total","nodata",paste0("class_",1:5))

codes_gfclimpio <- c(1,2,3,4,5,6,7,8,9,10,105,106,107,108,109,110,111,112,113,114,190,195,205,206,207,208,209,210,211,212,213,214,290,295,305,306,307,308,309,310,311,312,313,314,390,395)
df_perdida2 <- df_perdidas[,c(1,2,3,3+codes_gfclimpio)]
names(df_perdida2) <- c("clump_id","total","nodata",paste0("perd_",codes_gfclimpio))

codes_bosque  <- paste0("perd_",c(1,2,3,4,5,190,195,290,295,390,395))
codes_nobosq  <- paste0("perd_",c(6,7,8,9))
codes_agua    <- paste0("perd_",c(10))
codes_perd_bs <- paste0("perd_",c(105,106,107,108,109,110,111,112,113,114,305,306,307,308,309,310,311,312,313,314))
codes_perd_mg <- paste0("perd_",c(205,206,207,208,209,210,211,212,213,214))

dbf1 <- cbind(df_classif[,c(1,4:8)],df_perdida2)
head(dbf1)

dbf1$new_class <- max.col(dbf1[,paste0("class_",1:5)])

dbf1$new_perd  <- 1

dbf1[rowSums(dbf1[,codes_bosque]) <= 0.3*dbf1$total & rowSums(dbf1[,codes_perd_bs]) > 0,]$new_perd <- 2
dbf1[rowSums(dbf1[,codes_bosque]) >  0.3*dbf1$total & rowSums(dbf1[,codes_perd_bs]) > 0,]$new_perd <- 3

dbf1[rowSums(dbf1[,codes_bosque]) <= 0.3*dbf1$total & rowSums(dbf1[,codes_perd_mg]) > 0,]$new_perd <- 4
dbf1[rowSums(dbf1[,codes_bosque]) >  0.3*dbf1$total & rowSums(dbf1[,codes_perd_mg]) > 0,]$new_perd <- 5

dbf1[rowSums(dbf1[,codes_nobosq]) > 0.9*dbf1$total,]$new_perd <- 6
dbf1[dbf1[,codes_agua]            > 0.9*dbf1$total,]$new_perd <- 7

dbf1$combi <- dbf1$new_perd*10+dbf1$new_class

table(dbf1$combi)

write.table(dbf1[,c("clump_id","total","combi")],paste0(rootdir,"continente_reclass.txt"),row.names = F,col.names = F)

####### Reclassificar para la masquera
system(sprintf("(echo %s; echo 1; echo 1; echo 3; echo 0) | oft-reclass  -oi %s  -um %s %s",
               paste0(rootdir,"continente_reclass.txt"),
               paste0(rootdir,"tmp_continente_clump_mix_reclass.tif"),
               paste0(rootdir,"continente_clump_mix.tif"),
               paste0(rootdir,"continente_clump_mix.tif")
               
))

####### Reproject in Albers
system(sprintf("gdalwarp -t_srs %s -ot Byte -co COMPRESS=LZW %s %s",
               "../albers_equal_area_projection.txt",
               paste0(rootdir,"tmp_continente_clump_mix_reclass.tif"),
               paste0(rootdir,"aea_continente_map_dd_20171106.tif")
)
)

####### Comprimir el resultado
system(sprintf("gdal_translate -ot Byte -co COMPRESS=LZW %s %s",
               paste0(rootdir,"tmp_continente_clump_mix_reclass.tif"),
               paste0(rootdir,"geo_continente_map_dd_20171106.tif")
)
)


##########################################################################################
##########################################################################################
##########################################################################################
#### ANNOBON

##########################################################################################
#### Set no_data to zero for classification
system(sprintf("gdal_translate -a_nodata none -co COMPRESS=LZW %s %s",
               "../Clasificaciones_originales_LatLong/a2_5c_x.tif",
               paste0(rootdir,"tmp_annobon_classificacion_nd.tif")
))

system(sprintf("gdal_calc.py -A %s --outfile=%s --calc=\"%s\"",
               paste0(rootdir,"tmp_annobon_classificacion_nd.tif"),
               paste0(rootdir,"annobon_classificacion_nd.tif"),
               "(A==255)*0+(A<255)*((A<6)*A+(A==6)*7+(A==7)*6)"
))

#### Clip GFC clean to correct extent
system(sprintf("gdal_translate -a_nodata none -co COMPRESS=LZW %s %s",
               paste0(rootdir,"../gfc_cleaning/cleaning_20170612/annobongfcr.tif"),
               paste0(rootdir,"annobon_gfc_clean.tif")
))

##########################################################################################
#### Set no_data to zero for GFC clean
system(sprintf("gdal_calc.py -A %s --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               paste0(rootdir,"annobon_gfc_clean.tif"),
               paste0(rootdir,"annobon_gfc_clean_nd.tif"),
               "(A==65535)*0+(A<65535)*A"
               ))

##########################################################################################
#### Clump results
system(sprintf("oft-clump -i %s -o %s -um %s",
               paste0(rootdir,"annobon_classificacion_nd.tif"),
               paste0(rootdir,"annobon_classificacion_clump.tif"),
               paste0(rootdir,"annobon_classificacion_nd.tif")
))

##########################################################################################
#### Compute clump stats
system(sprintf("oft-stat -i %s -o %s -um %s",
               paste0(rootdir,"annobon_classificacion_clump.tif"),
               paste0(rootdir,"stat_annobon_classificacion_clump.txt"),
               paste0(rootdir,"annobon_classificacion_clump.tif")
))

##########################################################################################
#### Read and create a reclass code for clumps > 50 pixels
df <- read.table(paste0(rootdir,"stat_annobon_classificacion_clump.txt"))[,1:2]

names(df) <- c("clump_id","size")
summary(df)

df$new <- 1
df[df$size > 50,]$new <- 2

table(df$new)

write.table(df,paste0(rootdir,"annobon_clump_gt_50.txt"),row.names = F,col.names = F)

####### Reclassificar para la masquera
system(sprintf("(echo %s; echo 1; echo 1; echo 3; echo 0) | oft-reclass  -oi %s  -um %s %s",
               paste0(rootdir,"annobon_clump_gt_50.txt"),
               paste0(rootdir,"annobon_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"annobon_classificacion_clump.tif"),
               paste0(rootdir,"annobon_classificacion_clump.tif")

))

####### Segmentacion del imagen
system(sprintf("(echo 0; echo 0 ; echo 0)|oft-seg -region -ttest -automax %s %s",
               "../GFC_hansen_LSAT_mosaics_2014/Ann_14046.tif",
               paste0(rootdir,"seg_annobon_gfc.tif")
)
)

####### Reclip segments
system(sprintf("oft-clip.pl %s %s %s",
               paste0(rootdir,"annobon_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"seg_annobon_gfc.tif"),
               paste0(rootdir,"seg_annobon_gfc_clip.tif")
)
)

####### Reclip GFC_perdidas cleaned
system(sprintf("oft-clip.pl %s %s %s",
               paste0(rootdir,"annobon_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"annobon_gfc_clean_nd.tif"),
               paste0(rootdir,"annobon_gfc_clean_nd_clip.tif")
)
)

####### Apply mask to identify segments from both products
system(sprintf("gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc=\"%s\"",
               paste0(rootdir,"seg_annobon_gfc_clip.tif"),
               paste0(rootdir,"annobon_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"annobon_classificacion_clump.tif"),
               paste0(rootdir,"tmp_annobon_clump_mixed.tif"),
               "(B>0)*((B==1)*C+(B==2)*A)"
))


#################### SIEVE RESULTS x2
system(sprintf("gdal_sieve.py -st %s %s %s",
               2,
               paste0(rootdir,"tmp_annobon_clump_mixed.tif"),
               paste0(rootdir,"tmp_annobon_clump_mixed_sieve02.tif")
))

#################### SIEVE RESULTS x4
system(sprintf("gdal_sieve.py -st %s %s %s",
               4,
               paste0(rootdir,"tmp_annobon_clump_mixed_sieve02.tif"),
               paste0(rootdir,"tmp_annobon_clump_mixed_sieve04.tif")
))

#################### SIEVE RESULTS x8
system(sprintf("gdal_sieve.py -st %s %s %s",
               8,
               paste0(rootdir,"tmp_annobon_clump_mixed_sieve04.tif"),
               paste0(rootdir,"tmp_annobon_clump_mixed_sieve08.tif")
))


#################### SIEVE RESULTS x12
system(sprintf("gdal_sieve.py -st %s %s %s",
               12,
               paste0(rootdir,"tmp_annobon_clump_mixed_sieve08.tif"),
               paste0(rootdir,"tmp_annobon_clump_mixed_sieve12.tif")
))
####### Reclump results
system(sprintf("oft-clump -i %s -o %s -um %s",
               paste0(rootdir,"tmp_annobon_clump_mixed_sieve12.tif"),
               paste0(rootdir,"annobon_clump_mix.tif"),
               paste0(rootdir,"tmp_annobon_clump_mixed_sieve12.tif")
))

####### Compute zonal stats for LOSSES
system(sprintf("oft-zonal_large_list.py -i %s -o %s -um %s",
               paste0(rootdir,"annobon_gfc_clean_nd_clip.tif"),
               paste0(rootdir,"zonal_perdidas_annobon.txt"),
               paste0(rootdir,"annobon_clump_mix.tif")
))

####### Compute zonal stats for CLASSIFICATION
system(sprintf("oft-zonal_large_list.py -i %s -o %s -um %s",
               paste0(rootdir,"annobon_classificacion_nd.tif"),
               paste0(rootdir,"zonal_classif_annobon.txt"),
               paste0(rootdir,"annobon_clump_mix.tif")
))

####### Polygonize the clump results
# system(sprintf("gdal_polygonize.py -8 -f \"ESRI Shapefile\" %s %s",
#                paste0(rootdir,"annobon_clump_mix.tif"),
#                paste0(rootdir,"annobon_clump_mix.shp")
# ))


####### Read zonal stats and rename columns
df_perdidas <- read.table(paste0(rootdir,"zonal_perdidas_annobon.txt"))
df_classif  <- read.table(paste0(rootdir,"zonal_classif_annobon.txt"))

names(df_classif)  <- c("clump_id","total","nodata",paste0("class_",1:7))

codes_gfclimpio <- c(1,2,3,4,5,6,7,8,9,10,105,106,107,108,109,110,111,112,113,114,190,195)
df_perdida2 <- df_perdidas[,c(1,2,3,3+codes_gfclimpio)]
names(df_perdida2) <- c("clump_id","total","nodata",paste0("perd_",codes_gfclimpio))

codes_bosque  <- paste0("perd_",c(1,2,3,4,5,190,195))
codes_nobosq  <- paste0("perd_",c(6,7,8,9))
codes_agua    <- paste0("perd_",c(10))
codes_perd_bs <- paste0("perd_",c(105,106,107,108,109,110,111,112,113,114))

dbf1 <- cbind(df_classif[,c(1,4:10)],df_perdida2)
head(dbf1)
# ####### Read DBF and merge with Zonal stats
# dbf <- read.dbf(paste0(rootdir,"bioko_clump_mix_bckup.dbf"))
# dbf$unique <- row(dbf)[,1]
# 
# dbf1 <- merge(dbf,df_perdida2,by.x="DN",by.y="clump_id",all.x=TRUE)
# dbf1 <- merge(dbf1,df_classif[,c(1,4:7)],by.x="DN",by.y="clump_id",all.x=TRUE)
# 
# dbf1 <- arrange(dbf1,unique)
# dbf1[is.na(dbf1)] <- 0

dbf1$new_class <- max.col(dbf1[,paste0("class_",1:7)])

dbf1$new_perd  <- 1

dbf1[rowSums(dbf1[,codes_bosque]) <= 0.3*dbf1$total & rowSums(dbf1[,codes_perd_bs]) > 0,]$new_perd <- 2
dbf1[rowSums(dbf1[,codes_bosque]) >  0.3*dbf1$total & rowSums(dbf1[,codes_perd_bs]) > 0,]$new_perd <- 3

dbf1[rowSums(dbf1[,codes_bosque]) <= 0.3*dbf1$total & rowSums(dbf1[,codes_perd_mg]) > 0,]$new_perd <- 4
dbf1[rowSums(dbf1[,codes_bosque]) >  0.3*dbf1$total & rowSums(dbf1[,codes_perd_mg]) > 0,]$new_perd <- 5

dbf1[rowSums(dbf1[,codes_nobosq]) > 0.9*dbf1$total,]$new_perd <- 6
dbf1[dbf1[,codes_agua]            > 0.9*dbf1$total,]$new_perd <- 7

dbf1$combi <- dbf1$new_perd*10+dbf1$new_class

table(dbf1$combi)

write.table(dbf1[,c("clump_id","total","combi")],paste0(rootdir,"annobon_reclass.txt"),row.names = F,col.names = F)

####### Reclassificar para la masquera
system(sprintf("(echo %s; echo 1; echo 1; echo 3; echo 0) | oft-reclass  -oi %s  -um %s %s",
               paste0(rootdir,"annobon_reclass.txt"),
               paste0(rootdir,"tmp_annobon_clump_mix_reclass.tif"),
               paste0(rootdir,"annobon_clump_mix.tif"),
               paste0(rootdir,"annobon_clump_mix.tif")
               
))

####### Reproject in Albers
system(sprintf("gdalwarp -t_srs %s -ot Byte -co COMPRESS=LZW %s %s",
               "../albers_equal_area_projection.txt",
               paste0(rootdir,"tmp_annobon_clump_mix_reclass.tif"),
               paste0(rootdir,"aea_annobon_map_dd_20171106.tif")
)
)


####### Comprimir el resultado
system(sprintf("gdal_translate -ot Byte -co COMPRESS=LZW %s %s",
               paste0(rootdir,"tmp_annobon_clump_mix_reclass.tif"),
               paste0(rootdir,"geo_annobon_map_dd_20171106.tif")
)
)

system(sprintf("rm %s",
               paste0(rootdir,"tmp*")))


###############################################################
############################      COMPUTE AREAS FOR AEA MAPS
###############################################################

system(sprintf("oft-stat -i %s -o %s -um %s",
               paste0(rootdir,"aea_annobon_map_dd_20171106.tif"),
               paste0(rootdir,"aea_annobon_stat_dd_20171106.txt"),
               paste0(rootdir,"aea_annobon_map_dd_20171106.tif")
))

system(sprintf("oft-stat -i %s -o %s -um %s",
               paste0(rootdir,"aea_bioko_map_dd_20171106.tif"),
               paste0(rootdir,"aea_bioko_stat_dd_20171106.txt"),
               paste0(rootdir,"aea_bioko_map_dd_20171106.tif")
))

system(sprintf("oft-stat -i %s -o %s -um %s",
               paste0(rootdir,"aea_continente_map_dd_20171106.tif"),
               paste0(rootdir,"aea_continente_stat_dd_20171106.txt"),
               paste0(rootdir,"aea_continente_map_dd_20171106.tif")
))

stat_bioko <- read.table(paste0(rootdir,"aea_bioko_stat_dd_20171106.txt"))
stat_conti <- read.table(paste0(rootdir,"aea_continente_stat_dd_20171106.txt"))
stat_annob <- read.table(paste0(rootdir,"aea_annobon_stat_dd_20171106.txt"))


pix_a <- res(raster(paste0(rootdir,"aea_annobon_map_dd_20171106.tif")))[1]
pix_b <- res(raster(paste0(rootdir,"aea_bioko_map_dd_20171106.tif")))[1]
pix_c <- res(raster(paste0(rootdir,"aea_continente_map_dd_20171106.tif")))[1]

df_a <- tapply(stat_annob$V2,stat_annob$V3,sum)*pix_a*pix_a/10000
df_b <- tapply(stat_bioko$V2,stat_bioko$V3,sum)*pix_b*pix_b/10000
df_c <- tapply(stat_conti$V2,stat_conti$V3,sum)*pix_c*pix_c/10000

sum(df_a)+sum(df_b)+sum(df_c)

df_a <- data.frame(df_a)
names(df_a) <- "area_annobon"
df_a$class <- rownames(df_a)

df_b <- data.frame(df_b)
names(df_b) <- "area_bioko"
df_b$class <- rownames(df_b)

df_c <- data.frame(df_c)
names(df_c) <- "area_continente"
df_c$class <- rownames(df_c)

df <- merge(df_a,df_b,all=T)
df <- merge(df,df_c,all=T)

write.csv(df,"stats_gnq_dd_20171106.csv")
