####################################################################################
####### Object:  Segment and merge with existing classification               
####### Author:  remi.dannunzio@fao.org                               
####### Update:  2017/07/23                                      
####################################################################################
library(foreign)
library(plyr)
library(rgeos)
library(rgdal)
library(raster)
library(ggplot2)

options(stringsAsFactors = F)

rootdir <- "~/gnq_ws_20170726/data/"
setwd(rootdir)
start_time <- Sys.time()

##########################################################################################
#### Clump results
system(sprintf("oft-clump -i %s -o %s -um %s",
               paste0(rootdir,"bioko_classificacion_nd.tif"),
               paste0(rootdir,"tmp_bioko_classificacion_clump.tif"),
               paste0(rootdir,"bioko_classificacion_nd.tif")
))

##########################################################################################
#### Compute clump stats
system(sprintf("oft-stat -i %s -o %s -um %s",
               paste0(rootdir,"tmp_bioko_classificacion_clump.tif"),
               paste0(rootdir,"stat_bioko_classificacion_clump.txt"),
               paste0(rootdir,"tmp_bioko_classificacion_clump.tif")
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
               paste0(rootdir,"tmp_bioko_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"tmp_bioko_classificacion_clump.tif"),
               paste0(rootdir,"tmp_bioko_classificacion_clump.tif")
               
))

####### Segmentacion del imagen
system(sprintf("(echo 0; echo 0 ; echo 0)|oft-seg -region -ttest -automax %s %s",
               paste0(rootdir,"bioko_lsat_mosaic.tif"),
               paste0(rootdir,"tmp_seg_bioko_gfc.tif")
)
)

####### Dimension del producto final
e <- extent(raster(paste0(rootdir,"tmp_bioko_classificacion_clump_gt_50.tif")))

####### Reclip segmentos
system(sprintf("gdal_translate -a_nodata 0 -co COMPRESS=LZW -projwin %s %s %s %s %s %s",
               e@xmin,
               e@ymax,
               e@xmax,
               e@ymin,
               paste0(rootdir,"tmp_seg_bioko_gfc.tif"),
               paste0(rootdir,"seg_bioko_gfc_clip.tif")
))

####### Reclip producto
system(sprintf("gdal_translate -a_nodata 0 -co COMPRESS=LZW -projwin %s %s %s %s %s %s",
               e@xmin,
               e@ymax,
               e@xmax,
               e@ymin,
               paste0(rootdir,"bioko_gfc_clean_nd.tif"),
               paste0(rootdir,"bioko_gfc_clean_nd_clip.tif")
))


####### Apply mask to identify segments from both products
system(sprintf("gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc=\"%s\"",
               paste0(rootdir,"seg_bioko_gfc_clip.tif"),
               paste0(rootdir,"tmp_bioko_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"tmp_bioko_classificacion_clump.tif"),
               paste0(rootdir,"tmp_bioko_clump_mixed.tif"),
               "(B>0)*((B==1)*C+(B==2)*A)"
))


####### Reclump results
system(sprintf("oft-clump -i %s -o %s -um %s",
               paste0(rootdir,"tmp_bioko_clump_mixed.tif"),
               paste0(rootdir,"bioko_clump_mix.tif"),
               paste0(rootdir,"tmp_bioko_clump_mixed.tif")
))

####### Compute zonal stats for LOSSES
system(sprintf("../scripts/oft-zonal_large_list.py -i %s -o %s -um %s",
               paste0(rootdir,"bioko_gfc_clean_nd_clip.tif"),
               paste0(rootdir,"zonal_perdidas_bioko.txt"),
               paste0(rootdir,"bioko_clump_mix.tif")
))

####### Compute zonal stats for CLASSIFICATION
system(sprintf("../scripts/oft-zonal_large_list.py -i %s -o %s -um %s",
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

write.table(dbf1[,c("clump_id","total","combi")],paste0(rootdir,"bioko_reclass.txt"),row.names = F,col.names = F)

####### Reclassificar para la masquera
system(sprintf("(echo %s; echo 1; echo 1; echo 3; echo 0) | oft-reclass  -oi %s  -um %s %s",
               paste0(rootdir,"bioko_reclass.txt"),
               paste0(rootdir,"tmp_bioko_clump_mix_reclass.tif"),
               paste0(rootdir,"bioko_clump_mix.tif"),
               paste0(rootdir,"bioko_clump_mix.tif")
               
))



####### Comprimir el resultado
system(sprintf("gdal_calc.py -A %s -B %s --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               paste0(rootdir,"tmp_bioko_clump_mix_reclass.tif"),
               paste0(rootdir,"tmp_bioko_classificacion_clump_gt_50.tif"),
               paste0(rootdir,"bioko_clump_mix_reclass_20170726.tif"),
               "(B>0)*A"
)
)

####### Reproject in Albers
system(sprintf("gdalwarp -t_srs %s -ot Byte -overwrite -co COMPRESS=LZW %s %s",
               "albers_equal_area_projection.txt",
               paste0(rootdir,"bioko_clump_mix_reclass_20170726.tif"),
               paste0(rootdir,"aea_bioko_clump_mix_reclass_20170726.tif")
)
)

####### Limpiar la carpeta de los "tmp"
system(sprintf("rm %s",
               "tmp*.tif"))

####### Tiempo de processamiento
Sys.time()-start_time