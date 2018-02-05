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
workdir    <- "/media/dannunzio/OSDisk/Users/dannunzio/Documents/countries/guinea_ecuatorial/gis_data/degradation_option2_por_agente_poligonos/uni_version_20171108/"
setwd(workdir)


#### Replace 255 by NONE in UNI FILES
system(sprintf("gdal_translate -a_nodata none -co COMPRESS=LZW %s %s",
               paste0(workdir,"co_nd_uni_aea_continente_map_dd_20171106.tif"),
               paste0(workdir,"tmp_c.tif")
))

system(sprintf("gdal_translate -a_nodata none -co COMPRESS=LZW %s %s",
               paste0(workdir,"co_nd_uni_aea_bioko_map_dd_20171106.tif"),
               paste0(workdir,"tmp_b.tif")
))

system(sprintf("gdal_translate -a_nodata none -co COMPRESS=LZW %s %s",
               paste0(workdir,"co_nd_uni_aea_annobon_map_dd_20171106.tif"),
               paste0(workdir,"tmp_a.tif")
))


##########################################################################################
#### Rasterize boundaries
system(sprintf("oft-rasterize_attr.py -v %s -i %s -o %s -a %s",
               "../../shapefiles_GNQ/zone_c_contour_aea.shp",
               "co_nd_uni_aea_continente_map_dd_20171106.tif",
               "mask_c.tif",
               "id"
               ))

system(sprintf("oft-rasterize_attr.py -v %s -i %s -o %s -a %s",
               "../../shapefiles_GNQ/zone_a_contour_aea.shp",
               "co_nd_uni_aea_annobon_map_dd_20171106.tif",
               "mask_a.tif",
               "id"
))

system(sprintf("oft-rasterize_attr.py -v %s -i %s -o %s -a %s",
               "../../shapefiles_GNQ/zone_b_contour_aea.shp",
               "co_nd_uni_aea_bioko_map_dd_20171106.tif",
               "mask_b.tif",
               "id"
))

##########################################################################################
#### Create binary mask
system(sprintf("gdal_calc.py -A %s  --type=Byte --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               "mask_a.tif",
               "mask_a_bin.tif",
               "(A>0)"
))

system(sprintf("gdal_calc.py -A %s  --type=Byte --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               "mask_b.tif",
               "mask_b_bin.tif",
               "(A>0)"
))

system(sprintf("gdal_calc.py -A %s  --type=Byte --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               "mask_c.tif",
               "mask_c_bin.tif",
               "(A>0)"
))

##########################################################################################
#### Shrink by morphological filter
size_morpho <- 10

system(sprintf("otbcli_BinaryMorphologicalOperation -in %s -out %s -structype.ball.xradius %s -structype.ball.yradius %s -filter %s",
               "mask_c_bin.tif",
               paste0("mask_shrink",size_morpho,"_c.tif"),
               size_morpho,
               size_morpho,
               "erode"
))

system(sprintf("otbcli_BinaryMorphologicalOperation -in %s -out %s -structype.ball.xradius %s -structype.ball.yradius %s -filter %s",
               "mask_b_bin.tif",
               paste0("mask_shrink",size_morpho,"_b.tif"),
               size_morpho,
               size_morpho,
               "erode"
))

system(sprintf("otbcli_BinaryMorphologicalOperation -in %s -out %s -structype.ball.xradius %s -structype.ball.yradius %s -filter %s",
               "mask_a_bin.tif",
               paste0("mask_shrink",size_morpho,"_a.tif"),
               size_morpho,
               size_morpho,
               "erode"
))

##########################################################################################
#### Create a two layers boundary
system(sprintf("gdal_calc.py -A %s -B %s --type=Byte --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               "mask_a_bin.tif",
               paste0("mask_shrink",size_morpho,"_a.tif"),
               "mask_diff_a.tif",
               "A+B"
))

system(sprintf("gdal_calc.py -A %s -B %s --type=Byte --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               "mask_b_bin.tif",
               paste0("mask_shrink",size_morpho,"_b.tif"),
               "mask_diff_b.tif",
               "A+B"
))

system(sprintf("gdal_calc.py -A %s -B %s --type=Byte --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               "mask_c_bin.tif",
               paste0("mask_shrink",size_morpho,"_c.tif"),
               "mask_diff_c.tif",
               "A+B"
))

##########################################################################################
#### Replace NO DATA by WATER IN THE INNER BOUNDARIES AND LEAVE AS IT IS IN THE FRINGE
system(sprintf("gdal_calc.py -A %s -B %s --type=Byte --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               "tmp_c.tif",
               "mask_diff_c.tif",
               "tmp_uni_aea_map_dd_continente_20171206.tif",
               "(B==1)*((A<255)*A+(A==255)*0)+(B==2)*((A<255)*A+(A==255)*99)"
))

system(sprintf("gdal_calc.py -A %s -B %s --type=Byte --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               "tmp_b.tif",
               "mask_diff_b.tif",
               "tmp_uni_aea_map_dd_bioko_20171206.tif",
               "(B==1)*((A<255)*A+(A==255)*0)+(B==2)*((A<255)*A+(A==255)*99)"
))

system(sprintf("gdal_calc.py -A %s -B %s --type=Byte --outfile=%s --co=\"COMPRESS=LZW\" --calc=\"%s\"",
               "tmp_a.tif",
               "mask_diff_a.tif",
               "tmp_uni_aea_map_dd_annobon_20171206.tif",
               "(B==1)*((A<255)*A+(A==255)*0)+(B==2)*((A<255)*A+(A==255)*99)"
))

##########################################################################################
#### COMPRESS AND FINALIZE
system(sprintf("gdal_translate %s %s -ot byte -a_nodata none -co COMPRESS=LZW",
               "tmp_uni_aea_map_dd_annobon_20171206.tif",
               "uni_map_dd_annobon_aea_20171206.tif"))

system(sprintf("gdal_translate %s %s -ot byte -a_nodata none -co COMPRESS=LZW",
               "tmp_uni_aea_map_dd_bioko_20171206.tif",
               "uni_map_dd_bioko_aea_20171206.tif"))

system(sprintf("gdal_translate %s %s -ot byte -a_nodata none -co COMPRESS=LZW",
               "tmp_uni_aea_map_dd_continente_20171206.tif",
               "uni_map_dd_continente_aea_20171206.tif"))
