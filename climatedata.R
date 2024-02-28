## Climate variables and mapping ##
setwd("~/Desktop/Research/FireBiodiversity/Data")

library(sp)
library(sf)
library(raster)
library(geodata)

# Import site data, including lat/lon coordinates
sitedata <- read.csv("sitedata.csv")

# Switch Latitude/Longitude to x,y
sites <- sitedata[,3:2]
names(sites) <- c("x", "y")

# Set spatial coordinates and create spatial data frame
sts <- coordinates(sites)
sts <- SpatialPointsDataFrame(sts,data.frame(sav_sites = sitedata[,1]))

# Set coordinate reference system (CRS) to match worldclim data
crs(sts) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Repeat this process for savanna sites
sav_sites <- sitedata[sitedata$Biome == "Savanna",3:2]
names(sav_sites) <- c("x", "y")
sav_sts <- coordinates(sav_sites)
sav_sts <- SpatialPointsDataFrame(sav_sts,data.frame(sav_sites = sitedata[sitedata$Biome == "Savanna",1]))
crs(sav_sts) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# And forest sites
forest_sites <- sitedata[sitedata$Biome == "Forest",3:2]
names(forest_sites) <- c("x", "y")
forest_sts <- coordinates(forest_sites)
forest_sts <- SpatialPointsDataFrame(forest_sts,data.frame(forest_sites = sitedata[sitedata$Biome == "Forest",1]))
crs(forest_sts) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Use geodata library to extract worldclim bio data
r <- worldclim_global('worldclim', var='bio', res=10)
# And subset to the variables we are interested in
r <- r[[c(1,12, 15, 16, 17)]]
names(r) <- c("MAT", "MAP", "Seasonality", "WetQ", "DryQ")
values <- extract(r,st_as_sf(sts))

# Import and extract Global Aridity Index data from Trabucco and Zomer 2022
aridity_raster <- raster("ai_v3_yr.tif")
aridity <- extract(aridity_raster, st_as_sf(sts))

# Combine the climate data with the site data
clim_data <- cbind.data.frame(sitedata, values)
clim_data <- clim_data[,-7]
clim_data <- cbind.data.frame(clim_data, aridity)
write.csv(clim_data, file = "sitedata_full.csv", row.names = F)

# Mapping savanna and forest sites onto global map of MAP
colfunc <- colorRampPalette(c("white", "blue"))
plot(raster(r[[2]]), col=colfunc(30), zlim = c(0,4000), breaks = seq(0,4000,200))
# plot(aridity_raster, zlim = c(0,11000), col = colfunc(200))
plot(sav_sts, add = T, pch = 21, col = "black", bg = alpha("gold", 0.4))
plot(forest_sts,add=T, pch = 21, col = "black", bg = alpha("green4", 0.7))
