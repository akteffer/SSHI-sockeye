install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))

library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")


library(raster)
states <- getData(country="USA", level=1)
provinces <- getData(country="Canada", level=1)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("World map", subtitle = paste0("(", length(unique(world$NAME)), " countries)"))

# BC only
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-131, -115), ylim = c(48, 53), expand = FALSE)+



library("ggspatial")
state_prov <- rnaturalearth::ne_states(c("united states of america", "canada"))
ggplot(data = world) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-135, -115), ylim = c(48, 60))

install.packages("bcmaps")
library(bcmaps)

bc <- bc_bound()
plot(st_geometry(bc))
plot(fw.data$Longitude, fw.data$Latitude)

library(sf)
library(ggplot2)
head(fw.data)
arena2.map <- fw.data[fw.data$arena2>=0,]

my_sf <- st_as_sf(my_df, coords = c('LON', 'LAT'))
my_sf <- st_set_crs(my_sf, crs = 4326)
                  
#Plot it:
ggplot(my_sf) + 
  geom_sf(aes(color = cluster))


library(ggplot2)
library(RgoogleMaps)
# loading the required packages
library(ggplot2)
library(ggmap)
library("plotGoogleMaps")
library(Rcurl)

# creating a sample data.frame with your lat/lon points
lon <- fw.data$Longitude
lat <- fw.data$Latitude
df <- as.data.frame(cbind(lon,lat))
center = c(mean(lat), mean(lon))  #tell what point to center on
zoom <- 2 #zoom: 1 = furthest out (entire globe), larger numbers = closer in
terrmap <- GetMap(center=center, zoom=zoom, type= "satellite", destfile = "satellite.png")

# plot points and save image

png('sockeyemap.png')
PlotOnStaticMap(terrmap, lat = lat, lon = lon, pch = 20, col = 'red')
dev.off()

