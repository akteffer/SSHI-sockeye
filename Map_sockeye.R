## Map Sockeye manuscript 230324

#install.packages("bcmaps")
#install.packages("remotes")
#remotes::install_github("bcgov/bcmaps")

library(bcmaps)
available_layers()
library(remotes)
library(sf)
library(ggplot2)
library("sp")
library(grid)

world <- map_data("world")
dim(world)

us <- getData("GADM",country="USA",level=1)
canada <- getData("GADM",country="CAN",level=1)

# Load and plot the boundaries of B.C.

bc <- bc_bound()
plot(st_geometry(bc))


ggplot() + 
  geom_sf(data = bc_neighbours(), mapping = aes(fill = name)) + 
  geom_sf(data = bc_cities()) +
  coord_sf(datum = NA) +
  scale_fill_viridis_d(name = "Jurisdiction") +
  theme_minimal()


# Load watercourse data and plot with boundaries of B.C.
plot(bc_bound(class = "sp"))
plot(watercourses_15M(class = "sp"), add = TRUE)

## Start plotting sample locations on map
head(sw.data)
head(fw.data)
llsw <- sw.data[,c(22:23)]
dim(llsw)
llsw <- llsw[!duplicated(llsw), ]
llfw <- fw.data[,c(22:23)]
dim(llfw)
llfw <- llfw[!duplicated(llfw), ]
samplatlong <- rbind(llsw, llfw)
head(samplatlong) #this is the lat long of all samples in the study

# Load watercourse data and plot with boundaries of B.C.
plot(bc_bound(class = "sp"))
plot(watercourses_15M(class = "sp"), add = TRUE)

rect <- data.frame(
  x = c(-130, -115, -115, -130),
  y = c(47, 47, 53, 53)
)

insetmap <- ggplot() + 
  geom_sf(data = world) +
  geom_sf(data = bc_bound()) +
  xlim(-144,-52) +
  ylim(41,83) +
  labs(x=NULL,y=NULL) +
  geom_polygon(data=rect, aes(x, y, group = 1), alpha=.75)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.title=element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

sampmapsockeye <- ggplot() + 
  geom_sf(data = world) +
  geom_sf(data = bc_bound()) +
  xlim(-130,-115) +
  ylim(47,53) +
  geom_point(data = samplatlong, aes(y=Latitude, x=Longitude), shape=21, col="white", fill="black") 

jpeg(filename='figs/Fig_Sample Map.jpg', width=500, height=500, quality=300)
sampmapsockeye
dev.off()

jpeg(filename='figs/Fig_Sample Map.jpg', width=700, height=700, quality=300)
sampmapsockeye
print(insetmap, vp = viewport(0.835, 0.38, width = 0.3, height = 0.3))
dev.off()
