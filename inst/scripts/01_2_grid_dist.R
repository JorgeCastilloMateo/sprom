###########################################################
### Section 2 - GRID of Spain and DISTANCE to the coast ###
###########################################################

# input  data: "data/stations.rda"
# output data: (updated)stations, background, grid, limits

# Clear workspace
rm(list = setdiff(ls(), c("stations", "Y365")))

# Load packages
library("sf")
library("sp")
library("tidyverse")
library("rnaturalearth")
library("rnaturalearthdata")

# Read list of stations of interest
load("data/stations.rda")

# Adding geometry with projection 2062
stations <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = stations[c("LON", "LAT")], 
      data = stations[c("STAID", "STANAME", "LON", "LAT", "HGHT")],
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)



# DEFINE GRID OF SPAIN
# Boundaring box of Spain
limits <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
      data = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)

# Load the polygon of Spain and projection 2062
spain <- ne_countries(scale = "large", country = "Spain", returnclass = "sf")
spain <- st_transform(spain, 2062)

# Extract peninsular Spain
spain_coords <- Polygons(
  list(Polygon(st_coordinates(spain)[st_coordinates(spain)[,"L2"] == 3, 1:2])),
  ID = "spain")
spain_coords <- SpatialPolygons(list(spain_coords))
spain_coords <- as(spain_coords, "sf")
st_crs(spain_coords) <- st_crs(spain)

# Build squared grid 10 x 10 km and intersect with peninsular Spain
grid <- st_make_grid(spain, cellsize = 10000, what = "centers")
grid <- st_intersection(grid, spain_coords)



# OBTAIN DISTANCE TO THE COAST (in km)
# Load the polygon of neighboring countries, join them, and projection 2062
world_map <- ne_countries(scale = "large", returnclass = 'sf')
european_union <- c("Algeria", "Andorra", "France", "Gibraltar", "Morocco", "Portugal", "Spain")
european_union_map <- 
  world_map %>% 
  filter(name %in% european_union)
background <- st_transform(european_union_map, 2062)
bbox_europe <- st_bbox(c(xmin = -10, ymin = 20, xmax = 50, ymax = 80), crs = st_crs(european_union_map))
european_union_map_cropped <- st_crop(european_union_map, bbox_europe)
european_union_map_cropped <- st_union(european_union_map_cropped)
european_union_map_cropped <- st_transform(european_union_map_cropped, 2062)

# Distance from each grid cell to the boundary
dist <- st_distance(st_boundary(european_union_map_cropped), grid)
grid <- st_sf(dist = round(as.vector(dist) / 1000, 3), geometry = grid)

# Plot of distance to the coast
ggplot(data = background) + 
  geom_sf(fill = "antiquewhite") + 
  xlab("Longitude") + ylab("Latitude") + ggtitle("Distance to the coast") +
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text.x=element_text(size = 6),
        axis.text.y=element_text(size = 6, angle = 90),
        axis.title=element_text(size = 10, face = "bold")) + 
  geom_tile(data = grid, ggplot2::aes(x = st_coordinates(grid)[, 1], y = st_coordinates(grid)[, 2], fill = dist)) +
  scale_fill_gradient2(midpoint = 50, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), 
                       space = "Lab", limits = c(0, 360), name = "Distance (km)") +
  coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2])

# Distance from the stations to the boundary
dist <- st_distance(st_boundary(european_union_map_cropped), stations$geometry)
stations$DIST <- round(as.vector(dist) / 1000, 3)
