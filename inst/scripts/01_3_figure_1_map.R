######################
### Figure 1 - MAP ###
######################

# input  data: (updated)stations, limits
# output data: (updated)stations

# Clear workspace
rm(list = setdiff(ls(), c("background", "grid", "limits", "stations", "Y365")))

# Load packages
library("sf")
library("sp")
library("mapSpain")
library("tidyverse")
library("rnaturalearth")
library("rnaturalearthdata")

# Fix stations names
stations$NAME1 <- 
  c("BADAJOZ", "MADRID (RETIRO)", "MALAGA", "NAVACERRADA", "SALAMANCA", 
    "SAN SEBASTIAN", "TORTOSA", "VALENCIA", "ZARAGOZA", "BARCELONA (FABRA)",
    "ALBACETE", "BURGOS", "CIUDAD REAL", "CORUNA", "MURCIA", "SEVILLA", 
    "SORIA", "BILBAO", "SANTIAGO", "PONFERRADA", "LEON", "LOGRONO", "ZAMORA",
    "REUS", "BARCELONA (AEROPUERTO)", "MADRID (TORREJON)", "VITORIA", 
    "ALMERIA", "GIJON", "CACERES", "SANTANDER", "CASTELLON", "HUELVA",
    "LLEIDA", "MADRID (BARAJAS)", "MADRID (CUATROVIENTOS)", "MADRID (GETAFE)", 
    "MORON", "VALLADOLID", "DAROCA")

# Obtain features for the map 
hypsobath <- esp_get_hypsobath()
hypsobath <- st_transform(hypsobath, 2062)

bath_tints <- colorRampPalette( 
  rev( c( "#D8F2FE", "#C6ECFF", "#B9E3FF",
          "#ACDBFB", "#A1D2F7", "#96C9F0",
          "#8DC1EA", "#84B9E3", "#79B2DE",
          "#71ABD8" ) ) )

hyps_tints <- colorRampPalette( 
  rev( c( "#F5F4F2", "#E0DED8", "#CAC3B8", "#BAAE9A",
          "#AC9A7C", "#AA8753", "#B9985A", "#C3A76B",
          "#CAB982", "#D3CA9D", "#DED6A3", "#E8E1B6",
          "#EFEBC0", "#E1E4B5", "#D1D7AB", "#BDCC96",
          "#A8C68F", "#94BF8B", "#ACD0A5") ) )

levels <- sort(unique(hypsobath$val_inf))
br_bath <- length(levels[levels < 0])
br_terrain <- length(levels) - br_bath
pal <- c(bath_tints((br_bath)), hyps_tints((br_terrain) ))

hypsobath$val_inf[hypsobath$val_inf < 0] <- -1
hypsobath$val_inf <- as.factor(hypsobath$val_inf)
levels(hypsobath$val_inf)[1] <- "< 0"

aux.ggplot <- ggplot(hypsobath) +
  geom_sf(aes(fill = val_inf), color = NA) +
  coord_sf(xlim = st_coordinates(limits)[,1], 
           ylim = st_coordinates(limits)[,2]) + 
  scale_fill_manual(name = "Elevation", values = pal[c(7, 8:17)],
                    breaks = levels(hypsobath$val_inf),
                    guide = guide_legend(reverse = TRUE)) +
  xlab("Longitude") + ylab("Latitude") +
  geom_point(aes(x = X, y = Y), data = data.frame(st_coordinates(stations))) +
  ggrepel::geom_label_repel(aes(x = X, y = Y, label = stations$NAME1), 
                            nudge_y = 30000, label.size = 0.05,
                            max.time = 30, max.iter = 1000000, max.overlaps = 30,
                            data = data.frame(st_coordinates(stations)),
                            seed = 23)

world <- ne_countries(scale = "medium", returnclass = "sf")

map.europe <- ggplot(data = world) + 
  geom_sf(fill= "grey70") + 
  geom_rect(xmin = -10, xmax = 4, ymin = 35.5, ymax = 44, 
            fill = NA, colour = "black") +
  coord_sf(xlim = c(-11, 41), ylim = c(35, 70), expand = FALSE) + 
  xlab("") + ylab("") + 
  theme(panel.background = element_rect(fill = "grey95"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

limits2 <- data.frame(X = c(0, 4), Y = c(35.5, 35.5 + 3.15) - .25)
limits2 <- st_transform(
  as(
    SpatialPointsDataFrame(
      coords = limits2, 
      data = limits2,
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
    'sf'
  ),
  2062
)

aux.ggplot <- aux.ggplot + 
  annotation_custom(ggplotGrob(map.europe), 
                    xmin = st_coordinates(limits2)[1,1], 
                    xmax = st_coordinates(limits2)[2,1], 
                    ymin = st_coordinates(limits2)[1,2], 
                    ymax = st_coordinates(limits2)[2,2]) +
  geom_rect(xmin = st_coordinates(limits2)[1,1], 
            xmax = st_coordinates(limits2)[2,1], 
            ymin = st_coordinates(limits2)[1,2], 
            ymax = st_coordinates(limits2)[2,2], 
            fill = NA, colour = "black", size = 1.5)

# Write figure
ggplot2::ggsave("inst/img/MAIN_map.png", aux.ggplot, 
                width = 1.2 * 8.27, height = 1.2 * 11.69 / 2)
