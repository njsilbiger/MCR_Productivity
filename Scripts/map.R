# load libraries ####
library(here)
library(tidyverse)
library(lubridate)
library(broom)
library(patchwork)
library(wesanderson)
library(ggmap)
library(sf)
library(osmdata)
library(units)
library(ggspatial)

# 1) Get Mo'orea polygon from OSM and pad it by ~3 km so tiles include the coast
moorea_poly <- getbb("Moorea, Windward Islands, French Polynesia",
                     format_out = "sf_polygon")

stopifnot(!is.null(moorea_poly))
moorea_sf <- st_as_sf(moorea_poly)

# Buffer in meters (work in Web Mercator to buffer in planar meters)
moorea_buf <- moorea_sf |>
  st_transform(3857) |>
  st_buffer(set_units(3000, "m")) |>   # change 3000 for more/less padding
  st_transform(4326)

# 2) Build bbox (left, bottom, right, top) guaranteed numeric & ordered
bb <- st_bbox(moorea_buf)
bbox <- c(
  left   = unname(bb["xmin"]),
  bottom = unname(bb["ymin"]),
  right  = unname(bb["xmax"]),
  top    = unname(bb["ymax"])
)

# 3) Clip to Web Mercator lat limits (+/-85.0511) to satisfy tile math
bbox["bottom"] <- max(bbox["bottom"], -85.0511)
bbox["top"]    <- min(bbox["top"],     85.0511)

# 4) Safety swaps if a numeric issue ever inverts the box
if (bbox["left"] > bbox["right"])  bbox[c("left","right")]   <- bbox[c("right","left")]
if (bbox["bottom"] > bbox["top"])  bbox[c("bottom","top")]   <- bbox[c("top","bottom")]

# ---- Basemap via ggmap (Stamen; no API key required) ----
basemap <- get_stadiamap(
  bbox = bbox,
  zoom = 12,                      # tweak if tiles look too coarse/fine
  maptype = "stamen_terrain_background",
  color = "color"
)


## LTER1 site
site <- tibble(lat = -17.484,long = -149.8333)


map1<-ggmap(basemap) +
  geom_point(data = site, aes(x = long, y = lat), size = 2.5, color = "firebrick")+
  labs(x = "Longtitude",
       y = "Latitude")+
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size=12),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "grey"),
    panel.background = element_blank()
  ) +
 
 # coord_sf(
 #   xlim = c(bbox["left"], bbox["right"]),
  #  ylim = c(bbox["bottom"], bbox["top"]),
   # expand = FALSE,
  #  crs = 4326
  #) +
  annotation_scale(location = "bl",  # "bl" for bottom left, other options: "br", "tl", "tr"
                   width_hint = 0.2) + # Suggested proportion of the plot area the scale bar occupies
  annotation_north_arrow(location = "tr", # "tr" for top right, other options: "br", "tl", "bl"
                         which_north = "true", 
                         width = unit(0.75, "cm"),
                         height = unit(0.75, "cm") # Always points to true north
  ) + coord_sf(crs = 4326)


ggsave(filename = here("Output","map.pdf"),plot = map1)
