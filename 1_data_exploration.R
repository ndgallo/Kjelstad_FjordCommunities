## FjordCommunities
##
## Data exploration (updated 23 March 2024)

## Code sections overview ####
##  Map of stations - Figure 1
##  Maps of spatial distributions of env. variables - Figure 3
##  Maps of Fish and Crustacean CPUE and Periphylla CPUE - Figure 4
##  Map of stations and diversity - Figure 5

##  Comparison of trawl effect on biodiversity and biomass
##  Boxplot trawl comparison - Supplement 11

rm(list = ls()) #clear workspace

## Load relevant packages and data ####
source("0_setup.R")
load("_data/Statistical_analysis.rda")

### 1: Maps

### Map of stations - Figure 1 ####
map_stations <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = catch_df, aes(
      x = longitudestart, y = latitudestart,
      color = factor(startyear)
    ),
    size = 2.5, alpha = 0.7
  ) +
  scale_color_manual(
    values = c(
      "seagreen",
      "cyan1",
      "yellow",
      "royalblue4",
      "deeppink1",
      "darkolivegreen1",
      "darkorchid",
      "darkseagreen",
      "deeppink4",
      "darkgoldenrod1"
    ),
    name = "Year"
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 4)))
map_stations

map_Mas <- basemap(limits = c(5.25, 5.53, 60.81, 60.91), land.col = "gray90") +
  geom_spatial_point(
    data = catch_df, aes(
      x = longitudestart, y = latitudestart,
      color = factor(startyear)
    ),
    size = 2.2, alpha = 0.5
  ) +
  scale_color_manual(
    values = c(
      "seagreen",
      "cyan1",
      "yellow",
      "royalblue4",
      "deeppink1",
      "darkolivegreen1",
      "darkorchid",
      "darkseagreen",
      "deeppink4",
      "darkgoldenrod1"
    ),
    name = "Year"
  ) + rremove("xylab") + theme(legend.position = "none")
map_Mas

map_nor <- basemap(limits = c(3.8, 12, 58, 64), land.col = "gray90") + rremove("xylab")
map_nor

# combine
inset <- ggdraw(xlim = c(0, 60), ylim = c(0, 60)) +
  draw_plot(map_stations, x = 5, y = 0, width = 30, height = 60) +
  draw_plot(map_nor, x = 34.5, y = 21.5, width = 22, height = 35) +
  draw_plot(map_Mas, x = 37.2, y = 15.5, width = 19.3, height = 15)
inset
ggsave(inset, filename = "stations_map.png", width = 200, height = 250, unit = "mm")

# add names and arrow
finished_map <- inset + geom_segment(aes(x = 18, y = 34, xend = 37.5, yend = 25),
  arrow = arrow(length = unit(0.02, "npc")), size = 0.5
) +
  geom_rect(aes(xmin = 39, xmax = 45.5, ymin = 34, ymax = 42.5),
    fill = "cyan1", alpha = 0.25
  ) +
  annotate(geom = "text", x = 24, y = 36, label = "Sognefjord", color = "black") +
  annotate(geom = "text", x = 24, y = 45, label = "Nordfjord", color = "black") +
  annotate(geom = "text", x = 41, y = 24.5, label = "Masfjord", color = "black") +
  annotate(geom = "text", x = 28, y = 23, label = "Bømlafjord", color = "black") +
  annotate(geom = "text", x = 18.5, y = 29.5, label = "Bergen", color = "black") +
  annotate(geom = "text", x = 22, y = 41.3, label = "Førdefjord", color = "black") +
  annotate(geom = "text", x = 31.6, y = 42.7, label = "Lustrafjord", color = "black")

finished_map
ggsave(finished_map, filename = "Fig.1 Stations_map.png", width = 200, height = 250, unit = "mm", dpi = 600)

# map sill category
# catch_df$sill_category <- env_df$sill_category

# map_sillcat <- basemap(limits=c(4.5,7.5,59.4,62.1),land.col="gray90") +
# geom_spatial_point(data= catch_df,aes(x=longitudestart,y=latitudestart,
#                                      color=factor(sill_category)),
#                     size = 2.5, alpha=0.8) +
# scale_color_manual(values = c("seagreen",
#                              "deeppink1",
#                             "royalblue4",
#                            "darkorchid",
#                           "cyan1"),
#                name="Sill category")+
#  xlab("Longitude") + ylab("Latitude") +
# theme(text=element_text(size=14), legend.position="bottom")+
# guides(color = guide_legend(override.aes = list(size = 4)))
# map_sillcat

### Maps of spatial distributions of env. variables - Figure 3 ####

## A) Sill depth (fjord) or bottom depth (coastal), B) Station bottom depth, C) Temperature, D) ##Salinity, E) Oxygen, F) Aquaculture impact

bottomdepth_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = bottomdepth,
    shape = fjord_coast,
    size = 1
  )) +
  labs(title = "Bottom depth", tag = "a") +
  scale_color_gradientn(
    colors = c("cyan1", "darkslategray4", "royalblue4"),
    limits = c(100, 700),
    name = "(m)"
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 12)) +
  guides(size = "none", alpha = "none", shape = guide_legend("Station location"))
bottomdepth_map

ggsave("bottomdepth_map.tiff", plot = bottomdepth_map, width = 6, height = 4, units = "in", dpi = 600)

env_mod$fjord_sill <- FjordCommunities_variables_static$sill_depth_m
silldepth_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = fjord_sill,
    size = 0.1
  )) +
  labs(title = "Sill depth (fjord stations)", tag = "b") +
  scale_color_gradientn(
    colors = c("darkgoldenrod1", "deeppink", "deeppink4"),
    limits = c(30, 400),
    name = "(m)"
  ) +
  scale_shape_manual(values = c(15, 16)) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 12)) +
  guides(size = "none", alpha = "none")
silldepth_map
env_mod <- subset(env_mod, select = -(fjord_sill))
ggsave("silldepth_map.tiff", plot = silldepth_map, width = 6, height = 4, units = "in", dpi = 600)

temperature_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = Temperature,
    size = 0.1
  )) +
  labs(title = "Temperature", tag = "c") +
  scale_color_gradientn(
    colors = c("cyan1", "slateblue", "tomato"),
    limits = c(7, 9),
    name = "(°C)"
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 12)) +
  guides(size = "none", alpha = "none")
temperature_map
ggsave("temperature_map.tiff", plot = temperature_map, width = 6, height = 4, units = "in", dpi = 600)

salinity_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = Salinity,
    size = 0.1
  )) +
  labs(title = "Salinity", tag = "d") +
  scale_color_gradientn(
    colors = c("darkgoldenrod1", "darkslategray3", "royalblue4"),
    limits = c(33, 35.5),
    name = "(PSU)"
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 12)) +
  guides(size = "none", alpha = "none")
salinity_map
ggsave("salinity_map.tiff", plot = salinity_map, width = 6, height = 4, units = "in", dpi = 600)

oxygen_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = Oxygen,
    size = 0.1
  )) +
  labs(title = "Oxygen", tag = "e") +
  scale_color_gradientn(
    colors = c("tomato", "orchid", "royalblue4"),
    limits = c(2, 6),
    name = "(ml/L)"
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 12)) +
  guides(size = "none", alpha = "none")
oxygen_map
ggsave("oxygen_map.tiff", plot = oxygen_map, width = 6, height = 4, units = "in", dpi = 600)

aqua_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = aquaculture_impact,
    size = 0.1
  )) +
  labs(title = "Aquaculture impact score", tag = "f") +
  scale_color_gradientn(
    colors = c("darkslategray3", "goldenrod1", "tomato"),
    limits = c(0, 4100),
    name = ""
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 12)) +
  guides(size = "none", alpha = "none")
aqua_map
ggsave("aqua_map.tiff", plot = aqua_map, width = 6, height = 4, units = "in", dpi = 600)

## a) Station bottom depth, b) Sill depth (fjord) c) Temperature, d) Salinity, e) Oxygen, f) Aquaculture impact
library(cowplot)
multipanel <- plot_grid(bottomdepth_map, silldepth_map, temperature_map,
                        salinity_map, oxygen_map, aqua_map, labels = "")
multipanel

env_maps1 <- ggarrange(bottomdepth_map, silldepth_map, temperature_map,
                       salinity_map, oxygen_map, aqua_map,
                       ncol = 3, nrow = 2, common.legend = F,
                       heights = c(1, 1)
)
env_maps1
tiff("Figure3.tiff", width = 7, height = 5, units = "in", res = 300)
plot(env_maps1)
dev.off()

ggsave(env_maps1, filename = "env_maps.png", width = 7, height = 5)
env_maps2 <- ggarrange(oxygen_map, aqua_map, ncol = 2, nrow = 1, common.legend = F)
ggsave(env_maps2, filename = "env_maps2.png", width = 9, height = 5)

### Maps of Fish and Crustacean CPUE and Periphylla CPUE - Figure 4 ####
## tot catch minus periphylla
library(readxl)
library(dplyr)

comm_matrix <- read_excel("_data/FjordCommunities_community_matrix.xlsx", sheet = "Community matrix")
comm_matrix_minus_Periphylla <- comm_matrix %>% select(-catchweight_g_Periphylla, -ID)
CPUE_fish_and_crust_g <- rowSums(comm_matrix_minus_Periphylla)
CPUE_fish_and_crust_kg <- CPUE_fish_and_crust_g/1000
CPUE_Periphylla_g <- comm_matrix %>% select(catchweight_g_Periphylla)
CPUE_Periphylla_kg <- CPUE_Periphylla_g/1000
summary(CPUE_fish_and_crust_kg)
summary(CPUE_Periphylla_kg)

env_df <- read_excel("_data/FjordCommunities_env_df.xlsx", sheet = "Environmental data")
env_df$catchweight_tot_minusperiphylla_kg <- CPUE_fish_and_crust_kg
env_df$catchweight_kg_periphylla <- CPUE_Periphylla_kg$catchweight_g_Periphylla

tot_catch_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = env_df, aes(
      x = longitudestart,
      y = latitudestart,
      color = Trawl,
      size = catchweight_tot_minusperiphylla_kg
    ),
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "royalblue3",
    "deeppink4"
  )) +
  labs(title = "Total Catch minus Periphylla", tag = "a") +
  theme(text = element_text(size = 12)) +
  guides(
    size = guide_legend(title = "CPUE (kg/min)"),
    color = guide_legend(override.aes = list(size = 4))
  )
tot_catch_map

## periphylla catches
periphylla_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = env_df, aes(
      x = longitudestart,
      y = latitudestart,
      color = Trawl,
      size = catchweight_kg_periphylla
    ),
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "royalblue3",
    "deeppink4"
  )) +
  labs(title = "Catch Periphylla", tag = "b") +
  theme(text = element_text(size = 12)) +
  guides(
    size = guide_legend(title = "CPUE (kg/min)"),
    color = guide_legend(override.aes = list(size = 4))
  )
periphylla_map

catch_maps <- ggarrange(tot_catch_map, periphylla_map, ncol=2, nrow=1)
catch_maps
ggsave("Figure4.tiff", plot = catch_maps, width = 9, height = 5, units = "in", dpi = 600)

### Map of stations and diversity - Figure 5 ####
map_diversity <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = catch_df, aes(
      x = longitudestart, y = latitudestart,
      color = shannon_div,
      shape = Trawl
    ),
    size = 4
  ) +
  scale_color_gradientn(
    colors = c("deeppink4", "darkgoldenrod1", "royalblue3"),
    limits = c(0, 2.5),
    name = "Shannon-Wiener \n diversity index"
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 16))
map_diversity
ggsave("Figure5.tiff", plot = map_diversity, width = 10, height = 8, units = "in", dpi = 600)

### Comparison of trawl effect on biodiversity and biomass ####

### 4.1: Trawl effect on biomass 
#use env_df for biomass and env_mod for diversity
env_df$sill_category <- as.factor(env_df$sill_category) # categorical variable - factor
env_df$Trawl <- as.factor(env_df$Trawl) # categorical variable - factor

trawl_biomass <- ggplot(env_df, aes(
  x = Trawl, y = catchweight_tot_minusperiphylla_kg,
  fill = Trawl
)) +
  scale_fill_manual(values = alpha(c("royalblue3", "deeppink4"), .7)) +
  xlab("Trawl") +
  ylab("Fish and crustaceans catch rate (kg/min)") +
  geom_boxplot() +
  #stat_summary(
    #fun.data = get_box_stats, 
    #geom = "text", hjust = 0.5, vjust = -0.5) +
  labs(tag = "a") +
  theme_minimal() +
  theme(legend.position = "none")
trawl_biomass

summary_FishandCrusCPUE <- env_df %>%
  group_by(Trawl) %>%
  summarise(
    Mean = mean(catchweight_tot_minusperiphylla_kg, na.rm = TRUE),
    Median = median(catchweight_tot_minusperiphylla_kg, na.rm = TRUE)
  )
print(summary_FishandCrusCPUE)

trawl_periphylla <- ggplot(env_df, aes(
  x = Trawl, y = catchweight_kg_periphylla,
  fill = Trawl
)) +
  scale_fill_manual(values = alpha(c("royalblue3", "deeppink4"), .7)) +
  xlab("Trawl") +
  ylab("Periphylla catch rate (kg/min)") +
  geom_boxplot() +
  #stat_summary(fun.data = get_box_stats, geom = "text", hjust = 1, vjust = -2) +
  labs(tag = "b") +
  theme_minimal() +
  theme(legend.position = "none")
trawl_periphylla

summary_PeriphyllaCPUE <- env_df %>%
  group_by(Trawl) %>%
  summarise(
    Mean = mean(catchweight_kg_periphylla, na.rm = TRUE),
    Median = median(catchweight_kg_periphylla, na.rm = TRUE)
  )
print(summary_PeriphyllaCPUE)

### 4.2: Trawl effect on diversity 
trawl_diversity <- ggplot(env_mod, aes(
  x = Trawl, y = shannon_div,
  fill = Trawl
)) +
  scale_fill_manual(values = alpha(c("royalblue3", "deeppink4"), .7)) +
  xlab("Trawl") +
  ylab("Shannon-Wiener diversity") +
  geom_boxplot() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 1, vjust = 0.1) +
  labs(tag = "b") +
  theme_minimal() +
  theme(legend.position = "none")
trawl_diversity
# ggsave(trawl_boxplot, filename="trawl_boxplot.png", width=6)

summary_diversity <- env_mod %>%
  group_by(Trawl) %>%
  summarise(
    Mean = mean(shannon_div, na.rm = TRUE),
    Median = median(shannon_div, na.rm = TRUE)
  )
print(summary_diversity)

### Boxplot trawl comparison - Supplement 11 ####
trawl_catch <- ggarrange(trawl_biomass, trawl_diversity)
trawl_catch
ggsave("Supp11.tiff", plot = trawl_catch, width=8, height=5, units = "in", dpi = 600)
