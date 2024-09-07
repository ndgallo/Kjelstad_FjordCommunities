## FjordCommunities
##
## Data exploration (updated 2 September 2024)

## Code sections overview ####
##  1: Maps: Stations and total catch + Periphylla 
##  1.1: Map of stations - Figure 1
##  1.2: Maps of total catch and periphylla - Figure 5

##  2: Common species in total catch and per trawl
##  2.1: Common species in total CPUE
##  2.2: Common species in relative CPUE, average per trawl
##  2.3: Maps of the catches of the most common species

##  3: Diversity 
##  3.1: Map of stations and diversity - Figure 6

##  4: Compare trawl effect on biodiversity and biomass
##  4.1: Trawl effect on biomass
##  4.2: Trawl effect on biodiversity
##  4.3: Trawl effect - Figure 11

##  5.1 Div/CPUE responses to sill depth, CPUE vs diversity
##  5.2 TS plot

##  6.1 Maps of spatial distributions of variables
##

## Load relevant packages and data ####
source("0_setup.R")
source("1_data_loading.R")

### 1: Maps: stations + total catch/periphylla ####

### 1.1: Stations - Figure 1 ####
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
  annotate(geom = "text", x = 41, y = 27, label = "Masfjord", color = "black") +
  annotate(geom = "text", x = 28, y = 23, label = "Bømlafjord", color = "black") +
  annotate(geom = "text", x = 18.5, y = 29.5, label = "Bergen", color = "black") +
  annotate(geom = "text", x = 22, y = 41.3, label = "Førdefjord", color = "black") +
  annotate(geom = "text", x = 31.6, y = 42.7, label = "Lustrafjord", color = "black")

finished_map
ggsave(finished_map, filename = "Fig.1 Stations_map.png", width = 200, height = 250, unit = "mm")

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


### 1.2: Map of total catch and Periphylla - Figure 5  ####
## tot catch minus periphylla
tot_catch_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = CPUE_catchweight, aes(
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
  labs(title = "Total Catch minus Periphylla", tag = "A") +
  theme(text = element_text(size = 12)) +
  guides(
    size = guide_legend(title = "CPUE (kg/min)"),
    color = guide_legend(override.aes = list(size = 4))
  )
tot_catch_map

## periphylla catches
periphylla_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = CPUE_catchweight, aes(
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
  labs(title = "Catch Periphylla", tag = "B") +
  theme(text = element_text(size = 12)) +
  guides(
    size = guide_legend(title = "CPUE (kg/min)"),
    color = guide_legend(override.aes = list(size = 4))
  )
periphylla_map

# catch_maps <- ggarrange(tot_catch_map, periphylla_map, ncol=2, nrow=1)
# ggsave(catch_maps, filename="catch_maps.png", width=9, height=5)

## summary totals
summary(CPUE_catchweight$catchweight_total)
summary(CPUE_catchweight$catchweight_tot_minusperiphylla_kg)
summary(CPUE_catchweight$catchweight_kg_periphylla)

CPUE_catchweight %>%
  arrange(desc(catchweight_total)) %>%
  slice(1)
# lowest catchweight_total = ID 72 -->
# highest catchweight_total= ID 32 -->

CPUE_catchweight %>%
  arrange(desc(catchweight_tot_minusperiphylla_kg)) %>%
  slice(1)

### 2: Common species in total catch and per trawl ####

## most common species
## 1) sum up abundance: absolute CPUE weight + rel abundances - find 9 with highest value
## 2) group together all other species to new column 'other' in new df
## 3) transform to long format
## 4) make plots

### 2.1: Most common species in total catch ####
## 1) sum absolute (total) CPUE
str(sort(colSums(CPUE_catchweight_clean[, 1:length(CPUE_catchweight_clean)]), decreasing = TRUE)[1:3])
str(sort(colSums(CPUE_catchweight_clean[, 1:length(CPUE_catchweight_clean)]), decreasing = TRUE)[4:6])
str(sort(colSums(CPUE_catchweight_clean[, 1:length(CPUE_catchweight_clean)]), decreasing = TRUE)[7:9])
str(sort(colSums(CPUE_catchweight_clean[, 1:length(CPUE_catchweight_clean)]), decreasing = TRUE)[10:12])
str(sort(colSums(CPUE_catchweight_clean[, 1:length(CPUE_catchweight_clean)]), decreasing = TRUE)[13:15])
str(sort(colSums(CPUE_catchweight_clean[, 1:length(CPUE_catchweight_clean)]), decreasing = TRUE)[16:18])
# 1-3: Periphylla[29], Coryphaenoides rupestris [42], Chimaera[3],
# 4-6: Etmopterus[5], Benthosema glaciale[27],  Micromesistius poutassou [12],
# 7-9:  Galeus melastomus[7],Argentina silus[1],  Squalus acanthias[25]
# 10-12:Dendrobranchiata[58], Pollachius virens[16], Lophius piscatorius [21]
# 13-15: Molva dypterygia [13], Dipturus oxyrinchus[33], Pollachius pollachius [15],
# 16-18: Merluccius merluccius [11], Melanogrammus aeglefinus [22], Glyptocephalus cynoglossus[]
# which( colnames(CPUE_catchweight_clean)=="catchweight_g_Squalus acanthias" )

sum(CPUE_catchweight$catchweight_total)
# periphylla = 793182 / 1214243 = 65 % --> 65% in total catch
793181 / 1214243


## 2) new df with common species + others column
catchweight_common_biomass <- CPUE_catchweight_clean %>%
  mutate("c) other species/taxa" = rowSums(CPUE_catchweight_clean[, -c(29, 42, 3, 5, 27, 12, 7, 1, 25)])) %>%
  select(
    "c) other species/taxa", catchweight_g_Periphylla,
    `catchweight_g_Coryphaenoides rupestris`,
    `catchweight_g_Chimaera monstrosa`,
    `catchweight_g_Etmopterus spinax`,
    `catchweight_g_Benthosema glaciale`,
    `catchweight_g_Micromesistius poutassou`,
    `catchweight_g_Galeus melastomus`,
    `catchweight_g_Argentina silus`,
    `catchweight_g_Squalus acanthias`
  )
colnames(catchweight_common_biomass)[2] <- "a) Periphylla periphylla"
colnames(catchweight_common_biomass)[3] <- "b) roundnose grenadier (Coryphaenoides rupestris)"
colnames(catchweight_common_biomass)[4] <- "d) rabbit fish (Chimaera monstrosa)"
colnames(catchweight_common_biomass)[5] <- "e) velvet-belly lanternshark (Etmopterus spinax)"
colnames(catchweight_common_biomass)[6] <- "f) glacier lantern fish (Benthosema glaciale)"
colnames(catchweight_common_biomass)[7] <- "g) blue whiting (Micromesistius poutassou)"
colnames(catchweight_common_biomass)[8] <- "h) blackmouth catshark (Galeus melastomus)"
colnames(catchweight_common_biomass)[9] <- "i) greater argentine (Argentina silus)"
colnames(catchweight_common_biomass)[10] <- "j) spiny dogfish (Squalus acanthias)"
# colnames(catchweight_common_biomass)[11]="pelagic shrimp spp. (Pasiphaea and Eusergestes)"

# get values in kg/min
catchweight_common_biomass <- catchweight_common_biomass %>%
  mutate_if(is.numeric, ~ . / 1000)

## 3) transform to long format
catchweight_common_biomass_long <- as.data.frame(catchweight_common_biomass)
catchweight_common_biomass_long$Station <- rownames(catchweight_common_biomass_long)
catchweight_common_biomass_long <- melt(catchweight_common_biomass_long, id.vars = "Station")

## 4) create plots
barplot_absolute <- ggplot(catchweight_common_biomass_long, aes(
  x = value, y = reorder(variable, +value),
  fill = variable
)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c(
    "darkgoldenrod1",
    "darkolivegreen1",
    "darkorchid",
    "olivedrab4",
    "deeppink1",
    "darkturquoise",
    "plum1",
    "royalblue3",
    "deeppink3",
    "cyan1"
  )) +
  theme_minimal() +
  rremove("ylab") +
  xlab("CPUE (kg/min)") +
  labs(title = "Total CPUE", tag = "A") +
  guides(fill = "none") +
  theme(text = element_text(size = 14))
barplot_absolute

# ggsave(barplot_absolute, filename="barplot_abs.png", width=13, height=7)

### 2.2 Most common species in catch per trawl ####
## 1) sum relative - common species
str(sort(colSums(catchweight_relative[, 1:length(catchweight_relative)]), decreasing = TRUE)[1:3])
str(sort(colSums(catchweight_relative[, 1:length(catchweight_relative)]), decreasing = TRUE)[4:6])
str(sort(colSums(catchweight_relative[, 1:length(catchweight_relative)]), decreasing = TRUE)[7:9])
str(sort(colSums(catchweight_relative[, 1:length(catchweight_relative)]), decreasing = TRUE)[10:12])
str(sort(colSums(catchweight_relative[, 1:length(catchweight_relative)]), decreasing = TRUE)[13:15])
# 1-3: Periphylla[29, 1564 g/min], Coryphaenoides rupestris[42, 1363], Chimaera[3, 945],
# 4-6: Etmopterus[5, 657], Galeus melastomus[7, 463], Micromesistius poutassou [12, 409],
# 7-9: Benthosema glaciale[27, 377], Argentina silus[1, 363],  Squalus acanthias[25, 282]
# 10-12:Pollachius virens[16, 278], Dendrobranchiata[58, 274], Merluccius merluccius [11, 198]
# 13-15: Melanogrammus aeglefinus [22, 191], Pollachius pollachius [15, 167],
#        Maurolicus muelleri [10, 130]
which(colnames(catchweight_relative) == "catchweight_g_Squalus acanthias")


## 2) new df with common species + others column
catchweight_common_rel <- catchweight_relative %>%
  mutate(Other = rowSums(catchweight_relative[, -c(
    29, 42, 3, 5, 7, 12, 27,
    1, 25
  )])) %>%
  select(
    Other, catchweight_g_Periphylla,
    `catchweight_g_Coryphaenoides rupestris`,
    `catchweight_g_Chimaera monstrosa`,
    `catchweight_g_Etmopterus spinax`,
    `catchweight_g_Galeus melastomus`,
    `catchweight_g_Micromesistius poutassou`,
    `catchweight_g_Benthosema glaciale`,
    `catchweight_g_Argentina silus`,
    `catchweight_g_Squalus acanthias`
  )

colnames(catchweight_common_rel)[2] <- "Periphylla periphylla"
colnames(catchweight_common_rel)[3] <- "roundnose grenadier (Coryphaenoides rupestris)"
colnames(catchweight_common_rel)[4] <- "rabbit fish (Chimaera monstrosa)"
colnames(catchweight_common_rel)[5] <- "velvet-belly lanternshark (Etmopterus spinax)"
colnames(catchweight_common_rel)[6] <- "blackmouth catshark (Galeus melastomus)"
colnames(catchweight_common_rel)[7] <- "blue whiting (Micromesistius poutassou)"
colnames(catchweight_common_rel)[8] <- "glacier lantern fish (Benthosema glaciale)"
colnames(catchweight_common_rel)[9] <- "greater argentine (Argentina silus)"
colnames(catchweight_common_rel)[10] <- "spiny dogfish (Squalus acanthias)"

colSums(catchweight_common_rel)
# other = 2291 (26%), periphylla = 1564(18%), round = 1450 (16%), 945 (11%), 657 (8%),
# 463 (5%), 409 (5%), 377 (4%), 363 (4%), 282 (3%)
282 / 8800

## 3) transform to long format
catchweight_common_rel_long <- as.data.frame(catchweight_common_rel)
catchweight_common_rel_long$Station <- rownames(catchweight_common_rel_long)
catchweight_common_rel_long <- melt(catchweight_common_rel_long, id.vars = "Station")

catchweight_common_rel_long <- as.data.table(catchweight_common_rel_long)[, sum(value), by = .(variable)]

catchweight_common_rel_long <- catchweight_common_rel_long %>%
  mutate(perc = case_when(
    variable %in% "Other" ~ "c) 26%",
    variable %in% "Periphylla periphylla" ~ "a) 18%",
    variable %in% "roundnose grenadier (Coryphaenoides rupestris)" ~ "b) 16%",
    variable %in% "rabbit fish (Chimaera monstrosa)" ~ "d) 11%",
    variable %in% "velvet-belly lanternshark (Etmopterus spinax)" ~ "e) 8%",
    variable %in% "blackmouth catshark (Galeus melastomus)" ~ "h) 5%",
    variable %in% "blue whiting (Micromesistius poutassou)" ~ "g) 5%",
    variable %in% "glacier lantern fish (Benthosema glaciale)" ~ "f) 4%",
    variable %in% "greater argentine (Argentina silus)" ~ "i) 4%",
    variable %in% "spiny dogfish (Squalus acanthias)" ~ "j) 3%"
  ))


pie_biomass_rel <- ggplot(catchweight_common_rel_long, aes(x = "", y = V1, fill = variable)) +
  geom_bar(stat = "identity", width = 1, alpha = 0.8) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c(
    "darkgoldenrod1", # colors in order to match prev histogram
    "darkolivegreen1",
    "darkorchid",
    "olivedrab4",
    "deeppink1",
    "royalblue3",
    "plum1",
    "darkturquoise",
    "deeppink3",
    "cyan1"
  )) +
  guides(fill = "none") +
  rremove("xylab") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    rect = element_rect(fill = "transparent")
  ) +
  labs(title = "Species % on average per catch", tag = "B") +
  theme(text = element_text(size = 14)) +
  geom_text(aes(label = perc, x = 1.37), position = position_stack(vjust = 0.5))
pie_biomass_rel


# pair plots of abundance
barplots_abundance <- ggarrange(barplot_absolute, pie_biomass_rel, ncol = 2, nrow = 1)
barplots_abundance
# ggsave(barplots_abundance, filename="barplots_abundance.png", width=13, height=5)


## 100% stacked bar plot
# change row names so that it shows in correct order in stacked bar plot
# rownames(catchweight_common_rel)[rownames(catchweight_common_rel) == "1"] <- "01"
# rownames(catchweight_common_rel)[rownames(catchweight_common_rel) == "2"] <- "02"
# rownames(catchweight_common_rel)[rownames(catchweight_common_rel) == "3"] <- "03"
# rownames(catchweight_common_rel)[rownames(catchweight_common_rel) == "4"] <- "04"
# rownames(catchweight_common_rel)[rownames(catchweight_common_rel) == "5"] <- "05"
# rownames(catchweight_common_rel)[rownames(catchweight_common_rel) == "6"] <- "06"
# rownames(catchweight_common_rel)[rownames(catchweight_common_rel) == "7"] <- "07"
# rownames(catchweight_common_rel)[rownames(catchweight_common_rel) == "8"] <- "08"
# rownames(catchweight_common_rel)[rownames(catchweight_common_rel) == "9"] <- "09"

# catchweight_stackedbar <- ggplot(catchweight_common_long,
#                                mapping = aes(x = Station, y=value, fill=variable))+
# geom_bar(position="fill", stat = "identity")+
# ylab("Relative abundance (biomass)") +
# scale_fill_manual(values = c("darkgoldenrod1",
#                            "darkolivegreen1",
#                           "darkorchid",
#                          "olivedrab4",
#                         "deeppink1",
#                       "plum1",
#                       "darkseagreen",
#                     "deeppink4",
#                     "mintcream"),name="Common species/taxa")+
#  theme(panel.border = element_blank(),
#       panel.grid.major = element_blank(),
#      panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#    axis.text = element_text(size = 11),
#   axis.title = element_text(size = 11),
#  legend.title = element_text(size = 12),
# legend.text = element_text(size = 11),
# rect = element_rect(fill = "transparent"),
# axis.text.x = element_text(size = 10, angle = 90))

# catchweight_stackedbar

# ggsave(catchweight_stackedbar, filename = "catchweight_stackedbar.png", bg = "transparent", width=400, height = 150, units = "mm", limitsize = FALSE)


### 2.3: Map of catches of common species ####
CPUE_catchweight <- CPUE_catchweight %>%
  mutate(C_rupestris_kg = `catchweight_g_Coryphaenoides rupestris` / 1000) %>%
  mutate(C_monstrosa_kg = `catchweight_g_Chimaera monstrosa` / 1000) %>%
  mutate(E_spinax_kg = `catchweight_g_Etmopterus spinax` / 1000) %>%
  mutate(G_melastomus_kg = `catchweight_g_Galeus melastomus` / 1000)


## grenadier
grenadier_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = CPUE_catchweight, color = "royalblue3", aes(
      x = longitudestart,
      y = latitudestart,
      size = CPUE_catchweight$C_rupestris_kg
    ),
    alpha = 0.6
  ) +
  labs(title = "Roundnose grenadier \n (Coryphaenoides rupestris)", tag = "A") +
  theme(text = element_text(size = 14)) +
  guides(size = guide_legend(title = "CPUE (kg/min)"))
grenadier_map

## chimaera
chimaera_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = CPUE_catchweight, color = "royalblue3", aes(
      x = longitudestart,
      y = latitudestart,
      size = CPUE_catchweight$C_monstrosa_kg
    ),
    alpha = 0.6
  ) +
  labs(title = "Rabbit fish \n(Chimaera monstrosa)", tag = "B") +
  theme(text = element_text(size = 14)) +
  guides(size = guide_legend(title = "CPUE (kg/min)"))
chimaera_map


## etmopterus
velvetbelly_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = CPUE_catchweight, color = "royalblue3", aes(
      x = longitudestart,
      y = latitudestart,
      size = CPUE_catchweight$E_spinax_kg
    ),
    alpha = 0.6
  ) +
  labs(title = "Velvet-belly lanternshark \n(Etmopterus spinax)", tag = "C") +
  theme(text = element_text(size = 14)) +
  guides(size = guide_legend(title = "CPUE (kg/min)"))
velvetbelly_map

## blackmouth catchark
galeus_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = CPUE_catchweight, color = "royalblue3",
    aes(
      x = longitudestart,
      y = latitudestart,
      size = CPUE_catchweight$G_melastomus_kg
    ),
    alpha = 0.6
  ) +
  labs(title = "Blackmouth catshark \n(Galeus melastomus)", tag = "D") +
  theme(text = element_text(size = 14)) +
  guides(size = guide_legend(title = "CPUE (kg/min)"))
galeus_map


commonfish_maps <- ggarrange(grenadier_map, chimaera_map, velvetbelly_map,
  galeus_map,
  ncol = 2, nrow = 2, common.legend = F
)
# commonfish_maps

ggsave(commonfish_maps, filename = "commonfish_maps.png", width = 10, height = 10)

### 3: Diversity ####

### 3.1: Map of stations and diversity - Figure 6 ####
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
  theme(text = element_text(size = 16)) +
  labs(title = "Diversity with Periphylla", tag = "A")
map_diversity
# ggsave(map_diversity, filename="map_diversity.png", width=10, height=8)

## map diversity without periphylla
map_diversity_minusperi <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = catch_df, aes(
      x = longitudestart, y = latitudestart,
      color = shannon_div_minusperi,
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
  theme(text = element_text(size = 16)) +
  labs(title = "Diversity without Periphylla", tag = "B")
map_diversity_minusperi

div_maps <- ggarrange(map_diversity, map_diversity_minusperi,
  ncol = 2, nrow = 1, common.legend = T,
  legend = "right"
)
div_maps


### 4: Compare trawl effect on diversity, biomass (fish and crustaceans + periphylla)####

get_box_stats <- function(y, upper_limit = max(df$mpg) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

### 4.1: Trawl effect on biomass ####
trawl_biomass <- ggplot(CPUE_catchweight, aes(
  x = Trawl, y = catchweight_tot_minusperiphylla_kg,
  fill = Trawl
)) +
  scale_fill_manual(values = alpha(c("royalblue3", "deeppink4"), .7)) +
  xlab("Trawl") +
  ylab("Fish and crustaceans catch rate (kg/min)") +
  geom_boxplot() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = -3.5) +
  labs(tag = "A") +
  theme_minimal() +
  theme(legend.position = "none")


trawl_periphylla <- ggplot(CPUE_catchweight, aes(
  x = Trawl, y = catchweight_kg_periphylla,
  fill = Trawl
)) +
  scale_fill_manual(values = alpha(c("royalblue3", "deeppink4"), .7)) +
  xlab("Trawl") +
  ylab("Periphylla catch rate (kg/min)") +
  geom_boxplot() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 1, vjust = -2) +
  labs(tag = "B") +
  theme_minimal() +
  theme(legend.position = "none")
trawl_periphylla

### 4.2: Trawl effect on diversity ####
trawl_diversity <- ggplot(env_mod, aes(
  x = Trawl, y = shannon_div,
  fill = Trawl
)) +
  scale_fill_manual(values = alpha(c("royalblue3", "deeppink4"), .7)) +
  xlab("Trawl") +
  ylab("Shannon-Wiener diversity") +
  geom_boxplot() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 1, vjust = 0.1) +
  labs(tag = "C") +
  theme_minimal() +
  theme(legend.position = "none")
trawl_diversity
# ggsave(trawl_boxplot, filename="trawl_boxplot.png", width=6)

### 4.3: Trawl effect - Figure 11 ####
trawl_catch <- ggarrange(trawl_biomass, trawl_periphylla, trawl_diversity)
trawl_catch
# ggsave(trawl_catch, filename="trawl_boxplots.png", width=8, height=5)


### 5: Div and cpue responses to sill depth, div/CPUE relationships, and TS plot ####

### 5.1: Diversity responses to sill depth ####
plot_cont_sill_div <- ggplot(env_mod, aes(x = sill_depth_m, y = shannon_div)) +
  geom_point(aes(colour = env_mod$fjord_coast)) +
  geom_smooth() +
  scale_color_manual(values = alpha(c("royalblue3", "deeppink1"))) +
  xlab("Sill depth (fjord) or bottom depth (coastal) (m)") +
  ylab("Shannon-Wiener diversity index") +
  labs(color = "Station \nlocation") +
  theme_minimal()
plot_cont_sill_div

### 5.2: Periphylla CPUE response to sill depth ####
plot_cont_sill_per <- ggplot(env_mod, aes(x = sill_depth_m, y = Periphylla_kg)) +
  geom_point(aes(colour = env_mod$fjord_coast)) +
  geom_smooth() +
  scale_color_manual(values = alpha(c("royalblue3", "deeppink1"))) +
  xlab("Sill depth (fjord) or bottom depth (coastal) (m)") +
  ylab("Periphylla CPUE (kg/min)") +
  labs(color = "Station \nlocation") +
  theme_minimal()
plot_cont_sill_per

### 5.3: Fish and crustacean CPUE response to sill depth ####
plot_cont_sill_catch <- ggplot(env_mod, aes(x = sill_depth_m, y = catchweight_minus_Periphylla)) +
  geom_point(aes(colour = env_mod$fjord_coast)) +
  geom_smooth() +
  scale_color_manual(values = alpha(c("royalblue3", "deeppink1"))) +
  xlab("Sill depth (fjord) or bottom depth (coastal) (m)") +
  ylab("Fish abd crustacean CPUE (kg/min)") +
  labs(color = "Station \nlocation") +
  theme_minimal()
plot_cont_sill_catch

### 5.4: CPUE total catch vs diversity ####
CPUE_div_plot <- ggplot(env_mod, aes(x = tot_catch, y = shannon_div)) +
  geom_point(aes(colour = env_mod$fjord_coast)) +
  geom_smooth() +
  scale_color_manual(values = alpha(c("royalblue3", "deeppink1"))) +
  xlab("CPUE (kg/min)") +
  ylab("Shannon-Wiener diversity index") +
  labs(color = "Station \nlocation") +
  theme_minimal()
CPUE_div_plot

### 5.5: CPUE fish and crustaceans vs diversity ####
CPUE_div_plot_fish <- ggplot(env_mod, aes(x = catchweight_minus_Periphylla, y = shannon_div)) +
  geom_point() +
  geom_smooth() +
  xlab("Fish and crustacean CPUE (kg/min)") +
  ylab("Shannon-Wiener diversity index") +
  theme_minimal() +
  labs(tag = "A") +
  theme(text = element_text(size = 16))
CPUE_div_plot_fish

### 5.6: CPUE Periphylla vs diversity ####
CPUE_div_plot_peri <- ggplot(env_mod, aes(x = Periphylla_kg, y = shannon_div)) +
  geom_point() +
  geom_smooth() +
  xlab("Periphylla CPUE (kg/min)") +
  ylab("Shannon-Wiener diversity index") +
  theme_minimal() +
  labs(tag = "B") +
  theme(text = element_text(size = 16))
CPUE_div_plot_peri

### 5.7: CPUE versus diversity plots together ####
CPUE_div_plots <- ggarrange(CPUE_div_plot_fish, CPUE_div_plot_peri)
CPUE_div_plots

### 5.8: TS plot for periphylla catches ####
TS_plot <- ggplot(env_mod[which(env_mod$Periphylla_kg > 0), ],
  mapping = aes(
    x = Salinity, y = Temperature,
    size = Periphylla_kg,
    label = ID
  )
) +
  geom_point() +
  geom_text_repel(aes(label = ifelse(Periphylla_kg > 20, as.character(ID), "")),
    col = "red", size = 4, nudge_x = 0.1
  ) +
  xlab("Salinity (PSU)") +
  ylab("Temperature (ºC)") +
  theme_bw() +
  guides(size = guide_legend(title = "Periphylla \n(kg/min)"))

TS_plot
# Stations with >20 kg/min CPUE Periphylla:
# Lurefjorden: ID 32 (60.68483, 5.171333) and 88 (60.82400, 5.346000)
# Mangersfjord (close to Lurefjord): 29 (60.62517, 4.960667)
# Skåneviksfjord: 20 (59.69950, 5.796333), 25 (59.77267, 5.849833)
# Matersfjord (close to Skåneviksfjord): 24 (59.83050, 5.971500)

### 6: Maps of spatial distributions of variables - Figure 4 ####

## A) Sill depth (fjord) or bottom depth (coastal), B) Station bottom depth, C) Temperature, D) ##Salinity, E) Oxygen, F) Aquaculture impact

bottomdepth_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = bottomdepth,
    shape = fjord_coast,
    size = 1
  )) +
  labs(title = "Bottom depth", tag = "A") +
  scale_color_gradientn(
    colors = c("cyan1", "darkslategray4", "royalblue4"),
    limits = c(100, 700),
    name = "(m)"
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 12)) +
  guides(size = "none", alpha = "none", shape = guide_legend("Station location"))
bottomdepth_map

env_mod$fjord_sill <- static_variables$sill_depth_m
silldepth_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = fjord_sill,
    size = 0.1
  )) +
  labs(title = "Sill depth (fjord stations)", tag = "B") +
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


temperature_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = Temperature,
    size = 0.1,
    alpha = 0.9
  )) +
  labs(title = "Temperature", tag = "C") +
  scale_color_gradientn(
    colors = c("cyan1", "slateblue", "tomato"),
    limits = c(7, 9),
    name = "(°C)"
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 12)) +
  guides(size = "none", alpha = "none")
temperature_map


salinity_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = Salinity,
    size = 0.1
  )) +
  labs(title = "Salinity", tag = "D") +
  scale_color_gradientn(
    colors = c("darkgoldenrod1", "darkslategray3", "royalblue4"),
    limits = c(33, 35.5),
    name = "(PSU)"
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 12)) +
  guides(size = "none", alpha = "none")
salinity_map

oxygen_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = Oxygen,
    size = 0.1
  )) +
  labs(title = "Oxygen", tag = "E") +
  scale_color_gradientn(
    colors = c("tomato", "orchid", "royalblue4"),
    limits = c(2, 6),
    name = "(ml/L)"
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 12)) +
  guides(size = "none", alpha = "none")
oxygen_map

aqua_map <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(data = env_mod, aes(
    x = longitudestart,
    y = latitudestart,
    color = aquaculture_impact,
    size = 0.1
  )) +
  labs(title = "Aquaculture impact score", tag = "F") +
  scale_color_gradientn(
    colors = c("darkslategray3", "goldenrod1", "tomato"),
    limits = c(0, 4100),
    name = ""
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme(text = element_text(size = 12)) +
  guides(size = "none", alpha = "none")
aqua_map


## A) Sill depth (fjord) or bottom depth (coastal), B) Station bottom depth, C) Temperature, D) Salinity, E) Oxygen, F) Aquaculture impact

env_maps1 <- ggarrange(bottomdepth_map, silldepth_map, temperature_map,
  salinity_map,
  ncol = 2, nrow = 2, common.legend = F,
  heights = c(1, 1)
)
ggsave(env_maps1, filename = "env_maps1.png", width = 9, height = 10)
env_maps2 <- ggarrange(oxygen_map, aqua_map, ncol = 2, nrow = 1, common.legend = F)
ggsave(env_maps2, filename = "env_maps2.png", width = 9, height = 5)
