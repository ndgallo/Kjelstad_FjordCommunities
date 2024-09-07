## FjordCommunities
##
## Data analysis (updated 5 September 2024)
##
#### Code sections overview ####
##  1: GLM Regression
## 1.1 Diversity GLM
## 1.2 Fish and crustacean CPUE GLM
## 1.3 Periphylla CPUE GLM
##
## 2: Cluster analysis
##  2.1 Cluster
##  3.2 IndVal and common species
##
## ### Ordination
##  4.1 Ordination
##  4.2 Fit environmental variables

## Load relevant packages and data ####
source("0_setup.R")
source("1_data_loading.R")

## for regression model, use df env_mod:
## variables: bottomdepth, sill_category, Temperature, Salinity, Oxygen,
##            dist_coast_km, aquaculture_impact, Trawl
## test against shannon_div, Periphylla_kg, catchweight_minus_periphylla

# for environmental variables in ordination, use env_df (has only variables cleaned)

###  1. Regression models - biomass and diversity ####
# use env_mod
env_mod$sill_category <- as.factor(env_mod$sill_category) # categorical variable - factor
env_mod$Trawl <- as.factor(env_mod$Trawl) # categorical variable - factor

### 1.1 Fish and Crustacean CPUE GLM - Figure 6 ####

mod_catch_glm_log <- glm(
  log1p(catchweight_minus_Periphylla) ~
    Oxygen
    + Temperature
    + Salinity
    + dist_coast_km
    + bottomdepth
    + aquaculture_impact
    + sill_category
    + Trawl,
  data = env_mod
)
par(mfrow = c(2, 2))
plot(mod_catch_glm_log)

appraise(mod_catch_glm_log)
summary(mod_catch_glm_log)

par(mfrow = c(3, 3))
visreg(mod_catch_glm_log, "Oxygen",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Oxygen (ml/L)", ylab = "log1p(CPUE, kg/min)",
  main = "A) Fish and crustacean"
)

visreg(mod_catch_glm_log, "Temperature",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Temperature (ºC)", ylab = ""
)

visreg(mod_catch_glm_log, "Salinity",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Salinity (PSU)", ylab = ""
)

visreg(mod_catch_glm_log, "bottomdepth",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Bottom depth (m)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_catch_glm_log, "dist_coast_km",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Distance to coastline (km) **", ylab = ""
)

visreg(mod_catch_glm_log, "aquaculture_impact",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Aquaculture impact score", ylab = ""
)

visreg(mod_catch_glm_log, "Trawl",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Trawl", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_catch_glm_log, "sill_category",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Sill category", ylab = ""
)

### 1.2 Periphylla CPUE GLM - Figure 7 ####
# use env_mod

mod_peri_glm_log <- glm(
  log1p(Periphylla_kg) ~
    Oxygen
    + Temperature
    + Salinity
    + dist_coast_km
    + bottomdepth
    + aquaculture_impact
    + sill_category
    + Trawl,
  data = env_mod
)
par(mfrow = c(2, 2))
plot(mod_peri_glm_log)

summary(mod_peri_glm_log)
appraise(mod_peri_glm_log)
autoplot(mod_peri_glm_log)

par(mfrow = c(3, 3))
visreg(mod_peri_glm_log, "Oxygen",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Oxygen (ml/L)", ylab = "log1p(CPUE, kg/min)",
  main = "B) Periphylla"
)

visreg(mod_peri_glm_log, "Temperature",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Temperature (ºC) ***", ylab = ""
)

visreg(mod_peri_glm_log, "Salinity",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Salinity (PSU) *", ylab = ""
)

visreg(mod_peri_glm_log, "bottomdepth",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Bottom depth (m)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_peri_glm_log, "dist_coast_km",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Distance to coastline (km)", ylab = ""
)

visreg(mod_peri_glm_log, "aquaculture_impact",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Aquaculture impact score", ylab = ""
)

visreg(mod_peri_glm_log, "Trawl",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Trawl (F)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_peri_glm_log, "sill_category",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1, pch = 1), xlab = "Sill category *(3)", ylab = ""
)


### 1.3 Diversity GLM - Figure 8 ####
# uses env_mod

# sill depth included as a continuous variable here
mod_diversity_glm <- glm(
  shannon_div ~
    Oxygen
    + Temperature
    + Salinity
    + dist_coast_km
    + bottomdepth
    + aquaculture_impact
    + sill_category
    + Trawl,
  data = env_mod
)

summary(mod_diversity_glm)
anova(mod_diversity_glm)

par(mfrow = c(2, 2))
plot(mod_diversity_glm)

par(mfrow = c(3, 3))
visreg(mod_diversity_glm, "Oxygen",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1.5, pch = 1), xlab = "Oxygen (ml/L)", ylab = "H' diversity"
)

visreg(mod_diversity_glm, "Temperature",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1.5, pch = 1), xlab = "Temperature (ºC)", ylab = ""
)

visreg(mod_diversity_glm, "Salinity",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1.5, pch = 1), xlab = "Salinity (PSU)", ylab = ""
)

visreg(mod_diversity_glm, "bottomdepth",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1.5, pch = 1), xlab = "Bottom depth (m)  ***", ylab = "H' diversity"
)

visreg(mod_diversity_glm, "dist_coast_km",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1.5, pch = 1), xlab = "Distance to coastline (km)", ylab = ""
)

visreg(mod_diversity_glm, "aquaculture_impact",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1.5, pch = 1), xlab = "Aquaculture impact score", ylab = ""
)

visreg(mod_diversity_glm, "Trawl",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1.5, pch = 1), xlab = "Trawl (F)", ylab = "H' diversity"
)

visreg(mod_diversity_glm, "sill_category",
  line = list(col = "grey20"), fill = list(col = "lightblue"),
  points = list(cex = 1.5, pch = 1), xlab = "Sill category **(3)", ylab = ""
)

### 2. Cluster analysis and IndVal - Figure 9 ####

### 2.1 Cluster analysis ####
### Cluster -  Finding groups of similar observations

## hierarchical clustering
# calculate distance matrix
# use species_matrix_sqrt
species_vegdist_bray <- vegdist(species_matrix_sqrt, method = "bray")

# ward's method
species_hclust_ward <- hclust(species_vegdist_bray, method = "ward.D2")

## calculate how many clusters
# within-sum squares data - elbow method
fviz_nbclust(species_matrix_sqrt, kmeans, method = "wss") +
  labs(subtitle = "Elbow method")
# --> 4 clusters

library(dendextend)
dend_obj <- as.dendrogram(species_hclust_ward)
plot(dend_obj)
col_dend <- color_branches(dend_obj,
  h = 2,
  col = c("deeppink1", "seagreen", "darkgoldenrod1", "royalblue3"),
  groupLabels = TRUE
)
par(mfrow = c(1, 1))
dend_clust <- plot(col_dend)

cluster_stations <- data.frame(cutree(species_hclust_ward, k = 4))
colnames(cluster_stations)[1] <- "cluster"
cluster_stations <- cbind(ID = 1:nrow(cluster_stations), cluster_stations)

# numbers don't match:
# 1 in df is actually cluster 2
# 2 in df is actually cluster 1
# 3 in df is actually cluster 4
# 4 in df is actually cluster 3

cluster_stations$cluster[cluster_stations$cluster == "1"] <- "B"
cluster_stations$cluster[cluster_stations$cluster == "2"] <- "A"
cluster_stations$cluster[cluster_stations$cluster == "3"] <- "D"
cluster_stations$cluster[cluster_stations$cluster == "4"] <- "C"
cluster_stations$cluster[cluster_stations$cluster == "A"] <- "1"
cluster_stations$cluster[cluster_stations$cluster == "B"] <- "2"
cluster_stations$cluster[cluster_stations$cluster == "C"] <- "3"
cluster_stations$cluster[cluster_stations$cluster == "D"] <- "4"

# now it should match
# add to catch df
catch_df <- merge(catch_df, cluster_stations, by = "ID")

# map
map_cluster <- basemap(limits = c(4.5, 7.5, 59.4, 62.1), land.col = "gray90") +
  geom_spatial_point(
    data = catch_df, aes(
      x = longitudestart, y = latitudestart,
      color = factor(cluster)
    ),
    size = 4, alpha = 0.7
  ) +
  scale_color_manual(
    values = c("deeppink1", "seagreen", "darkgoldenrod1", "royalblue3"),
    name = "Cluster"
  ) +
  xlab("Longitude") + ylab("Latitude") +
  labs(title = "", tag = "") +
  theme(text = element_text(size = 14), legend.position = "none") +
  guides(color = guide_legend(override.aes = list(size = 4)))
map_cluster

# ggsave(map_cluster, filename="map_cluster.png", width=5, height=8)

### 2.2 IndVal analysis ####
# install.packages("labdsv")
library(labdsv)

# uses CPUE_catchweight
sp_iva <- CPUE_catchweight %>%
  select(starts_with("catchweight")) %>%
  subset(select = -c(catchweight_total, catchweight_tot_minusperiphylla_kg, catchweight_kg_periphylla))
sp_iva <- sp_iva[, (!apply(sp_iva == 0, 2, all))]

iva <- indval(sp_iva, cluster_stations$cluster)

summary(iva)

# most abundant species
catchweight_relative$cluster <- catch_df$cluster

catch_clust_1 <- catchweight_relative[catchweight_relative$cluster == "1", ] %>%
  subset(select = -cluster)
catch_clust_2 <- catchweight_relative[catchweight_relative$cluster == "2", ] %>%
  subset(select = -cluster)
catch_clust_3 <- catchweight_relative[catchweight_relative$cluster == "3", ] %>%
  subset(select = -cluster)
catch_clust_4 <- catchweight_relative[catchweight_relative$cluster == "4", ] %>%
  subset(select = -cluster)

str(sort(colSums(catch_clust_1[, 1:length(catch_clust_1)]), decreasing = TRUE)[1:3])
str(sort(colSums(catch_clust_1[, 1:length(catch_clust_1)]), decreasing = TRUE)[4:5])
# cluster 1: Chimaera monstrosa (596/2700 = 22%), Argentina silus(284 = 11%),
# Etmopterus spinax(266 = 10%), Galeus mealstomus (248 = 9%), Squalus acanthias (233 = 9%)
233 / 2700

str(sort(colSums(catch_clust_2[, 1:length(catch_clust_2)]), decreasing = TRUE)[1:3])
str(sort(colSums(catch_clust_2[, 1:length(catch_clust_2)]), decreasing = TRUE)[4:5])
# cluster 2: Etmopterus spinax (291/2300 = 13%), Pollachius virens (273 = 12%),
# BEnthosema glaciale (217 = 9%), Pollachius pollachius (167 = 7%), Periphylla (162 = 7%)
162 / 2300

str(sort(colSums(catch_clust_3[, 1:length(catch_clust_3)]), decreasing = TRUE)[1:3])
str(sort(colSums(catch_clust_3[, 1:length(catch_clust_3)]), decreasing = TRUE)[4:5])
# cluster 1: Coryphaenoides rupestris (1284/2200 = 58%), Chimaera monstrosa (221 = 10%),
# Galeus melastomus (186 = 8%), Periphylla (114 = 5%), Molva dypterygia (79 = 4%)
79 / 2200

str(sort(colSums(catch_clust_4[, 1:length(catch_clust_4)]), decreasing = TRUE)[1:3])
str(sort(colSums(catch_clust_4[, 1:length(catch_clust_4)]), decreasing = TRUE)[4:5])
# cluster 1: Periphylla (1200 = 75%), Chimaera monstrosa (74 = 5%),
# Micromesistius poutassou (73 = 5%), Coryphaenoides rupestris(52 = 3%), B glac (49 = 3%)
49 / 1600

mean(catch_df[catch_df$cluster == "1", "total_catch_CPUE"])
mean(catch_df[catch_df$cluster == "2", "total_catch_CPUE"])
mean(catch_df[catch_df$cluster == "3", "total_catch_CPUE"])
mean(catch_df[catch_df$cluster == "4", "total_catch_CPUE"])

range(catch_df[catch_df$cluster == "1", "bottomdepth"])
range(catch_df[catch_df$cluster == "2", "bottomdepth"])
range(catch_df[catch_df$cluster == "3", "bottomdepth"])
range(catch_df[catch_df$cluster == "4", "bottomdepth"])
mean(catch_df[catch_df$cluster == "1", "bottomdepth"])
mean(catch_df[catch_df$cluster == "2", "bottomdepth"])
mean(catch_df[catch_df$cluster == "3", "bottomdepth"])
mean(catch_df[catch_df$cluster == "4", "bottomdepth"])
median(catch_df[catch_df$cluster == "1", "bottomdepth"])
median(catch_df[catch_df$cluster == "2", "bottomdepth"])
median(catch_df[catch_df$cluster == "3", "bottomdepth"])
median(catch_df[catch_df$cluster == "4", "bottomdepth"])

#### 3. Ordination ####

# unimodal or linear distribution?
# uses species_matrix_sqrt
catch_DCA <- decorana(species_matrix_sqrt)
catch_DCA
# Axis lengths 3.2416 2.6714 2.3497 1.7177 -> Unimodal as axis lengths are > 3 (linear < 3)

## Since it is unimodal -> Correspondence Analysis CA & DCA
## --> if artifact (arch) in CA, use DCA. (if linear -> PCA)

## Choice or no underlying response model
##  -> Principal Co-ordinates Analysis PCoA
##  -> Non-metric Multidimensional Scaling NMDS

catch_CA <- cca(species_matrix_sqrt)
par(mfrow = c(2, 2))
plot(catch_CA)
# -> there is an arch/diamond shape -> artifact in the CA plot -> do not use

# DCA
plot(catch_DCA)

# NMDS
set.seed(10)
catch_NMDS <- metaMDS(species_matrix_sqrt, distance = "bray")
catch_NMDS
plot(catch_NMDS)
# stress = 0.2028 2D
stressplot(catch_NMDS)

catch_NMDS$points
catch_NMDS$species

nmds_tibble <- as_tibble(catch_NMDS$points) %>%
  mutate(cluster = catch_df$cluster)
nmds_tibble$station <- rownames(nmds_tibble)

species.scores <- as.data.frame(vegan::scores(catch_NMDS, "species"))
species.scores$species <- rownames(species.scores)

selected_sp <- species.scores[species.scores$species %in% c(
  "catchweight_g_Periphylla",
  "catchweight_g_Coryphaenoides rupestris",
  "catchweight_g_Caridea",
  "catchweight_g_Etmopterus spinax",
  "catchweight_g_Chimaera monstrosa",
  "catchweight_g_Galeus melastomus",
  "catchweight_g_Pollachius virens",
  "catchweight_g_Jellyfish spp",
  "catchweight_g_Argentina silus",
  "catchweight_g_Micromesistius poutassou",
  "catchweight_g_Benthosema glaciale",
  "catchweight_g_Maurolicus muelleri",
  "catchweight_g_Squalus acanthias",
  "catchweight_g_Galatheidae",
  "catchweight_g_Myxine glutinosa",
  "catchweight_g_Dendrobranchiata",
  "catchweight_g_Molva dypterygia",
  "catchweight_g_Cephalopoda"
), ]

head(selected_sp)
selected_sp$species[selected_sp$species == "catchweight_g_Periphylla"] <- "Periphylla"
selected_sp$species[selected_sp$species == "catchweight_g_Coryphaenoides rupestris"] <- "C.rupestris"
selected_sp$species[selected_sp$species == "catchweight_g_Caridea"] <- "Shrimp spp."
selected_sp$species[selected_sp$species == "catchweight_g_Etmopterus spinax"] <- "E.spinax"
selected_sp$species[selected_sp$species == "catchweight_g_Chimaera monstrosa"] <- "C.monstrosa"
selected_sp$species[selected_sp$species == "catchweight_g_Galeus melastomus"] <- "G.melastomus"
selected_sp$species[selected_sp$species == "catchweight_g_Pollachius virens"] <- "P.virens"
selected_sp$species[selected_sp$species == "catchweight_g_Jellyfish spp"] <- "Jellyfish spp"
selected_sp$species[selected_sp$species == "catchweight_g_Argentina silus"] <- "A.silus"
selected_sp$species[selected_sp$species == "catchweight_g_Micromesistius poutassou"] <- "M.poutassou"
selected_sp$species[selected_sp$species == "catchweight_g_Benthosema glaciale"] <- "B.glaciale"
selected_sp$species[selected_sp$species == "catchweight_g_Maurolicus muelleri"] <- "M.muelleri"
selected_sp$species[selected_sp$species == "catchweight_g_Squalus acanthias"] <- "S.acanthias"
selected_sp$species[selected_sp$species == "catchweight_g_Galatheidae"] <- "Galatheidae"
selected_sp$species[selected_sp$species == "catchweight_g_Myxine glutinosa"] <- "M.glutinosa"
selected_sp$species[selected_sp$species == "catchweight_g_Dendrobranchiata"] <- "pel. shrimp"
selected_sp$species[selected_sp$species == "catchweight_g_Molva dypterygia"] <- "M.dypterygia"
selected_sp$species[selected_sp$species == "catchweight_g_Cephalopoda"] <- "Cephalopoda"

nmds_tibble$cluster <- as.factor(nmds_tibble$cluster)

nmds_plot <- ggplot() +
  geom_point(
    data = nmds_tibble, aes(x = MDS1, y = MDS2, color = cluster, shape = cluster, fill = cluster),
    size = 4, alpha = 0.5
  ) +
  coord_fixed() +
  scale_shape_manual(values = c(15, 16, 17, 23)) +
  scale_color_manual(values = c("deeppink1", "seagreen", "darkgoldenrod1", "royalblue3")) +
  scale_fill_manual(values = c("deeppink1", "seagreen", "darkgoldenrod1", "royalblue3")) +
  ggrepel::geom_text_repel(
    data = selected_sp, aes(x = NMDS1, y = NMDS2, label = species),
    alpha = 0.8, size = 3.5, force = 0.02
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

nmds_plot

# ggsave(nmds_plot, filename="nmds_plot.png", width=6.5, height=6)

## Direct gradient analysis (not used)
# Detrended (CCA) = constrained DCA
par(mfrow = c(1, 1))
catch_RDA <- rda(species_matrix_sqrt ~ ., data = env_df)
catch_RDA
plot(catch_RDA)

### fit environmental variables to ordination
# uses catch_NMDS and env_fit
catch_envfit <- envfit(catch_NMDS, env_df,
  choices = 1:2,
  permutations = 999, na.rm = TRUE
)
catch_envfit

plot(catch_NMDS)
plot(catch_envfit, p.max = .05)
# dist_coast *, bottom depth ***, oxygen **, temp ***, salinity ***

# add env to NMDS ggplot
env.scores <- as.data.frame(vegan::scores(catch_envfit, display = "vectors"))
env.scores <- cbind(env.scores, Vectors = rownames(env.scores))
env.scores

env.scores$Vectors[env.scores$Vectors == "dist_coast_km"] <- "Coast.dist"
env.scores$Vectors[env.scores$Vectors == "aquaculture_impact"] <- "Aqua.imp"
env.scores$Vectors[env.scores$Vectors == "bottomdepth"] <- "Depth"
env.scores$Vectors[env.scores$Vectors == "latitudestart"] <- "Latitude"
env.scores$Vectors[env.scores$Vectors == "sill_depth_m"] <- "Sill.depth"

env.scores.fac <- as.data.frame(vegan::scores(catch_envfit, display = "factors"))
env.scores.fac <- cbind(env.scores.fac, factors.env = rownames(env.scores.fac))
env.scores.fac

env.scores.fac$factors.env[env.scores.fac$factors.env == "TrawlC"] <- "Trawl.C"
env.scores.fac$factors.env[env.scores.fac$factors.env == "TrawlF"] <- "Trawl.F"

nmds_plot_env <- ggplot() +
  geom_point(
    data = nmds_tibble, aes(x = MDS1, y = MDS2, color = cluster, shape = cluster, fill = cluster),
    size = 4, alpha = 0.5
  ) +
  coord_fixed() +
  scale_shape_manual(values = c(15, 16, 17, 23)) +
  scale_color_manual(values = c("deeppink1", "seagreen", "darkgoldenrod1", "royalblue3")) +
  scale_fill_manual(values = c("deeppink1", "seagreen", "darkgoldenrod1", "royalblue3")) +
  geom_segment(
    data = env.scores,
    aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
    arrow = arrow(length = unit(0.25, "cm")), colour = "grey20"
  ) +
  geom_text(
    data = env.scores, aes(x = NMDS1, y = NMDS2, label = Vectors),
    size = 3.5
  ) +
  geom_point(
    data = env.scores.fac, aes(x = NMDS1, y = NMDS2),
    shape = "+", size = 6, alpha = 0.6, colour = "brown4"
  ) +
  ggrepel::geom_text_repel(
    data = env.scores.fac, aes(x = NMDS1, y = NMDS2 + 0.05),
    label = env.scores.fac$factors.env, colour = "grey20", size = 3, force = 0.001
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())
nmds_plot_env

nmds_plot_env_cont <- ggplot() +
  geom_point(
    data = nmds_tibble, aes(x = MDS1, y = MDS2, color = cluster),
    size = 4, alpha = 0.5
  ) +
  coord_fixed() +
  scale_color_manual(values = c("deeppink1", "seagreen", "darkgoldenrod1", "royalblue3")) +
  geom_segment(
    data = env.scores,
    aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
    arrow = arrow(length = unit(0.25, "cm")), colour = "grey20"
  ) +
  geom_text(
    data = env.scores, aes(x = NMDS1, y = NMDS2, label = Vectors),
    size = 3.5
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())
nmds_plot_env_cont

# ggsave(nmds_plot_env, filename="nmds_plot_env.png", width=6.5, height=6)

## ordisurf - non-linear relationships between community data and env var

surf_ox <- ordisurf(catch_NMDS ~ env_df$Oxygen, knots = 10, isotropic = TRUE, main = "", )
surf_ox
summary(surf_ox)
# ***, dev ex = 33.3%

surf_depth <- ordisurf(catch_NMDS ~ env_df$bottomdepth, knots = 10, main = "")
summary(surf_depth)
# ***, dev ex = 75%

surf_sal <- ordisurf(catch_NMDS ~ env_df$Salinity, knots = 10, main = "")
summary(surf_sal)
# ***, dev ex = 63.4%

surf_temp <- ordisurf(catch_NMDS ~ env_df$Temperature, knots = 10, main = "")
summary(surf_temp)
# ***, dev ex = 49.5%

surf_aq <- ordisurf(catch_NMDS ~ env_df$aquaculture_impact, knots = 10, main = "")
summary(surf_aq)
# dev ex = 0.0%

surf_dist <- ordisurf(catch_NMDS ~ env_df$dist_coast_km, knots = 10, main = "")
summary(surf_dist)
## **, dev ex = 10.2%
head(surf_dist)
