## Additional supplementary checks

rm(list = ls())

## Load relevant packages and data ####
source("0_setup.R")
library(readxl)
library(dplyr)

# Rerun regression models with low salinity outliers removed ####

### Fish and Crustacean CPUE GLM  - Supplement 9.1 ####
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

env_df$sill_category <- as.factor(env_df$sill_category) # categorical variable - factor
env_df$Trawl <- as.factor(env_df$Trawl) # categorical variable - factor

# remove salinity outliers?
env_df_low_sal_removed <- filter(env_df,Salinity>34)

mod_catch_glm_log <- glm(
  log1p(catchweight_tot_minusperiphylla_kg) ~
    Oxygen
  + Temperature
  + Salinity
  + dist_coast_km
  + bottomdepth
  + (aquaculture_impact)
  + sill_category
  + Trawl,
  data = env_df_low_sal_removed
)
par(mfrow = c(2, 2))
plot(mod_catch_glm_log)

appraise(mod_catch_glm_log)
summary(mod_catch_glm_log)
vif(mod_catch_glm_log)

tiff(filename="_figures/Supplement9.1.tiff",width=4000,height=4000,
     units="px",bg="white",compression="lzw",pointsize=100)
par(mfrow = c(3, 3))
visreg(mod_catch_glm_log, "Oxygen",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Oxygen (ml/L)", ylab = "log1p(CPUE, kg/min)",
       main = "Fish and crustacean CPUE"
)

visreg(mod_catch_glm_log, "Temperature",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Temperature (ºC)", ylab = ""
)

visreg(mod_catch_glm_log, "Salinity",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Salinity (PSU) *", ylab = ""
)

visreg(mod_catch_glm_log, "bottomdepth",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Bottom depth (m)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_catch_glm_log, "dist_coast_km",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Distance to coastline (km) ***", ylab = ""
)

visreg(mod_catch_glm_log, "aquaculture_impact",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Aquaculture impact score", ylab = ""
)

visreg(mod_catch_glm_log, "Trawl",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Trawl *(F)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_catch_glm_log, "sill_category",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Sill category *(3)", ylab = ""
)
dev.off()

### Periphylla CPUE GLM - Supplement 9.3 ####
load("_data/Statistical_analysis.rda")
# uses env_mod
env_mod$sill_category <- as.factor(env_mod$sill_category) # categorical variable - factor
env_mod$Trawl <- as.factor(env_mod$Trawl) # categorical variable - factor
env_mod_low_sal_removed <- filter(env_mod,Salinity>34)

mod_peri_glm_log <- glm(
  log1p(Periphylla_kg) ~
    Oxygen
  + Temperature
  + Salinity
  + dist_coast_km
  + bottomdepth
  + (aquaculture_impact)
  + sill_category
  + Trawl,
  data = env_mod_low_sal_removed
)
par(mfrow = c(2, 2))
plot(mod_peri_glm_log)

summary(mod_peri_glm_log)
appraise(mod_peri_glm_log)
vif(mod_peri_glm_log)

tiff(filename="_figures/Supplement9.3.tiff",width=4000,height=4000,
     units="px",bg="white",compression="lzw",pointsize=100)
par(mfrow = c(3, 3))
visreg(mod_peri_glm_log, "Oxygen",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Oxygen (ml/L)", ylab = "log1p(CPUE, kg/min)",
       main = expression(paste(italic("P. periphylla"), " CPUE"))
)
visreg(mod_peri_glm_log, "Temperature",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Temperature (ºC) **", ylab = ""
)

visreg(mod_peri_glm_log, "Salinity",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Salinity (PSU)", ylab = ""
)

visreg(mod_peri_glm_log, "bottomdepth",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Bottom depth (m)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_peri_glm_log, "dist_coast_km",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Distance to coastline (km)", ylab = ""
)

visreg(mod_peri_glm_log, "aquaculture_impact",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Aquaculture impact score **", ylab = ""
)

visreg(mod_peri_glm_log, "Trawl",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Trawl (F)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_peri_glm_log, "sill_category",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Sill category **(3)", ylab = ""
)
dev.off()

### Diversity GLM - Supplement 9.5 ####
# uses env_mod_low_sal_removed
mod_diversity_glm <- glm(
  shannon_div ~
    Oxygen
  + Temperature
  + Salinity
  + dist_coast_km
  + bottomdepth
  + (aquaculture_impact)
  + sill_category
  + Trawl,
  data = env_mod_low_sal_removed
)

summary(mod_diversity_glm)
anova(mod_diversity_glm)

par(mfrow = c(2, 2))
plot(mod_diversity_glm)
vif(mod_diversity_glm)

tiff(filename="_figures/Supplement9.5.tiff",width=4000,height=4000,
     units="px",bg="white",compression="lzw",pointsize=100)
par(mfrow = c(3, 3))
visreg(mod_diversity_glm, "Oxygen",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Oxygen (ml/L)", ylab = "H' diversity", main = "Shannon-Wiener Diversity (H')"
)

visreg(mod_diversity_glm, "Temperature",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Temperature (ºC)", ylab = ""
)

visreg(mod_diversity_glm, "Salinity",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Salinity (PSU)", ylab = ""
)

visreg(mod_diversity_glm, "bottomdepth",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Bottom depth (m)  **", ylab = "H' diversity"
)

visreg(mod_diversity_glm, "dist_coast_km",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Distance to coastline (km)", ylab = ""
)

visreg(mod_diversity_glm, "aquaculture_impact",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Aquaculture impact score", ylab = ""
)

visreg(mod_diversity_glm, "Trawl",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Trawl (F)", ylab = "H' diversity"
)

visreg(mod_diversity_glm, "sill_category",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Sill category **(3)", ylab = ""
)
dev.off()

# Rerun regression models with aquaculture impact score log transformed ####

rm(list = ls())

## Load relevant packages and data ####
source("0_setup.R")
library(readxl)
library(dplyr)

### Fish and Crustacean CPUE GLM - Supplement 9.2 ####
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

env_df$sill_category <- as.factor(env_df$sill_category) # categorical variable - factor
env_df$Trawl <- as.factor(env_df$Trawl) # categorical variable - factor

# log10-transform aquaculture impact score
env_df$log_aquaculture_impact <- log10(env_df$aquaculture_impact)
hist(env_df$aquaculture_impact)
hist(env_df$log_aquaculture_impact)

mod_catch_glm_log <- glm(
  log1p(catchweight_tot_minusperiphylla_kg) ~
    Oxygen
  + Temperature
  + Salinity
  + dist_coast_km
  + bottomdepth
  + log_aquaculture_impact
  + sill_category
  + Trawl,
  data = env_df
)
par(mfrow = c(2, 2))
plot(mod_catch_glm_log)

appraise(mod_catch_glm_log)
summary(mod_catch_glm_log)
vif(mod_catch_glm_log)

tiff(filename="_figures/Supplement9.2.tiff",width=4000,height=4000,
     units="px",bg="white",compression="lzw",pointsize=100)
par(mfrow = c(3, 3))
visreg(mod_catch_glm_log, "Oxygen",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Oxygen (ml/L)", ylab = "log1p(CPUE, kg/min)",
       main = "Fish and crustacean CPUE"
)

visreg(mod_catch_glm_log, "Temperature",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Temperature (ºC)", ylab = ""
)

visreg(mod_catch_glm_log, "Salinity",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Salinity (PSU)", ylab = ""
)

visreg(mod_catch_glm_log, "bottomdepth",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Bottom depth (m)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_catch_glm_log, "dist_coast_km",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Distance to coastline (km)", ylab = ""
)

visreg(mod_catch_glm_log, "log_aquaculture_impact",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Aquaculture impact score (log10)", ylab = ""
)

visreg(mod_catch_glm_log, "Trawl",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Trawl *(F)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_catch_glm_log, "sill_category",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Sill category", ylab = ""
)
dev.off()

### Periphylla CPUE GLM - Supplement 9.4 ####
load("_data/Statistical_analysis.rda")
# uses env_mod
env_mod$sill_category <- as.factor(env_mod$sill_category) # categorical variable - factor
env_mod$Trawl <- as.factor(env_mod$Trawl) # categorical variable - factor
# log10-transform aquaculture impact score
env_mod$log_aquaculture_impact <- log10(env_mod$aquaculture_impact)
hist(env_mod$aquaculture_impact)
hist(env_mod$log_aquaculture_impact)

mod_peri_glm_log <- glm(
  log1p(Periphylla_kg) ~
    Oxygen
  + Temperature
  + Salinity
  + dist_coast_km
  + bottomdepth
  + log_aquaculture_impact
  + sill_category
  + Trawl,
  data = env_mod
)
par(mfrow = c(2, 2))
plot(mod_peri_glm_log)

summary(mod_peri_glm_log)
appraise(mod_peri_glm_log)
vif(mod_peri_glm_log)

tiff(filename="_figures/Supplement9.4.tiff",width=4000,height=4000,
     units="px",bg="white",compression="lzw",pointsize=100)
par(mfrow = c(3, 3))
visreg(mod_peri_glm_log, "Oxygen",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Oxygen (ml/L)", ylab = "log1p(CPUE, kg/min)",
       main = expression(paste(italic("P. periphylla"), " CPUE"))
)
visreg(mod_peri_glm_log, "Temperature",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Temperature (ºC) ***", ylab = ""
)

visreg(mod_peri_glm_log, "Salinity",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Salinity (PSU) *", ylab = ""
)

visreg(mod_peri_glm_log, "bottomdepth",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Bottom depth (m)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_peri_glm_log, "dist_coast_km",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Distance to coastline (km)", ylab = ""
)

visreg(mod_peri_glm_log, "log_aquaculture_impact",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Aquaculture impact score (log10)", ylab = ""
)

visreg(mod_peri_glm_log, "Trawl",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Trawl (F)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_peri_glm_log, "sill_category",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Sill category *(3)", ylab = ""
)
dev.off()

### Diversity GLM - Supplement 9.6 ####
# uses env_mod
mod_diversity_glm <- glm(
  shannon_div ~
    Oxygen
  + Temperature
  + Salinity
  + dist_coast_km
  + bottomdepth
  + log_aquaculture_impact
  + sill_category
  + Trawl,
  data = env_mod
)

summary(mod_diversity_glm)
anova(mod_diversity_glm)

par(mfrow = c(2, 2))
plot(mod_diversity_glm)
vif(mod_diversity_glm)

tiff(filename="_figures/Supplement9.6.tiff",width=4000,height=4000,
     units="px",bg="white",compression="lzw",pointsize=100)
par(mfrow = c(3, 3))
visreg(mod_diversity_glm, "Oxygen",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Oxygen (ml/L)", ylab = "H' diversity", main = "Shannon-Wiener Diversity (H')"
)

visreg(mod_diversity_glm, "Temperature",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Temperature (ºC)", ylab = ""
)

visreg(mod_diversity_glm, "Salinity",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Salinity (PSU) *", ylab = ""
)

visreg(mod_diversity_glm, "bottomdepth",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Bottom depth (m) ***", ylab = "H' diversity"
)

visreg(mod_diversity_glm, "dist_coast_km",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Distance to coastline (km) *", ylab = ""
)

visreg(mod_diversity_glm, "log_aquaculture_impact",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Aquaculture impact score (log10) *", ylab = ""
)

visreg(mod_diversity_glm, "Trawl",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Trawl (F)", ylab = "H' diversity"
)

visreg(mod_diversity_glm, "sill_category",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Sill category ***(3)", ylab = ""
)
dev.off()

# Rerun regression models with low salinity values removed and aquaculture impact score log10 transformed

rm(list = ls())

## Load relevant packages and data ####
source("0_setup.R")
library(readxl)
library(dplyr)

### Fish and Crustacean CPUE GLM - Supplement 9.7 ####
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

env_df$sill_category <- as.factor(env_df$sill_category) # categorical variable - factor
env_df$Trawl <- as.factor(env_df$Trawl) # categorical variable - factor

# log10-transform aquaculture impact score
env_df$log_aquaculture_impact <- log10(env_df$aquaculture_impact)
hist(env_df$aquaculture_impact)
hist(env_df$log_aquaculture_impact)

# remove salinity outliers?
env_df_low_sal_removed <- filter(env_df,Salinity>34)

mod_catch_glm_log <- glm(
  log1p(catchweight_tot_minusperiphylla_kg) ~
    Oxygen
  + Temperature
  + Salinity
  + dist_coast_km
  + bottomdepth
  + log_aquaculture_impact
  + sill_category
  + Trawl,
  data = env_df_low_sal_removed
)
par(mfrow = c(2, 2))
plot(mod_catch_glm_log)

appraise(mod_catch_glm_log)
summary(mod_catch_glm_log)
vif(mod_catch_glm_log)

tiff(filename="_figures/Supplement9.7.tiff",width=4000,height=4000,
     units="px",bg="white",compression="lzw",pointsize=100)
par(mfrow = c(3, 3))
visreg(mod_catch_glm_log, "Oxygen",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Oxygen (ml/L)", ylab = "log1p(CPUE, kg/min)",
       main = "Fish and crustacean CPUE"
)

visreg(mod_catch_glm_log, "Temperature",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Temperature (ºC)", ylab = ""
)

visreg(mod_catch_glm_log, "Salinity",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Salinity (PSU) *", ylab = ""
)

visreg(mod_catch_glm_log, "bottomdepth",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Bottom depth (m)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_catch_glm_log, "dist_coast_km",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Distance to coastline (km) *", ylab = ""
)

visreg(mod_catch_glm_log, "log_aquaculture_impact",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Aquaculture impact score (log10)", ylab = ""
)

visreg(mod_catch_glm_log, "Trawl",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Trawl *(F)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_catch_glm_log, "sill_category",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Sill category", ylab = ""
)
dev.off()

### Periphylla CPUE GLM - Supplement 9.8 ####
load("_data/Statistical_analysis.rda")
# uses env_mod
env_mod$sill_category <- as.factor(env_mod$sill_category) # categorical variable - factor
env_mod$Trawl <- as.factor(env_mod$Trawl) # categorical variable - factor
# log10-transform aquaculture impact score
env_mod$log_aquaculture_impact <- log10(env_mod$aquaculture_impact)
hist(env_mod$aquaculture_impact)
hist(env_mod$log_aquaculture_impact)

env_mod_low_sal_removed <- filter(env_mod,Salinity>34)

mod_peri_glm_log <- glm(
  log1p(Periphylla_kg) ~
    Oxygen
  + Temperature
  + Salinity
  + dist_coast_km
  + bottomdepth
  + log_aquaculture_impact
  + sill_category
  + Trawl,
  data = env_mod_low_sal_removed
)
par(mfrow = c(2, 2))
plot(mod_peri_glm_log)

summary(mod_peri_glm_log)
appraise(mod_peri_glm_log)
vif(mod_peri_glm_log)

tiff(filename="_figures/Supplement9.8.tiff",width=4000,height=4000,
     units="px",bg="white",compression="lzw",pointsize=100)
par(mfrow = c(3, 3))
visreg(mod_peri_glm_log, "Oxygen",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Oxygen (ml/L)", ylab = "log1p(CPUE, kg/min)",
       main = expression(paste(italic("P. periphylla"), " CPUE"))
)
visreg(mod_peri_glm_log, "Temperature",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Temperature (ºC) **", ylab = ""
)

visreg(mod_peri_glm_log, "Salinity",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Salinity (PSU)", ylab = ""
)

visreg(mod_peri_glm_log, "bottomdepth",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Bottom depth (m)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_peri_glm_log, "dist_coast_km",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Distance to coastline (km) *", ylab = ""
)

visreg(mod_peri_glm_log, "log_aquaculture_impact",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Aquaculture impact score (log10) **", ylab = ""
)

visreg(mod_peri_glm_log, "Trawl",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Trawl (F)", ylab = "log1p(CPUE, kg/min)"
)

visreg(mod_peri_glm_log, "sill_category",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Sill category **(3)", ylab = ""
)
dev.off()

### Diversity GLM - Supplement 9.9 ####
# uses env_mod_low_sal_removed
mod_diversity_glm <- glm(
  shannon_div ~
    Oxygen
  + Temperature
  + Salinity
  + dist_coast_km
  + bottomdepth
  + log_aquaculture_impact
  + sill_category
  + Trawl,
  data = env_mod_low_sal_removed
)

summary(mod_diversity_glm)
anova(mod_diversity_glm)

par(mfrow = c(2, 2))
plot(mod_diversity_glm)
vif(mod_diversity_glm)

tiff(filename="_figures/Supplement9.9.tiff",width=4000,height=4000,
     units="px",bg="white",compression="lzw",pointsize=100)
par(mfrow = c(3, 3))
visreg(mod_diversity_glm, "Oxygen",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Oxygen (ml/L)", ylab = "H' diversity", main = "Shannon-Wiener Diversity (H')"
)

visreg(mod_diversity_glm, "Temperature",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Temperature (ºC)", ylab = ""
)

visreg(mod_diversity_glm, "Salinity",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Salinity (PSU)", ylab = ""
)

visreg(mod_diversity_glm, "bottomdepth",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Bottom depth (m) **", ylab = "H' diversity"
)

visreg(mod_diversity_glm, "dist_coast_km",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Distance to coastline (km) **", ylab = ""
)

visreg(mod_diversity_glm, "log_aquaculture_impact",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Aquaculture impact score (log10) **", ylab = ""
)

visreg(mod_diversity_glm, "Trawl",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Trawl (F)", ylab = "H' diversity"
)

visreg(mod_diversity_glm, "sill_category",
       line = list(col = "grey20"), fill = list(col = "lightblue"),
       points = list(cex = .7, pch = 16), xlab = "Sill category ***(3)", ylab = ""
)
dev.off()
