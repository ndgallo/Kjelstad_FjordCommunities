## Data loading for FjordCommunities project
##
##
## Sep 2023
##
## Code sections overview ####
## 1. IMR catch and CTD
## 1.1: load IMR biotic data
## 1.2: load IMR CTD data
## 1.3: create biotic df
## 1.4: create CTD dfs w bottom depth variables
## 1.5: combine biotic + ctd + clean out benthos
##
## 2. IMR catch and CTD
## 2.1: load UiB catch data
## 2.2: load UiB CTD data, added bottom depth variables to catch dfs 2011-2018
## 2.3: merge UiB data and remove benthos
## 2.4: rename all columns to match IMR data
##
## 3. Merge UiB and IMR catch data and standardize
## 3.1: merge IMR and UiB
## 3.2: standardize catch data CPUE + make relative df
## 3.3: make combined dataframe with all variables
##

## Load relevant packages (source script installs missing packages) ####
source("0_setup.R")

# other packages:
library(cowplot)
library(ggpubr)
library(plyr)
library(janitor)
library(vegan)
library(cluster)
library(densityClust)
library(ggplot2)
library(GGally)
library(dplyr)
library(ggrepel)
library(factoextra)
library(mgcv)
library(gratia)
library(visreg)
library(ggordiplots)

## Prepare datasets ####

#### 1: IMR catch and CTD ####

###### 1.1: IMR biotic data ####

# load files
biotic_files <- list.files("./_data/biotic", pattern = ".xml", full.names = TRUE)
biotic <- ReadBiotic(biotic_files) # load files

# get out relevant tables
for (i in seq_along(biotic)) {
  mission_df <- rbindlist(lapply(biotic, "[[", "mission"), fill = TRUE) %>% as_tibble()
  station_df <- rbindlist(lapply(biotic, "[[", "fishstation"), fill = TRUE) %>% as_tibble()
  catch_df <- rbindlist(lapply(biotic, "[[", "catchsample"), fill = TRUE) %>% as_tibble()
  individ_df <- rbindlist(lapply(biotic, "[[", "individual"), fill = TRUE) %>% as_tibble()
}

# check contents
table(station_df$startyear)
table(station_df$gear, exclude = FALSE)
# 2154 3296 3513 3548 4401 <NA>
# 11   18   11    7   13   35
# important to use only the stations with the right gear types
# 2154 = mammoth 32xx = bottom trawls 35xx = pelagic trawls 4401 = dredge
# NB: the shrimp trawl from 2021 is currently NA, we need to fix that

# combining data: for the analysis the relevant tables need to be merged (e.g. station and catch)
# to do so each table includes common identifiers - serialnumber is a unique station code that in all
# tables from station downward

###### 1.2: load IMR CTD data ####

# NB: because we used two different CTD systems on the two IMR cruises,
# loading and combining the data is different for each

# Filter out only bottom values

# Survey 2021854
# in this survey we have one CTD for each trawl station, so we can easily map them together sequentially
files <- list.files(path = "./_data/ctd/2021854", pattern = ".xlsx")


ctd_df1 <- data.frame()
for (i in seq_along(files)) {
  df <- readxl::read_xlsx(path = paste0("./_data/ctd/2021854/", files[i]), sheet = "Data", skip = 1)
  df$serialnumber <- unique(subset(station_df, missionnumber == 54)$serialnumber)[i]
  ctd_df1 <- rbind(ctd_df1, df)
}
names(ctd_df1)[c(3, 6, 7, 13)] <- c("Temperature", "Temperature2", "O2concentration", "O2saturation")

ctd_df1 <- ctd_df1 %>% filter(Depth > 2)


# Need to convert oxygen from micromol/L to ml/L
# using molar volume at standard pressure and temperature
# 1mol=22.391L

ctd_df1 <- ctd_df1 %>%
  mutate(oxygen = ((O2concentration * 22.391) / 1000))


# Survey 2022603
# in this survey we had separate CTD stations, so we need to map them together by spatial proximity

# read.cnv
# install package to use read.ctd
# install.packages("oce")
library(oce)

files <- list.files(path = "./_data/ctd/2022603", pattern = ".cnv")
ctd_df2 <- data.table()
for (i in seq_along(files)) {
  df <- read.ctd(paste0("./_data/ctd/2022603/", files[i]))
  df <- ctdTrim(df, method = "downward")
  df <- df@data %>%
    as.data.table() %>%
    mutate(
      lon = df@metadata$longitude, lat = df@metadata$latitude,
      date = df@metadata$startTime, ctdstation = df@metadata$station,
      ship = df@metadata$ship, filename = files[i]
    )
  ctd_df2 <- rbind(df, ctd_df2)
}

ctd_df2$Depth <- ctd_df2$pressure

# code snipped to find closest ctd-trawl combination
zone <- floor((mean(subset(station_df, missionnumber == 3)$longitudestart, na.rm = T) + 180) / 6) + 1
utmCRS <- CRS(paste0("+proj=utm +zone=", zone, " +datum=WGS84 +units=km +no_defs"))

stationSP <- SpatialPoints(
  subset(station_df, missionnumber == 3) %>% select(longitudestart, latitudestart),
  CRS("+proj=longlat +datum=WGS84")
)
stationSP <- spTransform(stationSP, utmCRS)

ctdSP <- SpatialPoints(
  ctd_df2 %>% distinct(lon, lat),
  CRS("+proj=longlat +datum=WGS84")
)
ctdSP <- spTransform(ctdSP, utmCRS)


dist <- gDistance(ctdSP, stationSP, byid = T) # distances between stations
minDist <- apply(dist, 1, function(x) order(x, decreasing = F)[1]) # find closest
ctd_match <- data.frame(
  serialnumber = unique(subset(station_df, missionnumber == 3)$serialnumber),
  ctdstation = unique(ctd_df2$ctdstation)[minDist],
  distance = dist[minDist]
)

ctd_df2 <- inner_join(ctd_df2, ctd_match, by = "ctdstation")

#plots to look at CTD and trawl locations
#p1 <- basemap(limits=c(4.5,7.5,59.5,62.3),land.col="lightgrey") +
#geom_spatial_text(data=ctd_df2,aes(x=lon,y=lat,label=ctdstation)) +
#xlab("Lengdegrad") + ylab("Breddegrad")

#p2 <- basemap(limits=c(4.5,7.5,59.5,62.3),land.col="lightgrey") +
#geom_spatial_text(data=subset(station_df,missionnumber==3),aes(x=longitudestart,y=latitudestart,label=serialnumber))+
#xlab("Lengdegrad") + ylab("Breddegrad")

#ggpubr::ggarrange(p1,p2)

###### 1.3: Creating new biotic_df with station and catch dfs ####

# I change 2021 NA gear type to 3200

station_df <- station_df %>% mutate(gear = ifelse(is.na(gear), 3200, gear))
table(station_df$gear, exclude = FALSE)

station_df <- station_df %>%
  filter(gear %in% c("3296", "3200")) # filtered out for gear type


# combine catch and station dfs by serial number
# keeping variables I might need or find interesting
biotic_df <- inner_join(catch_df, station_df, by = c("startyear", "serialnumber", "platform", "missionnumber")) %>% # FZ: good to have all common columns in there as fail safe
  mutate(
    bottomdepth = (bottomdepthstop + bottomdepthstart) / 2,
    fishingtime_min = soaktime * 60
  )

# there's multiple stations with >1 subsamples of the same species. This will cause problems down the road, good to clean it up by summing to total catch per station
detach(package:plyr) # had to detach plyr (to be able to use dplyr..) for grouped summarise
biotic_df2 <- biotic_df %>%
  group_by(
    serialnumber, scientificname, commonname, startyear,
    agingstructure, stationstartdate, stationstarttime,
    stationstoptime, latitudestart, longitudestart, latitudeend,
    longitudeend, fishingdepthmax, bottomdepth,
    fishingdepthmin, vesselspeed, distance, trawldoorspread,
    wirelength, fishingtime_min
  ) %>%
  summarise(
    catchcount = sum(catchcount, na.rm = T),
    catchweight = sum(catchweight, na.rm = T)
  ) %>%
  ungroup()
library(plyr) # reload plyr

biotic_df3 <- biotic_df2 %>%
  mutate(catchweight_g = catchweight * 1000) %>%
  dplyr::select(serialnumber, scientificname, catchcount, catchweight_g, latitudestart, longitudestart, bottomdepth, fishingtime_min, startyear, distance)


### change to wide format
biotic_df_wide <- pivot_wider(biotic_df3,
  names_from = scientificname,
  values_from = c(catchcount, catchweight_g),
  values_fill = 0
)
# add month to df
month <- c(
  2, 2, 2, 2, 2, 2, 2, 2, 2,
  2, 2, 2, 2, 2, 2, 2, 2,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11, 11, 11
)
biotic_df_wide$month <- month

###### 1.4: Bottom CTD data IMR ####

# Find CTD depths equal to bottom depths
bottomdepths <- biotic_df_wide %>%
  select(serialnumber, bottomdepth)

bottomdepths$bottomdepth <- round(bottomdepths$bottomdepth, 0) # removing decimals

# 2021
bottomdepths_2021 <- bottomdepths[19:53, ] # choose only from 2021 survey

ctd_df1$Depth <- round(ctd_df1$Depth, 0) # removing decimals


ctd_bottom_df1 <- ctd_df1 %>%
  group_by(serialnumber) %>%
  dplyr::summarise(
    Bottom_depth_CTD = max(Depth),
    Temperature = mean(Temperature[Depth > 0.99 * max(Depth)]),
    Salinity = mean(Salinity[Depth > 0.99 * max(Depth)]),
    Oxygen = mean(oxygen[Depth > 0.99 * max(Depth)])
  )

ctd_bottom_df1 <- inner_join(ctd_bottom_df1, bottomdepths_2021, by = "serialnumber")
# depths match up pretty good


# 2022
bottomdepths_2022 <- bottomdepths[1:18, ] # choose only from 2022 survey

ctd_bottom_df2 <- ctd_df2 %>%
  group_by(serialnumber) %>%
  dplyr::summarise(
    Bottom_depth_CTD = max(Depth),
    Temperature = mean(temperature[Depth > 0.99 * max(Depth)]),
    Salinity = mean(salinity[Depth > 0.99 * max(Depth)]),
    Oxygen = mean(oxygen[Depth > 0.99 * max(Depth)])
  )

ctd_bottom_df2 <- inner_join(ctd_bottom_df2, bottomdepths_2022, by = "serialnumber")


# combine CTD dfs
ctd_bottom_IMR <- rbind(ctd_bottom_df1, ctd_bottom_df2) %>%
  select(-bottomdepth)

###### 1.5: combine CTD and biotic + clean ####

# Add CTD to biotic df
biotic_df_wide <- left_join(biotic_df_wide, ctd_bottom_IMR, by = "serialnumber")

# find empty columns and remove them
which(colSums(biotic_df_wide == 0) == nrow(biotic_df_wide))

IMR_catch_clean <- biotic_df_wide %>% # remove empty columns
  subset(select = -c(
    catchcount_Lophogaster, catchcount_Sergestes,
    catchcount_Echinoidea, catchcount_Fucales,
    catchcount_Polychaeta, catchcount_Priapulus,
    `catchcount_Psilaster andromeda`, catchcount_Alcyonacea,
    catchcount_NA, catchcount_Laminariaceae,
    `catchcount_Echinocardium cordatum`, catchcount_Pandalina,
    catchcount_Phaeophyceae, catchcount_Crustacea
  ))

### remove benthos
IMR_catch_clean <- IMR_catch_clean %>% # remove empty columns
  subset(select = -c(
    `catchcount_Parastichopus Tremulus`,
    `catchweight_g_Parastichopus Tremulus`,
    catchcount_Ophiuroidea, catchweight_g_Ophiuroidea,
    catchcount_Porifera, catchweight_g_Porifera,
    catchcount_Echinidae, catchweight_g_Echinidae,
    catchcount_Hydroidolina, catchweight_g_Hydroidolina,
    catchcount_Onchidorididae, catchweight_g_Onchidorididae,
    catchcount_Phakellia, catchweight_g_Phakellia,
    catchcount_Echinidea, catchweight_g_Echinidea,
    catchcount_Pennatulacea, catchweight_g_Pennatulacea,
    catchcount_Hydrozoa, catchweight_g_Hydrozoa,
    catchcount_Asteronyx, catchweight_g_Asteronyx,
    catchcount_Bivalvia, catchweight_g_Bivalvia,
    catchweight_g_Fucales, catchweight_g_Polychaeta,
    catchweight_g_Priapulus, `catchweight_g_Psilaster andromeda`,
    catchcount_Psilaster, catchweight_g_Psilaster,
    `catchcount_Scaphander lignarius`, `catchweight_g_Scaphander lignarius`,
    catchcount_Spatangoida, catchweight_g_Spatangoida,
    catchweight_g_Alcyonacea, catchweight_g_Laminariaceae,
    catchcount_Anthozoa, catchweight_g_Anthozoa,
    `catchweight_g_Echinocardium cordatum`,
    catchcount_Holothuroidea, catchweight_g_Holothuroidea,
    catchcount_Sepietta, catchweight_g_Sepietta,
    `catchcount_Sepietta Neglecta`, `catchweight_g_Sepietta Neglecta`,
    catchweight_g_Phaeophyceae,
    catchcount_Actiniaria, catchweight_g_Actiniaria,
    catchcount_Asteroidea, catchweight_g_Asteroidea,
    catchweight_g_NA, catchweight_g_Echinoidea
  ))


# pool invertebrates
IMR_catch_clean <- IMR_catch_clean %>%
  mutate(
    "catchcount_Cephalopoda" = (IMR_catch_clean$catchcount_Cephalopoda
      + IMR_catch_clean$`catchcount_Todaropsis eblanae`
      + IMR_catch_clean$`catchcount_Rossia macrosoma`),
    "catchweight_g_Cephalopoda" = (IMR_catch_clean$catchweight_g_Cephalopoda
      + IMR_catch_clean$`catchweight_g_Todaropsis eblanae`
      + IMR_catch_clean$`catchweight_g_Rossia macrosoma`),
    "catchcount_Caridea" = (IMR_catch_clean$`catchcount_Atlantopandalus propinqvus`
      + IMR_catch_clean$`catchcount_Dichelopandalus bonnieri`
      + IMR_catch_clean$`catchcount_Pandalus borealis`
      + IMR_catch_clean$catchcount_Pontophilus),
    "catchweight_g_Caridea" = (IMR_catch_clean$`catchweight_g_Atlantopandalus propinqvus`
      + IMR_catch_clean$`catchweight_g_Dichelopandalus bonnieri`
      + IMR_catch_clean$`catchweight_g_Pandalus borealis`
      + IMR_catch_clean$catchweight_g_Pandalina
      + IMR_catch_clean$catchweight_g_Pontophilus),
    "catchcount_Dendrobranchiata" = (IMR_catch_clean$catchcount_Pasiphaea
      + IMR_catch_clean$`catchcount_Pasiphaea multidentata`
      + IMR_catch_clean$`catchcount_Pasiphaea sivado`
      + IMR_catch_clean$`catchcount_Pasiphaea tarda`
      + IMR_catch_clean$catchcount_Pasiphaeidae
      + IMR_catch_clean$`catchcount_Sergestes arcticus`),
    "catchweight_g_Dendrobranchiata" = (IMR_catch_clean$catchweight_g_Pasiphaea
      + IMR_catch_clean$`catchweight_g_Pasiphaea multidentata`
      + IMR_catch_clean$`catchweight_g_Pasiphaea tarda`
      + IMR_catch_clean$`catchweight_g_Pasiphaea sivado`
      + IMR_catch_clean$catchweight_g_Pasiphaeidae
      + IMR_catch_clean$catchweight_g_Sergestes
      + IMR_catch_clean$`catchweight_g_Sergestes arcticus`),
    "catchcount_Galatheidae" = (IMR_catch_clean$catchcount_Munida
      + IMR_catch_clean$`catchcount_Munida rugosa`),
    "catchweight_g_Galatheidae" = (IMR_catch_clean$catchweight_g_Munida
      + IMR_catch_clean$`catchweight_g_Munida rugosa`),
    "catchcount_Isopoda" = (IMR_catch_clean$catchcount_Isopoda
      + IMR_catch_clean$catchcount_Aegidae),
    "catchweight_g_Isopoda" = (IMR_catch_clean$catchweight_g_Isopoda
      + IMR_catch_clean$catchweight_g_Aegidae),
    "catchcount_Benthosema glaciale" = (IMR_catch_clean$catchcount_Myctophiformes
      + IMR_catch_clean$`catchcount_Benthosema glaciale`),
    "catchweight_g_Benthosema glaciale" = (IMR_catch_clean$catchweight_g_Myctophiformes
      + IMR_catch_clean$`catchweight_g_Benthosema glaciale`),
    "catchcount_Euphausiacea" = (IMR_catch_clean$catchcount_Euphausiacea
      + IMR_catch_clean$`catchcount_Meganyctiphanes norvegica`),
    "catchweight_g_Euphausiacea" = (IMR_catch_clean$catchweight_g_Euphausiacea
      + IMR_catch_clean$`catchweight_g_Meganyctiphanes norvegica`),
    "catchweight_g_Other Crustaceans" = (IMR_catch_clean$catchweight_g_Crustacea
      + IMR_catch_clean$catchweight_g_Decapoda)
  )


# rename Cyanea capillata to Jellyfish sp.
colnames(IMR_catch_clean)[colnames(IMR_catch_clean) == "catchcount_Cyanea capillata"] <- "catchcount_Jellyfish spp"
colnames(IMR_catch_clean)[colnames(IMR_catch_clean) == "catchweight_g_Cyanea capillata"] <- "catchweight_g_Jellyfish spp"
# rename catchcount_Decapoda to catchcount_Other Crustaceans
colnames(IMR_catch_clean)[colnames(IMR_catch_clean) == "catchcount_Decapoda"] <- "catchcount_Other Crustaceans"

IMR_catch_clean[, grepl("Jellyfish", names(IMR_catch_clean))]


# remove columns that are merged:
IMR_catch_clean <- IMR_catch_clean %>%
  subset(select = -c(
    `catchcount_Todaropsis eblanae`, `catchweight_g_Todaropsis eblanae`,
    `catchcount_Rossia macrosoma`, `catchweight_g_Rossia macrosoma`,
    `catchcount_Atlantopandalus propinqvus`,
    `catchweight_g_Atlantopandalus propinqvus`,
    `catchcount_Dichelopandalus bonnieri`, `catchweight_g_Dichelopandalus bonnieri`, catchcount_Myctophiformes, catchweight_g_Myctophiformes,
    `catchcount_Pandalus borealis`, `catchweight_g_Pandalus borealis`,
    catchcount_Pontophilus, catchweight_g_Pontophilus,
    catchweight_g_Pandalina,
    `catchcount_Sergestes arcticus`, `catchweight_g_Sergestes arcticus`,
    catchweight_g_Sergestes,
    catchcount_Pasiphaea, `catchcount_Pasiphaea multidentata`,
    `catchcount_Pasiphaea sivado`, `catchcount_Pasiphaea tarda`,
    catchcount_Pasiphaeidae, catchweight_g_Pasiphaea,
    `catchweight_g_Pasiphaea multidentata`, `catchweight_g_Pasiphaea sivado`,
    `catchweight_g_Pasiphaea tarda`, catchweight_g_Pasiphaeidae,
    catchcount_Munida, `catchcount_Munida rugosa`,
    catchweight_g_Munida, `catchweight_g_Munida rugosa`,
    `catchcount_Meganyctiphanes norvegica`,
    `catchweight_g_Meganyctiphanes norvegica`,
    catchcount_Aegidae, catchweight_g_Aegidae,
    catchweight_g_Crustacea, catchweight_g_Decapoda
  ))


###
#### 2: UiB catch and CTD ####
###

###### 2.1: load catch data 2011-2018 + 2022 ####

## 2011
uib_2011 <- read_excel("_data/biotic/2011_BIO325_Benthic_data_NDG.xlsx") %>%
  subset(select = -c(starttime.local, endtime.local)) # remove to be able to merge later
uib_2011 <- uib_2011[-7, ] # remove NA row

# add lat,lon and sp. to df. only get it to work with as.numeric
extracolumns_2011 <- tibble(
  Ship.station =
    as.numeric(c("127", "117", "115", "111", "109", "107")),
  latitude = as.numeric(c(
    60.87268, 60.87194, 60.87246,
    60.87272, 60.87267, 60.87255
  )),
  longitude = as.numeric(c(
    5.39864, 5.40193, 5.40270,
    5.39716, 5.39774, 5.39354
  )),
  M.norvegica.g = as.numeric(c(16, 3, 30, 2.5, 26, 0)),
  M.norvegica.no = as.numeric(c(41, 10, 79, 7, 76, 0)),
  S.arcticus.g = as.numeric(c(168, 1150, 218, 184, 314, 0)),
  S.arcticus.no = as.numeric(c(154, 817, 206, 179, 266, 0)),
  Shrimp.sp.g = as.numeric(c(21, 7, 2.5, 0, 0, 0)),
  Shrimp.sp.no = as.numeric(c(6, 2, 2, 0, 0, 0))
)

uib_2011 <- inner_join(uib_2011, extracolumns_2011, by = "Ship.station")

# insert values in already existing column
uib_2011$P.tarda.g[uib_2011$Ship.station == "127"] <- 795
uib_2011$P.tarda.no[uib_2011$Ship.station == "127"] <- 120

uib_2011$P.tarda.g[uib_2011$Ship.station == "117"] <- 1340
uib_2011$P.tarda.no[uib_2011$Ship.station == "117"] <- 197

uib_2011$P.tarda.g[uib_2011$Ship.station == "115"] <- 170
uib_2011$P.tarda.no[uib_2011$Ship.station == "115"] <- 29

uib_2011$P.tarda.g[uib_2011$Ship.station == "109"] <- 271
uib_2011$P.tarda.no[uib_2011$Ship.station == "109"] <- 55


## 2012
uib_2012 <- read_excel("_data/biotic/2012_BIO325_Benthic_data_NDG.xlsx") %>%
  subset(select = -c(starttime.local, endtime.local)) # remove to be able to merge later - wrong format
uib_2012 <- uib_2012[-4, ]
colnames(uib_2012)[colnames(uib_2012) == "month"] <- "mnd"

# add lat lon to df.
latlon_2012 <- tibble(
  Ship.station =
    as.numeric(c("366", "363", "372")),
  latitude = as.numeric(c(60.87183, 60.87317, 60.87233)),
  longitude = as.numeric(c(5.39017, 5.39683, 5.38533))
)

uib_2012 <- inner_join(uib_2012, latlon_2012)

## 2013
uib_2013 <- read_excel("_data/biotic/2013_BIO325_Benthic_data_NDG.xlsx") %>%
  subset(select = -c(starttime.local, endtime.local)) # remove to be able to merge later
uib_2013 <- uib_2013[-3, ]

# add lat lon to df.
latlon_2013 <- tibble(
  Ship.station =
    as.numeric(c("188", "198")),
  latitude = as.numeric(c(60.87314, 60.87356)),
  longitude = as.numeric(c(5.40993, 5.41459))
)

uib_2013 <- inner_join(uib_2013, latlon_2013)


## 2014
uib_2014 <- read_excel("_data/biotic/2014_BIO325_Benthic_data_NDG.xlsx") %>%
  subset(select = -c(starttime.local, endtime.local)) # remove to be able to merge later
uib_2014 <- uib_2014[-c(3, 4, 5), ]

# add lat lon to df.
latlon_2014 <- tibble(
  Ship.station =
    as.numeric(c("259", "258")),
  latitude = as.numeric(c(60.872, 60.87233)),
  longitude = as.numeric(c(5.410, 5.41733))
)

uib_2014 <- inner_join(uib_2014, latlon_2014)

## 2015 - NEED TO ADD LAT LON
uib_2015 <- read_excel("_data/biotic/2015_BIO325_Benthic_data_NDG.xlsx") %>%
  subset(select = -c(starttime.local, endtime.local)) # remove to be able to merge later

# add lat lon to df.
latlon_2015 <- tibble(
  Ship.station =
    as.numeric(c("111", "115", "118", "120", "128")),
  latitude = as.numeric(c(
    60.8771759,
    60.8773495,
    60.8738653,
    60.8749217,
    60.8740067
  )),
  longitude = as.numeric(c(
    5.4466318,
    5.4506083,
    5.4220284,
    5.4359760,
    5.4263521
  ))
)

uib_2015 <- inner_join(uib_2015, latlon_2015)

## 2016 (previously col named Ship.station now only station)
uib_2016 <- read_excel("_data/biotic/2016_BIO325_Benthic_data_NDG.xlsx") %>%
  subset(select = -fishing.time.min) %>%
  subset(select = -c(starttime.local, endtime.local)) # remove to be able to merge later
colnames(uib_2016)[colnames(uib_2016) == "station"] <- "Ship.station"
colnames(uib_2016)[colnames(uib_2016) == "Comments"] <- "comments.to.others"

# new column for fishing time
newcolumns_2016 <- tibble(
  Ship.station = as.numeric(c(
    "155", "156", "157",
    "160", "161", "163",
    "164", "173", "180"
  )),
  fishing.time.min = as.numeric(c(
    30, 28, 40, 26, 11,
    28, 41, 31, 20
  )),
  bottomdepth = as.numeric(c(
    635, 373, 374,
    373, 58, 374,
    651, 646, 372
  ))
)
uib_2016 <- inner_join(uib_2016, newcolumns_2016, by = "Ship.station")

## 2017
uib_2017 <- read_excel("_data/biotic/2017_BIO325_Benthic_data_NDG.xlsx") %>%
  subset(select = -c(starttime.local, endtime.local)) # remove to be able to merge later
colnames(uib_2017)[colnames(uib_2017) == "station"] <- "Ship.station"
colnames(uib_2017)[colnames(uib_2017) == "Comments"] <- "comments.to.others"

bottomdepth_2017 <- tibble(
  Ship.station = as.numeric(c("426", "430", "436")),
  bottomdepth = as.numeric(c(373, 373, 373))
)
uib_2017 <- inner_join(uib_2017, bottomdepth_2017, by = "Ship.station")

## 2018
uib_2018 <- read_excel("_data/biotic/2018_BIO325_Benthic_data_NDG.xlsx") %>%
  subset(select = -c(starttime.local, endtime.local)) # remove to be able to merge later
colnames(uib_2018)[colnames(uib_2018) == "station"] <- "Ship.station"
colnames(uib_2018)[colnames(uib_2018) == "Comments"] <- "comments.to.others"

bottomdepth_2018 <- tibble(
  Ship.station = as.numeric(c("564", "567", "575")),
  bottomdepth = as.numeric(c(366, 373, 372))
)
uib_2018 <- inner_join(uib_2018, bottomdepth_2018, by = "Ship.station")

## 2022
uib_2022 <- read_excel("_data/biotic/2022_CLIFORD_Benthic_data_NDG.xlsx") %>%
  subset(select = -c(starttime.local, endtime.local)) # remove to be able to merge later
colnames(uib_2022)[colnames(uib_2022) == "station"] <- "Ship.station"
# changing Mean.depth to bottomdepth to compare to other dfs
colnames(uib_2022)[colnames(uib_2022) == "Mean.depth"] <- "bottomdepth"
colnames(uib_2022)[colnames(uib_2022) == "month"] <- "mnd"
colnames(uib_2022)[colnames(uib_2022) == "data.type"] <- "trawl"
colnames(uib_2022)[colnames(uib_2022) == "totalcatch.g"] <- "tot.catch.g"
colnames(uib_2022)[colnames(uib_2022) == "fjord"] <- "area"
# changing CTD varaible names to match IMR data
colnames(uib_2022)[colnames(uib_2022) == "temperature.C"] <- "Temperature"
colnames(uib_2022)[colnames(uib_2022) == "salinity"] <- "Salinity"
colnames(uib_2022)[colnames(uib_2022) == "oxygen.ml_l"] <- "Oxygen"

###### 2.2: UiB CTD data 2011-2018 ####

## 2011
files_2011 <- list.files(path = "./_data/ctd/2011", pattern = ".cnv")
ctd_df_2011 <- data.table()
for (i_2011 in seq_along(files_2011)) {
  df_2011 <- read.ctd(paste0("./_data/ctd/2011/", files_2011[i_2011]))
  df_2011 <- ctdTrim(df_2011, method = "downward")
  df_2011 <- df_2011@data %>%
    as.data.table() %>%
    mutate(
      lon = df_2011@metadata$longitude,
      lat = df_2011@metadata$latitude,
      date = df_2011@metadata$startTime,
      ctdstation = df_2011@metadata$station,
      ship = df_2011@metadata$ship,
      filename = files_2011[i_2011]
    )
  ctd_df_2011 <- rbind(df_2011, ctd_df_2011)
}

ctd_df_2011$Depth <- ctd_df_2011$pressure

ctd_bottom_2011 <- ctd_df_2011 %>%
  group_by(ctdstation) %>%
  dplyr::summarise(
    Bottom_depth_CTD = max(Depth),
    Temperature = mean(temperature[Depth > 0.99 * max(Depth)]),
    Salinity = mean(salinity[Depth > 0.99 * max(Depth)]),
    Oxygen = mean(oxygen[Depth > 0.99 * max(Depth)])
  )

# closest CTD stations 2011
zone_2011 <- floor((mean(subset(uib_2011)$longitude, na.rm = T) + 180) / 6) + 1
utmCRS_2011 <- CRS(paste0("+proj=utm +zone=", zone_2011, " +datum=WGS84 +units=km +no_defs"))

stationSP_2011 <- SpatialPoints(
  uib_2011 %>% select(longitude, latitude),
  CRS("+proj=longlat +datum=WGS84")
)
stationSP_2011 <- spTransform(stationSP_2011, utmCRS_2011)

ctdSP_2011 <- SpatialPoints(
  ctd_df_2011 %>% distinct(lon, lat),
  CRS("+proj=longlat +datum=WGS84")
)
ctdSP_2011 <- spTransform(ctdSP_2011, utmCRS_2011)


dist_2011 <- gDistance(ctdSP_2011, stationSP_2011, byid = T) # distances between stations
minDist_2011 <- apply(dist_2011, 1, function(x) order(x, decreasing = F)[1]) # find closest
ctd_match_2011 <- data.frame(
  Ship.station = unique(uib_2011$Ship.station),
  ctdstation = unique(ctd_df_2011$ctdstation)[minDist_2011],
  distance = dist_2011[minDist_2011]
)

ctd_bottom_2011 <- inner_join(ctd_bottom_2011, ctd_match_2011, by = "ctdstation")

uib_2011 <- left_join(uib_2011, ctd_bottom_2011, by = "Ship.station")


## 2012
files_2012 <- list.files(path = "./_data/ctd/2012", pattern = ".cnv")
ctd_df_2012 <- data.table()
for (i_2012 in seq_along(files_2012)) {
  df_2012 <- read.ctd(paste0("./_data/ctd/2012/", files_2012[i_2012]))
  df_2012 <- ctdTrim(df_2012, method = "downward")
  df_2012 <- df_2012@data %>%
    as.data.table() %>%
    mutate(
      lon = df_2012@metadata$longitude,
      lat = df_2012@metadata$latitude,
      date = df_2012@metadata$startTime,
      ctdstation = df_2012@metadata$station,
      ship = df_2012@metadata$ship,
      filename = files_2012[i_2012]
    )
  ctd_df_2012 <- rbind(df_2012, ctd_df_2012)
}

ctd_df_2012$Depth <- ctd_df_2012$pressure

ctd_bottom_2012 <- ctd_df_2012 %>%
  group_by(ctdstation) %>%
  dplyr::summarise(
    Bottom_depth_CTD = max(Depth),
    Temperature = mean(temperature[Depth > 0.99 * max(Depth)]),
    Salinity = mean(salinity[Depth > 0.99 * max(Depth)]),
    Oxygen = mean(oxygen[Depth > 0.99 * max(Depth)])
  )

# closest CTD stations 2012
zone_2012 <- floor((mean(subset(uib_2012)$longitude, na.rm = T) + 180) / 6) + 1
utmCRS_2012 <- CRS(paste0("+proj=utm +zone=", zone_2012, " +datum=WGS84 +units=km +no_defs"))

stationSP_2012 <- SpatialPoints(
  uib_2012 %>% select(longitude, latitude),
  CRS("+proj=longlat +datum=WGS84")
)
stationSP_2012 <- spTransform(stationSP_2012, utmCRS_2012)

ctdSP_2012 <- SpatialPoints(
  ctd_df_2012 %>% distinct(lon, lat),
  CRS("+proj=longlat +datum=WGS84")
)
ctdSP_2012 <- spTransform(ctdSP_2012, utmCRS_2012)


dist_2012 <- gDistance(ctdSP_2012, stationSP_2012, byid = T) # distances between stations
minDist_2012 <- apply(dist_2012, 1, function(x) order(x, decreasing = F)[1]) # find closest
ctd_match_2012 <- data.frame(
  Ship.station = unique(uib_2012$Ship.station),
  ctdstation = unique(ctd_df_2012$ctdstation)[minDist_2012],
  distance = dist_2012[minDist_2012]
)


ctd_bottom_2012 <- inner_join(ctd_bottom_2012, ctd_match_2012, by = "ctdstation")

uib_2012 <- left_join(uib_2012, ctd_bottom_2012, by = "Ship.station")


## 2013
files_2013 <- list.files(path = "./_data/ctd/2013", pattern = ".cnv")
ctd_df_2013 <- data.table()
for (i_2013 in seq_along(files_2013)) {
  df_2013 <- read.ctd(paste0("./_data/ctd/2013/", files_2013[i_2013]))
  df_2013 <- ctdTrim(df_2013, method = "downward")
  df_2013 <- df_2013@data %>%
    as.data.table() %>%
    mutate(
      lon = df_2013@metadata$longitude,
      lat = df_2013@metadata$latitude,
      date = df_2013@metadata$startTime,
      ctdstation = df_2013@metadata$station,
      ship = df_2013@metadata$ship,
      filename = files_2013[i_2013]
    )
  ctd_df_2013 <- rbind(df_2013, ctd_df_2013)
}

ctd_df_2013$Depth <- ctd_df_2013$pressure

ctd_bottom_2013 <- ctd_df_2013 %>%
  group_by(ctdstation) %>%
  dplyr::summarise(
    Bottom_depth_CTD = max(Depth),
    Temperature = mean(temperature[Depth > 0.99 * max(Depth)]),
    Salinity = mean(salinity[Depth > 0.99 * max(Depth)]),
    Oxygen = mean(oxygen[Depth > 0.99 * max(Depth)])
  )

# closest CTD stations 2013
zone_2013 <- floor((mean(subset(uib_2013)$longitude, na.rm = T) + 180) / 6) + 1
utmCRS_2013 <- CRS(paste0("+proj=utm +zone=", zone_2013, " +datum=WGS84 +units=km +no_defs"))

stationSP_2013 <- SpatialPoints(
  uib_2013 %>% select(longitude, latitude),
  CRS("+proj=longlat +datum=WGS84")
)
stationSP_2013 <- spTransform(stationSP_2013, utmCRS_2013)

ctdSP_2013 <- SpatialPoints(
  ctd_df_2013 %>% distinct(lon, lat),
  CRS("+proj=longlat +datum=WGS84")
)
ctdSP_2013 <- spTransform(ctdSP_2013, utmCRS_2013)


dist_2013 <- gDistance(ctdSP_2013, stationSP_2013, byid = T) # distances between stations
minDist_2013 <- apply(dist_2013, 1, function(x) order(x, decreasing = F)[1]) # find closest
ctd_match_2013 <- data.frame(
  Ship.station = unique(subset(uib_2013)$Ship.station),
  ctdstation = unique(ctd_df_2013$ctdstation)[minDist_2013],
  distance = dist_2013[minDist_2013]
)

ctd_bottom_2013 <- inner_join(ctd_bottom_2013, ctd_match_2013, by = "ctdstation")

uib_2013 <- left_join(uib_2013, ctd_bottom_2013, by = "Ship.station")



## 2014
files_2014 <- list.files(path = "./_data/ctd/2014", pattern = ".cnv")
ctd_df_2014 <- data.table()
for (i_2014 in seq_along(files_2014)) {
  df <- read.ctd(paste0("./_data/ctd/2014/", files_2014[i_2014]))
  df <- ctdTrim(df, method = "downward")
  df <- df@data %>%
    as.data.table() %>%
    mutate(
      lon = df@metadata$longitude, lat = df@metadata$latitude,
      date = df@metadata$startTime, ctdstation = df@metadata$station,
      ship = df@metadata$ship, filename = files_2014[i_2014]
    )
  ctd_df_2014 <- rbind(df, ctd_df_2014)
}

ctd_df_2014$Depth <- ctd_df_2014$pressure

ctd_bottom_2014 <- ctd_df_2014 %>%
  group_by(ctdstation) %>%
  dplyr::summarise(
    Bottom_depth_CTD = max(Depth),
    Temperature = mean(temperature[Depth > 0.99 * max(Depth)]),
    Salinity = mean(salinity[Depth > 0.99 * max(Depth)]),
    Oxygen = mean(oxygen[Depth > 0.99 * max(Depth)])
  )

# closest CTD stations 2014
zone_2014 <- floor((mean(subset(uib_2014)$longitude, na.rm = T) + 180) / 6) + 1
utmCRS_2014 <- CRS(paste0("+proj=utm +zone=", zone_2014, " +datum=WGS84 +units=km +no_defs"))

stationSP_2014 <- SpatialPoints(
  uib_2014 %>% select(longitude, latitude),
  CRS("+proj=longlat +datum=WGS84")
)
stationSP_2014 <- spTransform(stationSP_2014, utmCRS_2014)

ctdSP_2014 <- SpatialPoints(
  ctd_df_2014 %>% distinct(lon, lat),
  CRS("+proj=longlat +datum=WGS84")
)
ctdSP_2014 <- spTransform(ctdSP_2014, utmCRS_2014)


dist_2014 <- gDistance(ctdSP_2014, stationSP_2014, byid = T) # distances between stations
minDist_2014 <- apply(dist_2014, 1, function(x) order(x, decreasing = F)[1]) # find closest
ctd_match_2014 <- data.frame(
  Ship.station = unique(uib_2014$Ship.station),
  ctdstation = unique(ctd_df_2014$ctdstation)[minDist_2014],
  distance = dist_2014[minDist_2014]
)

ctd_bottom_2014 <- inner_join(ctd_bottom_2014, ctd_match_2014, by = "ctdstation")

uib_2014 <- left_join(uib_2014, ctd_bottom_2014, by = "Ship.station")



## 2015
files <- list.files(path = "./_data/ctd/2015", pattern = ".cnv")
ctd_df_2015 <- data.table()
for (i in seq_along(files)) {
  df <- read.ctd(paste0("./_data/ctd/2015/", files[i]))
  df <- ctdTrim(df, method = "downward")
  df <- df@data %>%
    as.data.table() %>%
    mutate(
      lon = df@metadata$longitude, lat = df@metadata$latitude,
      date = df@metadata$startTime, ctdstation = df@metadata$station,
      ship = df@metadata$ship, filename = files[i]
    )
  ctd_df_2015 <- rbind(df, ctd_df_2015)
}

ctd_df_2015$Depth <- ctd_df_2015$pressure

ctd_bottom_2015 <- ctd_df_2015 %>%
  group_by(ctdstation) %>%
  dplyr::summarise(
    Bottom_depth_CTD = max(Depth),
    Temperature = mean(temperature[Depth > 0.99 * max(Depth)]),
    Salinity = mean(salinity[Depth > 0.99 * max(Depth)]),
    Oxygen = mean(oxygen[Depth > 0.99 * max(Depth)])
  )

# closest CTD stations 2015 - NEED TO FIND LON + LAT FOR STATION DF
zone_2015 <- floor((mean(subset(uib_2015)$longitude, na.rm = T) + 180) / 6) + 1
utmCRS_2015 <- CRS(paste0("+proj=utm +zone=", zone_2015, " +datum=WGS84 +units=km +no_defs"))

stationSP_2015 <- SpatialPoints(
  uib_2015 %>% select(longitude, latitude),
  CRS("+proj=longlat +datum=WGS84")
)
stationSP_2015 <- spTransform(stationSP_2015, utmCRS_2015)

ctdSP_2015 <- SpatialPoints(
  ctd_df_2015 %>% distinct(lon, lat),
  CRS("+proj=longlat +datum=WGS84")
)
ctdSP_2015 <- spTransform(ctdSP_2015, utmCRS_2015)


dist_2015 <- gDistance(ctdSP_2015, stationSP_2015, byid = T) # distances between stations
minDist_2015 <- apply(dist_2015, 1, function(x) order(x, decreasing = F)[1]) # find closest
ctd_match_2015 <- data.frame(
  Ship.station = unique(uib_2015$Ship.station),
  ctdstation = unique(ctd_df_2015$ctdstation)[minDist_2015],
  distance = dist_2015[minDist_2015]
)

ctd_bottom_2015 <- inner_join(ctd_bottom_2015, ctd_match_2015, by = "ctdstation")

uib_2015 <- left_join(uib_2015, ctd_bottom_2015, by = "Ship.station")



## 2016
files <- list.files(path = "./_data/ctd/2016", pattern = ".cnv")
ctd_df_2016 <- data.table()
for (i in seq_along(files)) {
  df <- read.ctd(paste0("./_data/ctd/2016/", files[i]))
  df <- ctdTrim(df, method = "downward")
  df <- df@data %>%
    as.data.table() %>%
    mutate(
      lon = df@metadata$longitude, lat = df@metadata$latitude,
      date = df@metadata$startTime, ctdstation = df@metadata$station,
      ship = df@metadata$ship, filename = files[i]
    )
  ctd_df_2016 <- rbind(df, ctd_df_2016)
}

ctd_df_2016$Depth <- ctd_df_2016$pressure

ctd_bottom_2016 <- ctd_df_2016 %>%
  group_by(ctdstation) %>%
  dplyr::summarise(
    Bottom_depth_CTD = max(Depth),
    Temperature = mean(temperature[Depth > 0.99 * max(Depth)]),
    Salinity = mean(salinity[Depth > 0.99 * max(Depth)]),
    Oxygen = mean(oxygen[Depth > 0.99 * max(Depth)])
  )

# closest CTD stations 2016

uib_2016$longitude <- as.numeric(as.character(uib_2016$longitude))
uib_2016$latitude <- as.numeric(as.character(uib_2016$latitude))
# sapply(uib_2016, class)

zone_2016 <- floor((mean(subset(uib_2016)$longitude, na.rm = T) + 180) / 6) + 1
utmCRS_2016 <- CRS(paste0("+proj=utm +zone=", zone_2016, " +datum=WGS84 +units=km +no_defs"))

stationSP_2016 <- SpatialPoints(
  uib_2016 %>% select(longitude, latitude),
  CRS("+proj=longlat +datum=WGS84")
)
stationSP_2016 <- spTransform(stationSP_2016, utmCRS_2016)

ctdSP_2016 <- SpatialPoints(
  ctd_df_2016 %>% distinct(lon, lat),
  CRS("+proj=longlat +datum=WGS84")
)
ctdSP_2016 <- spTransform(ctdSP_2016, utmCRS_2016)


dist_2016 <- gDistance(ctdSP_2016, stationSP_2016, byid = T) # distances between stations
minDist_2016 <- apply(dist_2016, 1, function(x) order(x, decreasing = F)[1]) # find closest
ctd_match_2016 <- data.frame(
  Ship.station = unique(uib_2016$Ship.station),
  ctdstation = unique(ctd_df_2016$ctdstation)[minDist_2016],
  distance = dist_2016[minDist_2016]
)

ctd_bottom_2016 <- inner_join(ctd_bottom_2016, ctd_match_2016, by = "ctdstation")

uib_2016 <- left_join(uib_2016, ctd_bottom_2016, by = "Ship.station")


## 2017
files <- list.files(path = "./_data/ctd/2017", pattern = ".cnv")
ctd_df_2017 <- data.table()
for (i in seq_along(files)) {
  df <- read.ctd(paste0("./_data/ctd/2017/", files[i]))
  df <- ctdTrim(df, method = "downward")
  df <- df@data %>%
    as.data.table() %>%
    mutate(
      lon = df@metadata$longitude, lat = df@metadata$latitude,
      date = df@metadata$startTime, ctdstation = df@metadata$station,
      ship = df@metadata$ship, filename = files[i]
    )
  ctd_df_2017 <- rbind(df, ctd_df_2017)
}

ctd_df_2017$Depth <- ctd_df_2017$pressure

ctd_bottom_2017 <- ctd_df_2017 %>%
  group_by(ctdstation) %>%
  dplyr::summarise(
    Bottom_depth_CTD = max(Depth),
    Temperature = mean(temperature[Depth > 0.99 * max(Depth)]),
    Salinity = mean(salinity[Depth > 0.99 * max(Depth)]),
    Oxygen = mean(oxygen[Depth > 0.99 * max(Depth)])
  )

# closest CTD stations 2017
zone_2017 <- floor((mean(subset(uib_2017)$longitude, na.rm = T) + 180) / 6) + 1
utmCRS_2017 <- CRS(paste0("+proj=utm +zone=", zone_2017, " +datum=WGS84 +units=km +no_defs"))

stationSP_2017 <- SpatialPoints(
  uib_2017 %>% select(longitude, latitude),
  CRS("+proj=longlat +datum=WGS84")
)
stationSP_2017 <- spTransform(stationSP_2017, utmCRS_2017)

ctdSP_2017 <- SpatialPoints(
  ctd_df_2017 %>% distinct(lon, lat),
  CRS("+proj=longlat +datum=WGS84")
)
ctdSP_2017 <- spTransform(ctdSP_2017, utmCRS_2017)


dist_2017 <- gDistance(ctdSP_2017, stationSP_2017, byid = T) # distances between stations
minDist_2017 <- apply(dist_2017, 1, function(x) order(x, decreasing = F)[1]) # find closest
ctd_match_2017 <- data.frame(
  Ship.station = unique(uib_2017$Ship.station),
  ctdstation = unique(ctd_df_2017$ctdstation)[minDist_2017],
  distance = dist_2017[minDist_2017]
)

ctd_bottom_2017 <- inner_join(ctd_bottom_2017, ctd_match_2017, by = "ctdstation")

uib_2017 <- left_join(uib_2017, ctd_bottom_2017, by = "Ship.station")


## 2018
files <- list.files(path = "./_data/ctd/2018", pattern = ".cnv")
ctd_df_2018 <- data.table()
for (i in seq_along(files)) {
  df <- read.ctd(paste0("./_data/ctd/2018/", files[i]))
  df <- ctdTrim(df, method = "downward")
  df <- df@data %>%
    as.data.table() %>%
    mutate(
      lon = df@metadata$longitude, lat = df@metadata$latitude,
      date = df@metadata$startTime, ctdstation = df@metadata$station,
      ship = df@metadata$ship, filename = files[i]
    )
  ctd_df_2018 <- rbind(df, ctd_df_2018)
}

ctd_df_2018$Depth <- ctd_df_2018$pressure

ctd_bottom_2018 <- ctd_df_2018 %>%
  group_by(ctdstation) %>%
  dplyr::summarise(
    Bottom_depth_CTD = max(Depth),
    Temperature = mean(temperature[Depth > 0.99 * max(Depth)]),
    Salinity = mean(salinity[Depth > 0.99 * max(Depth)]),
    Oxygen = mean(oxygen[Depth > 0.99 * max(Depth)])
  )

# closest CTD stations 2018
zone_2018 <- floor((mean(subset(uib_2018)$longitude, na.rm = T) + 180) / 6) + 1
utmCRS_2018 <- CRS(paste0("+proj=utm +zone=", zone_2018, " +datum=WGS84 +units=km +no_defs"))

stationSP_2018 <- SpatialPoints(
  uib_2018 %>% select(longitude, latitude),
  CRS("+proj=longlat +datum=WGS84")
)
stationSP_2018 <- spTransform(stationSP_2018, utmCRS_2018)

ctdSP_2018 <- SpatialPoints(
  ctd_df_2018 %>% distinct(lon, lat),
  CRS("+proj=longlat +datum=WGS84")
)
ctdSP_2018 <- spTransform(ctdSP_2018, utmCRS_2018)


dist_2018 <- gDistance(ctdSP_2018, stationSP_2018, byid = T) # distances between stations
minDist_2018 <- apply(dist_2018, 1, function(x) order(x, decreasing = F)[1]) # find closest
ctd_match_2018 <- data.frame(
  Ship.station = unique(uib_2018$Ship.station),
  ctdstation = unique(ctd_df_2018$ctdstation)[minDist_2018],
  distance = dist_2018[minDist_2018]
)

ctd_bottom_2018 <- inner_join(ctd_bottom_2018, ctd_match_2018, by = "ctdstation")

uib_2018 <- left_join(uib_2018, ctd_bottom_2018, by = "Ship.station")


## 2022 - already have paired ctd stations to the trawl df
##
files <- list.files(path = "./_data/ctd/2022_CLIFORD_Benthic_data_CTDs", pattern = ".cnv")
ctd_df_2022 <- data.table()
for (i in seq_along(files)) {
  df <- read.ctd(paste0("./_data/ctd/2022_CLIFORD_Benthic_data_CTDs/", files[i]))
  df <- ctdTrim(df, method = "downward")
  df <- df@data %>%
    as.data.table() %>%
    mutate(
      lon = df@metadata$longitude, lat = df@metadata$latitude,
      date = df@metadata$startTime, ctdstation = df@metadata$station,
      ship = df@metadata$ship, filename = files[i]
    )
  ctd_df_2022 <- rbind(df, ctd_df_2022)
}


###### 2.3: merging UiB catch data ####

uib_catch <- rbind.fill(
  uib_2011, uib_2012, uib_2013, uib_2014,
  uib_2015, uib_2016, uib_2017, uib_2018, uib_2022
) %>%
  subset(select = -c(station, series, startdepth, enddepth, Data.Flag)) # remove as not                                                              present in 2011-2018 data

uib_catch <- cbind(ID = 1:nrow(uib_catch), uib_catch) # add ID to row - some stations is same nr

sapply(uib_catch, class)
# change Galatheidae, Pandalidae, lat and lon to numeric
uib_catch$Galatheidae.g <- as.numeric(as.character(uib_catch$Galatheidae.g))
uib_catch$Galatheidae.no <- as.numeric(as.character(uib_catch$Galatheidae.no))
uib_catch$Pandalidae.g <- as.numeric(as.character(uib_catch$Pandalidae.g))
uib_catch$Pandalidae.no <- as.numeric(as.character(uib_catch$Pandalidae.no))
uib_catch$latitude <- as.numeric(as.character(uib_catch$latitude))
uib_catch$longitude <- as.numeric(as.character(uib_catch$longitude))

uib_catch[is.na(uib_catch)] <- 0 # make NA zero

# calculate differences in CTD and trawling depth
a <- uib_catch$bottomdepth
b <- uib_catch$Bottom_depth_CTD
print(a - b)

## bottom depths trawl and CTD do not add up !!
## ID: 14 (-60), 15 (-120), 23 (-302!!) EDIT: ID 23 have later been removed - outlier
## (34, 35, 36 has no original in column Bottom_depth_CTD - so is ok)
## will find and change values manually in ctd data for 14, 15 and 23:
uib_catch$Salinity[uib_catch$ID == "14"] <- 35.0583
uib_catch$Salinity[uib_catch$ID == "15"] <- 35.0551
uib_catch$Salinity[uib_catch$ID == "23"] <- 31.5725
uib_catch$Temperature[uib_catch$ID == "14"] <- 8.1933
uib_catch$Temperature[uib_catch$ID == "15"] <- 8.1839
uib_catch$Temperature[uib_catch$ID == "23"] <- 11.7187
uib_catch$Oxygen[uib_catch$ID == "14"] <- 2.25638
uib_catch$Oxygen[uib_catch$ID == "15"] <- 2.33401
uib_catch$Oxygen[uib_catch$ID == "23"] <- 5.84603
uib_catch$Bottom_depth_CTD[uib_catch$ID == "14"] <- 400
uib_catch$Bottom_depth_CTD[uib_catch$ID == "15"] <- 340
uib_catch$Bottom_depth_CTD[uib_catch$ID == "23"] <- 58


# fix misspelled names
colnames(uib_catch)[colnames(uib_catch) == "H.hippglossus.g)"] <- "H.hippoglossus.g"

colnames(uib_catch)[colnames(uib_catch) == "C.intestinalis.g)"] <- "C.intestinalis.g."

colnames(uib_catch)[colnames(uib_catch) == "M.dypterygia .no"] <- "M.dypterygia.no."

# remove benthos
uib_catch <- uib_catch %>%
  subset(select = -c(
    Seastar.g, Seastar.no,
    P.tremulus.g, P.tremulus.no, S.tremulus.g, S.tremulus.no,
    Bivalve.g, Bivalve.no, Gastropod.g, Gastropod.no,
    Ascidiella.sp.g, Ascidiella.sp.no, Halecium.g, Halecium.no,
    Psolus.sp.g, Psolus.sp.no, P.maximus.g, P.maximus.no,
    P.maxima.g, P.maxima.no, Holothurian.g, Holothurian.no,
    O.albida.g, O.albida.no, P.andromeda.g, P.andromeda.no,
    C.intestinalis.g, C.intestinalis.no, C.intestinalis.g.,
    Brisaster.sp.g, Brisaster.sp.no, P.septemradiatum.g, P.septemradiatum.no,
    Priapulida.sp.g, Priapulida.sp.no, Ceramaster.sp.g, Ceramaster.sp.no,
    C.granularis.g, C.granularis.no, G.acutus.g, G.acutus.no,
    Anemone.orange.g, Anemone.orange.no, Echinoidea.g, Echinoidea.no,
    Actiniaria.g, Actiniaria.no, B.natans.g, B.natans.no,
    Asteroidea.g, Asteroidea.no
  ))

# find empty columns and remove them
which(colSums(uib_catch == 0) == nrow(uib_catch))

uib_catch <- uib_catch %>% # remove empty columns
  subset(select = -c(
    Panserekke.g, Panserekke.no, ...55, ...56,
    M.intestinalis.g, M.intestinalis.no,
    M.merluccius.g, M.merluccius.no, H.hippoglossus.g.
  ))

which(colSums(uib_catch == 0) == nrow(uib_catch))

## need to merge - done
# uib_catch[ , grepl( "dypterygia" , names( uib_catch ) ) ]
# uib_catch[ , grepl( "dipterygia" , names( uib_catch ) ) ]

## need to merge - done
# uib_catch[ , grepl( "argenteus" , names( uib_catch ) ) ]
# uib_catch[ , grepl( "thori" , names( uib_catch ) ) ]

## need to merge - done
# uib_catch[ , grepl( "R.cimbrius" , names( uib_catch ) ) ]
# uib_catch[ , grepl( "E.cimbrius" , names( uib_catch ) ) ]

## merge - done
# uib_catch[ , grepl( "Jellyfish" , names( uib_catch ) ) ]
# uib_catch[ , grepl( "Red.jellies" , names( uib_catch ) ) ]
# uib_catch[ , grepl( "capillata" , names( uib_catch ) ) ]

# uib_catch[ , grepl( "E.elegans" , names( uib_catch ) ) ] #wrong catchweight - 170 g (mysid)?

# merge columns that needs to be merged
uib_catch_clean <- uib_catch %>%
  mutate(
    "E.cimbrius.g" = (uib_catch$E.cimbrius.g + uib_catch$R.cimbrius.g),
    "E.cimbrius.no" = (uib_catch$E.cimbrius.no + uib_catch$R.cimbrius.no),
    "M.dypterygia.g" = (uib_catch$M.dypterygia.g + uib_catch$M.dipterygia.g),
    "M.dypterygia.no" = (uib_catch$M.dypterygia.no +
      uib_catch$M.dipterygia.no + uib_catch$M.dypterygia.no.),
    "G.thori.g" = (uib_catch$G.thori.g + uib_catch$G.argenteus.g),
    "G.thori.no" = (uib_catch$G.thori.no + uib_catch$G.argenteus.no),
    "catchweight_g_Jellyfish spp" = (uib_catch$Jellyfish.g + uib_catch$Red.jellies.g) +
      (uib_catch$C.capillata.g),
    "catchweight_g_Galatheidae" = (uib_catch$Galatheidae.g + uib_catch$Munida.g +
      uib_catch$M.tenuimana.g + uib_catch$M.sarsi.g),
    "catchcount_Galatheidae" = (uib_catch$Galatheidae.no + uib_catch$Munida.no +
      uib_catch$M.tenuimana.no + uib_catch$M.sarsi.no),
    "catchweight_g_Dendrobranchiata" = (uib_catch$Pasiphaea.g...106 +
      uib_catch$Pasiphaea.g...107
      + uib_catch$P.tarda.g
      + uib_catch$S.arcticus.g),
    "catchcount_Dendrobranchiata" = (uib_catch$P.tarda.no
      + uib_catch$S.arcticus.no),
    "catchweight_g_Caridea" = (uib_catch$Pandalidae.g
      + uib_catch$Pandalus.sp.g
      + uib_catch$P.borealis.g
      + uib_catch$C.crangon.g
      + uib_catch$Pansereke.g
      + uib_catch$Shrimp.sp.g),
    "catchcount_Caridea" = (uib_catch$Pandalidae.no
      + uib_catch$Pandalus.sp.no
      + uib_catch$P.borealis.no
      + uib_catch$C.crangon.no
      + uib_catch$Pansereke.no
      + uib_catch$Shrimp.sp.no)
  ) # missing no for 2022 uib catch


# remove the columns that are merged into new
uib_catch_clean <- uib_catch_clean %>%
  subset(select = -c(
    Galatheidae.g, Galatheidae.no,
    Munida.g, Munida.no, M.tenuimana.g,
    M.tenuimana.no, M.sarsi.g, M.sarsi.no,
    Pasiphaea.g...106, Pasiphaea.g...107, P.tarda.g,
    S.arcticus.g, S.arcticus.no,
    P.tarda.no, Jellyfish.g, Red.jellies.g, C.capillata.g,
    Jellyfish, Jellyfish.no, Red.jellies.no, C.capillata.no,
    Pandalidae.g, Pandalidae.no, Pandalus.sp.no, Pandalus.sp.g,
    P.borealis.g, P.borealis.no, C.crangon.g, C.crangon.no,
    Pansereke.g, Pansereke.no, Shrimp.sp.no, Shrimp.sp.g,
    G.argenteus.g, G.argenteus.no,
    M.dypterygia.no., M.dipterygia.g, M.dipterygia.no,
    R.cimbrius.g, R.cimbrius.no
  ))


### comments: fix values that are wrong
## ID 8, M.dypterygia.g and M.dypterygia.no values wrong - 349 g and 1 no
uib_catch_clean$M.dypterygia.g[uib_catch_clean$ID == "8"] <- 349
uib_catch_clean$M.dypterygia.no[uib_catch_clean$ID == "8"] <- 1
## ID 8, M.glutinosa.no missing, should be 1 (21g)?
uib_catch_clean$M.glutinosa.no[uib_catch_clean$ID == "8"] <- 1
## ID 17, Shrimp-like spp. 18 no, 19,3 g
uib_catch_clean$Shrimp.sp.g[uib_catch_clean$ID == "17"] <- 19.3
uib_catch_clean$Shrimp.sp.no[uib_catch_clean$ID == "17"] <- 18
## ID 27, 25 unknown shrimp 95g
uib_catch_clean$Shrimp.sp.g[uib_catch_clean$ID == "27"] <- 95
uib_catch_clean$Shrimp.sp.no[uib_catch_clean$ID == "27"] <- 25

###### 2.4: rename columns to match IMR ####

### rename column names to match IMR data!
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "ID"] <- "serialnumber"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "year"] <- "startyear"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "mnd"] <- "month"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "latitude"] <- "latitudestart"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "longitude"] <- "longitudestart"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "fishing.time.min"] <- "fishingtime_min"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "A.silus.no"] <- "catchcount_Argentina silus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "A.silus.g"] <- "catchweight_g_Argentina silus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "C.rupestris.no"] <- "catchcount_Coryphaenoides rupestris"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "C.rupestris.g"] <- "catchweight_g_Coryphaenoides rupestris"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "A.sphyraena.no"] <- "catchcount_Argentina sphyraena"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "A.sphyraena.g"] <- "catchweight_g_Argentina sphyraena"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "G.melastomus.no"] <- "catchcount_Galeus melastomus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "G.melastomus.g"] <- "catchweight_g_Galeus melastomus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "C.monstrosa.no"] <- "catchcount_Chimaera monstrosa"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "C.monstrosa.g"] <- "catchweight_g_Chimaera monstrosa"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "E.spinax.no"] <- "catchcount_Etmopterus spinax"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "E.spinax.g"] <- "catchweight_g_Etmopterus spinax"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "B.glaciale.no"] <- "catchcount_Benthosema glaciale"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "B.glaciale.g"] <- "catchweight_g_Benthosema glaciale"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.muelleri.no"] <- "catchcount_Maurolicus muelleri"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.muelleri.g"] <- "catchweight_g_Maurolicus muelleri"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "P.periphylla.no"] <- "catchcount_Periphylla"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "P.periphylla.g"] <- "catchweight_g_Periphylla"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "G.cynoglossus.no"] <- "catchcount_Glyptocephalus cynoglossus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "G.cynoglossus.g"] <- "catchweight_g_Glyptocephalus cynoglossus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.glutinosa.no"] <- "catchcount_Myxine glutinosa"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.glutinosa.g"] <- "catchweight_g_Myxine glutinosa"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Pandalidae.no"] = "catchcount_Pandalina"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Pandalidae.g"] = "catchweight_g_Pandalina"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "D.oxyrinchus.no"] <- "catchcount_Dipturus oxyrinchus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "D.oxyrinchus.g"] <- "catchweight_g_Dipturus oxyrinchus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.poutassou.no"] <- "catchcount_Micromesistius poutassou"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.poutassou.g"] <- "catchweight_g_Micromesistius poutassou"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.norvegica.no"] <- "catchcount_Euphausiacea"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.norvegica.g"] <- "catchweight_g_Euphausiacea"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "S.arcticus.no"] = "catchcount_Sergestes arcticus"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "S.arcticus.g"] = "catchweight_g_Sergestes arcticus"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Shrimp.sp.no"] = "catchcount_Shrimp sp"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Shrimp.sp.g"] = "catchweight_g_Shrimp sp"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "C.crangon.no"] = "catchcount_Crangon crangon"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "C.crangon.g"] = "catchweight_g_Crangon crangon"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Crustacea.no"] <- "catchcount_Other Crustaceans"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Crustacea.g"] <- "catchweight_g_Other Crustaceans"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "L.piscatorius.no"] <- "catchcount_Lophius piscatorius"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "L.piscatorius.g"] <- "catchweight_g_Lophius piscatorius"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Isopod.no"] <- "catchcount_Isopoda"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Isopod.g"] <- "catchweight_g_Isopoda"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "B.brosme.no"] <- "catchcount_Brosme brosme"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "B.brosme.g"] <- "catchweight_g_Brosme brosme"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "S.scombrus.no"] <- "catchcount_Scomber scombrus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "s.scombrus.g"] <- "catchweight_g_Scomber combrus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "T.trachurus.no"] <- "catchcount_Trachurus trachurus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "T.trachurus.g"] <- "catchweight_g_Trachurus trachurus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.merlangus.no"] <- "catchcount_Merlangius merlangus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.merlangus.g"] <- "catchweight_g_Merlangius merlangus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "P.norvegicus.no"] <- "catchcount_Pomatoschistus norvegicus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "P.norvegicus.g"] <- "catchweight_g_Pomatoschistus norvegicus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Hyperiidea.no"] <- "catchcount_Amphipoda"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Hyperiidea.g"] <- "catchweight_g_Amphipoda"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "S.sprattus.no"] <- "catchcount_Sprattus sprattus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "S.sprattus.g"] <- "catchweight_g_Sprattus sprattus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.dypterygia.no"] <- "catchcount_Molva dypterygia"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.dypterygia.g"] <- "catchweight_g_Molva dypterygia"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "R.nidarosiensis.g"] <- "catchweight_g_Dipturus nidarosiensis"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "R.nidarosiensis.no"] <- "catchcount_Dipturus nidarosiensis"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "H.hippoglossus.g)"] <- "catchweight_g_Hippoglossus hippoglossus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "H.hippoglossus.no"] <- "catchcount_Hippoglossus hippoglossus"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "P.borealis.g"] = "catchweight_g_Pandalus b<orealis"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "P.borealis.no"] = "catchcount_Pandalus borealis"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.aeglefinus.no"] <- "catchcount_Melanogrammus aeglefinus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.aeglefinus.g"] <- "catchweight_g_Melanogrammus aeglefinus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "C.harengus.no"] <- "catchcount_Clupea harengus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "C.harengus.g"] <- "catchweight_g_Clupea harengus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "C.lumpus.no"] <- "catchcount_Cyclopterus lumpus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "C.lumpus.g"] <- "catchweight_g_Cyclopterus lumpus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "L.limanda.no"] <- "catchcount_Limanda limanda"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "L.limanda.g"] <- "catchweight_g_Limanda limanda"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "P.flesus.no"] <- "catchcount_Platichthys flesus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "P.flesus.g"] <- "catchweight_g_Platichthys flesus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "G.morhua.no"] <- "catchcount_Gadus morhua"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "G.morhua.g"] <- "catchweight_g_Gadus morhua"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.kitt.no"] <- "catchcount_Microstomus kitt"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.kitt.g"] <- "catchweight_g_Microstomus kitt"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "T.esmarkii.no"] <- "catchcount_Trisopterus esmarkii"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "T.esmarkii.g"] <- "catchweight_g_Trisopterus esmarkii"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.merluccius.no"] <- "catchcount_Merluccius merluccius"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "M.merluccius.g"] <- "catchweight_g_Merluccius merluccius"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "E.elegans.no"] <- "catchcount_Mysid"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "E.elegans.g"] <- "catchweight_g_Mysid"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "S.norvegicus.no"] <- "catchcount_Sebastes norvegicus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "S.norvegicus.g"] <- "catchweight_g_Sebastes norvegicus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "G.aculeatus.no"] <- "catchcount_Gasterosteus aculeatus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "G.aculeatus.g"] <- "catchweight_g_Gasterosteus aculeatus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "L.whiffiagonis.no"] <- "catchcount_Lepidorhombus whiffiagonis"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "L.whiffiagonis.g"] <- "catchweight_g_Lepidorhombus whiffiagonis"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Pandalus.sp.g"] = "catchweight_g_Pandalus sp"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Pandalus.sp.no"] = "catchcount_Pandalus sp"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Pansereke.g"] = "catchweight_g_Pontophilus"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Pansereke.no"] = "catchcount_Pontophilus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "E.cimbrius.g"] <- "catchweight_g_Enchelyopus cimbrius"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "E.cimbrius.no"] <- "catchcount_Enchelyopus cimbrius"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "G.thori.g"] <- "catchweight_g_Gadiculus argenteus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "G.thori.no"] <- "catchcount_Gadiculus argenteus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "A.radiata.g"] <- "catchweight_g_Amblyraja radiata"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "A.radiata.no"] <- "catchcount_Amblyraja radiata"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Pasiphaea.no"] = "catchcount_Pasiphaea"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "Pasiphaea.g"] = "catchweight_g_Pasiphaea"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "N.norvegicus.no"] = "catchcount_Nephrops norvegicus"
# colnames(uib_catch_clean)[colnames(uib_catch_clean) == "N.norvegicus.g"] = "catchweight_g_Nephrops norvegicus"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "L.maja.no"] <- "catchcount_Lithodes maja"
colnames(uib_catch_clean)[colnames(uib_catch_clean) == "L.maja.g"] <- "catchweight_g_Lithodes maja"

# remove columns not in IMR data
uib_catch_clean_2 <- uib_catch_clean %>%
  subset(select = -c(
    area, day, Ship.station, day.or.night, comments.to.others,
    distance, CTD.station, Data.Flag.Comments, others.g, others.no
  ))



###
#### 3: merge UiB and IMR catch data and standardize ####
###

###### 3.1: merge data ####

catch_df <- rbind.fill(IMR_catch_clean, uib_catch_clean_2)

sapply(catch_df, class)
catch_df$`catchweight_g_Benthosema glaciale` <- as.numeric(as.character(catch_df$`catchweight_g_Benthosema glaciale`))
catch_df$`catchweight_g_Maurolicus muelleri` <- as.numeric(as.character(catch_df$`catchweight_g_Maurolicus muelleri`))

catch_df[is.na(catch_df)] <- 0 # make NA zero

# remove station with depth 58 m and no close CTD
catch_df <- subset(catch_df, bottomdepth != "58")
catch_df <- cbind(ID = 1:nrow(catch_df), catch_df)

catch_df <- catch_df %>%
  mutate(catchweight_total = (rowSums(select(., starts_with("catchweight"))))) %>%
  mutate(catchweight_tot_minusperiphylla_kg = (catchweight_total - catchweight_g_Periphylla) / 1000) %>%
  mutate(catchweight_kg_periphylla = catchweight_g_Periphylla / 1000) %>%
  mutate(total_catch_CPUE = (catchweight_total / fishingtime_min) / 1000)


catch_df["Trawl"] <- "C" # new column with all C values = Campelen
catch_df[catch_df$startyear == "2021", "Trawl"] <- "F" # new values for 2021 cruise stations with all Flekkerøy trawl
# print(catch_df$Trawl)

###### 3.2 standardize catch ####

species_catchweight <- catch_df %>%
  select(
    ID, fishingtime_min, longitudestart, latitudestart, startyear, Trawl,
    starts_with("catchweight")
  )

which(colnames(species_catchweight) == "catchweight_total")
which(colnames(species_catchweight) == "catchweight_tot_minusperiphylla_kg")
which(colnames(species_catchweight) == "catchweight_kg_periphylla")

species_catchweight$empty_cols <- rowSums(species_catchweight[-c(1, 2, 3, 4, 5, 6, 74, 75, 76)] == 0) # find empty columns in a row minus columns ID, fishingtime, total catch, tot catch min peri, peri catch kg
species_catchweight$number_sp <- 67 - species_catchweight$empty_cols # substract empty columns (number of species not present in the station). 68 possible species. output = species present in station
catch_df$number_sp <- species_catchweight$number_sp


## standardize by fishingtime
catchweight <- startsWith(names(species_catchweight), "catchweight")
CPUE_catchweight <- replace(species_catchweight, catchweight, species_catchweight[catchweight] / species_catchweight$fishingtime_min)

rownames(CPUE_catchweight) <- CPUE_catchweight$ID

# remove columns from df to make clean df
CPUE_catchweight_clean <- CPUE_catchweight %>%
  subset(select = -c(
    ID, fishingtime_min, latitudestart, longitudestart,
    startyear, Trawl, empty_cols, number_sp,
    catchweight_total, catchweight_tot_minusperiphylla_kg,
    catchweight_kg_periphylla
  ))

# make relative df
catchweight_relative <- (CPUE_catchweight_clean / rowSums(CPUE_catchweight_clean) * 100)
# print(CPUE_catchweight_clean)

###### 3.3: create bottom depth env_df with all variables ####

library(here)
here()
## load static variables
static_variables <- read_excel("FjordCommunities_variables_static.xlsx") %>%
  select(
    ID, dist_shallowest_sill, dist_coast_km, shallowest_sill_depth, sill_depth_m, Trawl,
    sill_category, dist_aquaculture, aquaculture_impact
  )

# bottom env
env_bottom <- catch_df %>%
  select(ID, bottomdepth, latitudestart, longitudestart, Oxygen, Temperature, Salinity)

# make env df
env_df <- inner_join(static_variables, env_bottom, by = c("ID"))
env_df <- as.data.frame(env_df)
rownames(env_df) <- env_df$ID
env_df <- env_df %>%
  subset(select = -c(ID, sill_depth_m))

# fill NAs with bottom depth to make continous variable of the sill depths
env_df <- env_df %>%
  mutate(sill_depth_m = coalesce(shallowest_sill_depth, bottomdepth)) %>%
  subset(select = -c(shallowest_sill_depth))

sapply(env_df, class)
env_df$dist_shallowest_sill <- as.numeric(env_df$dist_shallowest_sill)
env_df$dist_coast_km <- as.numeric(env_df$dist_coast_km)
env_df[is.na(env_df)] <- 0

# library(writexl)
# write_xlsx(env_df, "C:/Users/47971/Documents/FjordCommunities/FjordCommunities\\env_df.xlsx")

# catchweight and diversity to env df to use for regression models:
env_mod <- cbind(env_df, CPUE_catchweight$catchweight_tot_minusperiphylla_kg)
colnames(env_mod)[colnames(env_mod) == "CPUE_catchweight$catchweight_tot_minusperiphylla_kg"] <- "catchweight_minus_Periphylla"
env_mod$Periphylla_kg <- CPUE_catchweight$catchweight_kg_periphylla
env_mod$tot_catch <- CPUE_catchweight$catchweight_total
env_mod$ID <- CPUE_catchweight$ID
env_mod$shannon_div <- catch_df$shannon_div
env_mod$fjord_coast <- "Fjord"
env_mod <- within(env_mod, fjord_coast[fjord_coast == "Fjord" & sill_category == "1"] <- "Coastal")

## check env for correlation and outliers
# env_check <- ggpairs(env_df, cardinality_threshold = 100)
# env_check
# -> remove distance to sill and aquaculture site
env_df <- env_df %>%
  subset(select = -c(
    dist_shallowest_sill, dist_aquaculture,
    latitudestart, longitudestart, sill_category
  ))
# env_df$sill_category <- as.factor(env_df$sill_category) # categorical variable - factor


### data preparation for analysis
species_matrix <- as.matrix(catchweight_relative) # matrix of spp. relative abundance
species_matrix_nonrel <- as.matrix(CPUE_catchweight_clean)
species_matrix_sqrt <- sqrt(species_matrix) # square-root transformed relative abundance

# diversity
shannon_div <- diversity(species_matrix_nonrel, index = "shannon")
catch_df <- dplyr::mutate(catch_df, shannon_div = shannon_div)
summary(catch_df$shannon_div) # summary stat on diversity
env_mod$shannon_div <- shannon_div

#### diversity without periphylla?
CPUE_catchweight_minusperi <- subset(CPUE_catchweight_clean, select = -catchweight_g_Periphylla)

catchweight_relative_minusperi <- (CPUE_catchweight_minusperi / rowSums(CPUE_catchweight_minusperi) * 100)

species_matrix_minusperi <- as.matrix(catchweight_relative_minusperi) # matrix of spp. relative abundance
species_matrix_minusperi_sqrt <- sqrt(species_matrix_minusperi)
species_matrix_minusperi_nonrel <- as.matrix(CPUE_catchweight_minusperi)

# diversity
shannon_div_minusperi <- diversity(species_matrix_minusperi_nonrel, index = "shannon")
catch_df <- dplyr::mutate(catch_df, shannon_div_minusperi = shannon_div_minusperi)
summary(catch_df$shannon_div_minusperi) # summary stat on diversity

# save input data for analysis

save(env_mod, catch_df, file = "_data/Statistical_analysis.rda")
