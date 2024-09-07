

#Load required packages
packages <- c("tidyverse","rgeos","rgdal","reshape2","readODS","sp","lubridate",
              "mapdata","marmap","mapplots","gridExtra","ggforce","stringr",
              "bookdown","oce",
              "readxl","lme4","devtools","mgcv","glmmTMB","jtools","data.table","sjstats",
              "RstoxBase","RstoxData","readxl","ggOceanMaps","ggpubr",
              "knitr", "kableExtra","splines")

# Install missing packages
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  if("RstoxBase" %in% packages[!installed_packages])  {
    install.packages("RstoxBase", repos = c("https://stoxproject.github.io/repo", "https://cloud.r-project.org"))## Install Rstox
  } else if ("RstoxData" %in% packages[!installed_packages] ) {
    install.packages("RstoxData", repos = c("https://stoxproject.github.io/repo/", "https://cloud.r-project.org/")) ## Install RstoxData
  } else if ("ggOceanMaps" %in% packages[!installed_packages]) {
    devtools::install_github("MikkoVihtakari/ggOceanMapsData") # Install ggOceanMaps
    devtools::install_github("MikkoVihtakari/ggOceanMaps")
  } else {
    install.packages(packages[!installed_packages])
  }
}
sapply(packages, require, character.only = TRUE)

