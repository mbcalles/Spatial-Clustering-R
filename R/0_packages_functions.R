#################### Load Packages ####################
library(doParallel)
library(foreach)
library(igraph)
library(tidygraph)
library(stplanr)
library(spdep)
library(sf)
library(dplyr)
library(stringr)
library(cancensus)
library(lubridate)
library(readr)
library(tidyr)

options(cancensus.api_key = "CensusMapper_96346368de96ca68bd3f0127e25013d0")

#################### Load Functions ####################
file.sources = list.files(path = paste0(getwd(),"/R/functions/"),pattern="*.R$",full.names = TRUE)
sapply(file.sources,source,.GlobalEnv)
