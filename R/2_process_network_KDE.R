
source("R/load_packages_functions.R")

####################  read in processed data ####################

bikemaps <- read_sf("input/processed/inc_bm_201601_201709_snp_30m.gpkg")
network <- read_sf("input/processed/edge_ec_201601_201709_total.gpkg")
list_list_10m <- load("input/processed/lixel_list_edge_ec_201601_201709_total.RData")

n_cores <- bigstatsr::nb_cores()

#KDE for BikeMaps incidents

system.time({
  bm_network_kde <- network_kde(lixel_list = lixel_list_10m,
                                 point_process=bikemaps,
                                 bandwidth = c(50,100,150,250),
                                 n_cores = n_cores,
                                 point_process_is_lixel_midpoint = FALSE)
})


#KDE for BikeMaps incidents, with Moran's I significance

# system.time({
#   bm_network_kde_lisa <- network_kde_lisa(lixel_list = lixel_list_10m,
#                                  point_process=bikemaps,
#                                  bandwidth = c(50,100,250,500),
#                                  n_cores = n_cores,
#                                  point_process_is_lixel_midpoint = FALSE,
#                                  nsim=99,
#                                  spatial_lag=10
#                                  )
# })

# plot(bm_network_kde_lisa["pvalue_kde_bw_50"
#                          ])

#to do - implement different spatial lags and distances

### KDE for Bicycling Counts ###

exposure <- st_line_midpoints(lixel_list_10m$target_lixel) %>%
  left_join(network %>% st_drop_geometry(),by="ID")

system.time({
  exp_network_kde <- network_kde(
    lixel_list = lixel_list_10m,
    point_process=exposure,
    bandwidth = c(50,100,150,250),
    n_cores = n_cores,
    attribute = exposure$CAADB,
    point_process_is_lixel_midpoint = TRUE
      )
})

#write bikemaps kernel density to geopackage
write_sf(bm_network_kde$network_kde,paste0(getwd(),"/output/","bikemaps_kde_10m_lixel_50m_500m_bw.gpkg"))
#write nearest neighbour list to .RData
saveRDS(bm_network_kde$neighbours, file=paste0(getwd(),"/output/","bikemaps_neighbour_list_10m_lixel_500m.rdata"))


#write exposure kernel density to geopackage
write_sf(exp_network_kde$network_kde,paste0(getwd(),"/output/","strava_kde_10m_lixel_50m_500m_bw.gpkg"))
#write nearest neighbour list to .RData
saveRDS(exp_network_kde$neighbours, file=paste0(getwd(),"/output/","strava_neighbour_list_10m_lixel_500m.rdata"))

