
### Load data

network_vic <- read_sf(paste0(getwd(),"/input/processed/edge_ec_201601_201709_total.gpkg")) %>%
  distinct(.keep_all = TRUE)
bikemaps <- read_sf(paste0(getwd(),"/input/processed/inc_bm_201601_201709_snp_30m.gpkg"))
study_area <- read_sf(paste0(getwd(),"/input/processed/census_ct_ec_201601_201709_total.gpkg"))

### Clean Network Data

network_sln = SpatialLinesNetwork(network_vic)
network_sln <- sln_clean_graph(network_sln)# Remove unconnected roads

### Lixelize cleaned network data

lixel_list_10m <- lixelize_network(
  sf_network = network_sln@sl,
  max_lixel_length = 10
  )


n_cores <- bigstatsr::nb_cores()

#KDE for BikeMaps incidents

system.time({
  bm_network_kde <- network_kde(lixel_list = lixel_list_10m,
                                 point_process=bikemaps,
                                 bandwidth = c(50,100,250,500),
                                 n_cores = n_cores,
                                 point_process_is_lixel_midpoint = FALSE)
})

#KDE for BikeMaps incidents, with Moran's I significance

system.time({
  bm_network_kde_lisa <- network_kde_lisa(lixel_list = lixel_list_10m,
                                 point_process=bikemaps,
                                 bandwidth = c(50,100,250,500),
                                 n_cores = n_cores,
                                 point_process_is_lixel_midpoint = FALSE,
                                 nsim=99,
                                 spatial_lag=1
                                 )
})



### KDE for Bicycling Counts ###

exposure <- st_line_midpoints(lixel_list_10m$target_lixel) %>%
  left_join(network_vic %>% st_drop_geometry(),by="ID")

system.time({
  exp_network_kde <- network_kde(
    lixel_list = lixel_list_10m,
    point_process=exposure,
    bandwidth = c(50,100,250,500),
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

