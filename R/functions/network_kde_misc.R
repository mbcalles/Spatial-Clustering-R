prank <- function(x) ifelse(is.na(x), NA, rank(x,ties.method = "min") / sum(!is.na(x)))

detect_planar_hotspot <- function(kde,percentiles){


  kde@data@values <- prank(kde@data@values)#convert density values to rank
  n_intervals <- sapply(percentiles,function(x) 1/(1-x)) # calculate number of intervals needed to isolate top x%
  hotspots <-  lapply(1:length(percentiles), function(i) raster::cut(kde,breaks = n_intervals[i]+1)) #cut raster

  for(i in 1:length(hotspots)){

    hotspots[[i]]@data@values <- ifelse(hotspots[[i]]@data@values<max(hotspots[[i]]@data@values,na.rm=TRUE),
                                        NA,
                                        hotspots[[i]]@data@values) #isolate the top x% of cells

  }

  hotspots <-  lapply(1:length(hotspots), function(i) rasterToPolygons(hotspots[[i]],dissolve = TRUE) %>% st_as_sf()) #convert to polygons

  for(i in 1:length(hotspots)){

    hotspots[[i]]$layer <- as.character(percentiles[i]) #change value to original input

  }

  hotspots <- do.call(rbind,hotspots)
  hotspots$layer <- factor(hotspots$layer, levels = as.character(percentiles))
  return(hotspots) #combine polygons

}

detect_network_hotspot <- function(network_kde,density_col,percentiles){

  if(typeof(density_col)=="character"){

    enquo_density_col <- as.name(density_col)

  }else{

    enquo_density_col <- enquo(density_col)

  }

  hotspots <- lapply(1:length(percentiles),function(x) network_kde %>% filter(prank(!!enquo_density_col)>=percentiles[x]) %>% st_union %>% st_as_sf())
  names(hotspots) <- percentiles
  hotspots <- hotspots %>% purrr::map(~mutate(.,density_col = quo_name(enquo_density_col)))
  hotspots <- do.call(what = sf:::rbind.sf, args = hotspots) %>%
    mutate(thresholds = row.names(.))
  return(hotspots)

}


