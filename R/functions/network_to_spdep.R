
sgbp <- function(x, predicate, region.id, ncol) {
  structure(x,
            predicate = predicate,
            region.id = region.id,
            ncol = ncol,
            class = "sgbp")
}


sgbp_to_nb <- function(x, ...){
  attrs <- attributes(x)
  x <- lapply(x, function(i) { if(length(i) == 0L) 0L else i } )
  attributes(x) <- attrs
  class(x) <- "nb"
  x
}

# adjacency spatial relationship
sfl_to_nb_spadj <- function(sf_network,spatial_lag = 1){

  sf_sgbp <- st_touches(sf_network,sf_network) #shared vertex (y/n)
  #convert sgbp fomrat to nb class for spdep
  sf_nb <- sgbp_to_nb(sf_sgbp)

  if(spatial_lag>1){
  sf_nb <- nblag_cumul(nblag(sf_nb,maxlag = spatial_lag))
  }
  return(sf_nb)

  }

# network distance spatial relationship - midpoint to midpoint.
sfl_to_nb_ntwrkdist <- function(sf_network,network_dist = 50,precision=10){

  sf_network_chopped <- split_lines(sf_network,max_length = precision)
  mp <- st_line_midpoints(sf_network)

  nodes <- sf_network_chopped %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(LIXID = L1) %>%
    group_by(LIXID) %>%
    slice(c(1, n())) %>%
    ungroup() %>%
    mutate(start_end = rep(c('start', 'end'), times = n()/2))

  # Give each node a unique index
  nodes <- nodes %>%
    mutate(xy = paste(.$X, .$Y)) %>%
    mutate(nodeID = group_indices(., factor(xy, levels = unique(xy)))) %>%
    select(-xy)

  # Combine the node indices with the edges
  start_nodes <- nodes %>%
    filter(start_end == 'start') %>%
    pull(nodeID)
  end_nodes <- nodes %>%
    filter(start_end == 'end') %>%
    pull(nodeID)
  sf_network_chopped = sf_network_chopped %>%
    mutate(from = start_nodes, to = end_nodes)
  # Remove duplicate nodes
  nodes <- nodes %>%
    distinct(nodeID, .keep_all = TRUE) %>%
    select(-c(LIXID, start_end)) %>%
    st_as_sf(coords = c('X', 'Y')) %>%
    st_set_crs(st_crs(sf_network_chopped))

  # Convert to tbl_graph
  graph <- tbl_graph(nodes = nodes, edges = as_tibble(sf_network_chopped), directed = FALSE)

  graph <- graph %>%
    activate(edges) %>%
    mutate(length = st_length(geom))

  nearest_node_to_midpoint <- st_nearest_feature(mp,nodes)#find nodes from shortest distance network associated with each source lixel

  cl <- makeCluster(n_cores)
  registerDoParallel(cores=cl) #parallel computing

  network_dist <- foreach::foreach(i = 1:length(nearest_node_to_midpoint),.packages = c("magrittr","igraph","tidygraph")) %dopar% {
    temp <- distances(
      graph = graph,
      weights = graph %>% activate(edges) %>% pull(length),
      v = nearest_node_to_midpoint[i]
    )[,nearest_node_to_midpoint]

    which(temp<=network_dist&temp>0)

  }
  stopCluster(cl)

  sf_sgbp <- sgbp(network_dist,"network distance",1:length(network_dist),length(network_distance))

  sf_nb <- sgbp_to_nb(sf_sgbp)

  return(sf_nb)

  }




