network_kde <- function(lixel_list,point_process,bandwidth = 100,n_cores=1,attribute=1,point_process_is_lixel_midpoint=FALSE){
  #1. Create segement base linear reference system from original road network

  #Done outside function using "lixelize_network"

  #2. Divide each segment in basic linear units of defined length (lixel_length)

  #Done outside function using "lixelize_network"
  #Our approach varies from Xie and Yan (2008) in that we define a "maximum" length for lixels, and  segments that are larger than the specified length
  #are split into equal sized lixels. Instead of having x lixels with 1 extra "residual lixel" we instead have x+1 lixels of equal length, under the lixel limit
  #E.g. a 45 lixel with a maximum lixel length of 10m  would be split into 5 segments of 9m instead of 4 lixels of 10m and a residual lixel of 5m.

  #3. Create a network of lixels by establishing the network topology between lixels as well as between lixels and lxnodes.

  #Define topology for calculating shortest distances by converting lixel with half the length to calculate distances
  #The shorter the lixel length the more accurate the calculation of network distance from nearest node of source lixel to nearest node of target lixel
  print("Defining network topology...")

  # Create nodes at the start and end point of each edge

  nodes <- lixel_list$shortest_distance_network %>%
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

  lixel_list$shortest_distance_network = lixel_list$shortest_distance_network %>%
    mutate(from = start_nodes, to = end_nodes)

  # Remove duplicate nodes
  nodes <- nodes %>%
    distinct(nodeID, .keep_all = TRUE) %>%
    select(-c(LIXID, start_end)) %>%
    st_as_sf(coords = c('X', 'Y')) %>%
    st_set_crs(st_crs(lixel_list$shortest_distance_network))

  # Convert to tbl_graph
  graph <- tbl_graph(nodes = nodes, edges = as_tibble(lixel_list$shortest_distance_network), directed = FALSE)

  graph <- graph %>%
    activate(edges) %>%
    mutate(length = st_length(geometry))

  lixel_list$shortest_distance_network <- NULL

  #4. Create the center points of all the lixels for the target lixel (lxcenters)
  print("Calculating lixel midpoints (lxcenters)...")

  lxcenters <- st_line_midpoints(lixel_list$target_lixel)

  #5. Select a point process occuring within the road network

  #input as function parameter

  #6. For each point find its nearest lixel. Count the number of points nearest to a lixel and assigned as property of lixel.

  #Points should be snapped to network within some distance threshold prior to

  if(point_process_is_lixel_midpoint==FALSE){

  print("Counting number of events on each lixel...")

  point_process <- st_join(point_process,lixel_list$target_lixel["LIXID"],
                           join = st_is_within_distance, dist = 0.001) #for each point assign the LIXID that it falls on

  source_lixels <- point_process %>% #summarize the number of points by LIXID (e.g. count the points for each LIXID)
    filter(!is.na(LIXID)) %>%
    group_by(LIXID) %>%
    summarise(n_events = n()) %>%
    ungroup() %>%
    st_drop_geometry()

  source_lixels <- inner_join(lxcenters,source_lixels,by="LIXID") #define geometry for source lixels
  }

  if(point_process_is_lixel_midpoint==TRUE){

    source_lixels <- lxcenters %>% mutate(n_events = 1)
  }

  print(paste0(sum(source_lixels$n_events)," events on ",nrow(source_lixels)," source lixels"))

  #7. Define a search bandwidth, measured with the short path network distance

  #input as function parameter

  #8. Calculate the shortest-path network distance from each source lixel to lxcenter of all its neighouring lixels within the seach bandwidth


  nearest_node_to_source_lixel <- st_nearest_feature(source_lixels,nodes)#find nodes from shortest distance network associated with each source lixel
  nearest_node_to_lxcentre <- st_nearest_feature(lxcenters,nodes)#find nodes from shortest distance network associated with each lxcenters

  print("Calculating distances from each source lixel to all other lixel centers... ")


  cl <- makeCluster(n_cores)
  registerDoParallel(cores=cl) #parallel computing

  start_time <- Sys.time()
  distances <- foreach::foreach(i = 1:length(nearest_node_to_source_lixel),.packages = c("magrittr","igraph","tidygraph")) %dopar% {
    temp <- distances(
      graph = graph,
      weights = graph %>% activate(edges) %>% pull(length),
      v = nearest_node_to_source_lixel[i]
    )[,nearest_node_to_lxcentre]

    data.frame(LIXID = lxcenters[temp<=max(bandwidth),]$LIXID,
               dist = temp[temp<=max(bandwidth)])
  }
  end_time <- Sys.time()

  print(end_time - start_time)

  stopCluster(cl)

  rm("graph")

  # gauss <-function(x) {
  #
  #   t1 <- 1/(sqrt(2*pi))
  #   t2 <- exp(-((x^2)/(2*bandwidth^2)))
  #   n <- (1/bandwidth)*(t1*t2)
  #   n <- ifelse(x>bandwidth,0,n)
  #
  #
  #   return(n)
  #
  # }
  #
  # plot(gauss(0:(bandwidth*2)))

  quartic <-function(x,r) {
    K <- 3/pi
    t1 <- 1-(x^2/r^2)
    q <- ifelse(x>r,0,K*t1)
    q <- (1/r)*q

    return(q)
  }

  LIXID <- unlist(lapply(lapply(distances,`[`,"LIXID"),function(x) pull(x)))
  distances_list <- lapply(lapply(distances,`[`,"dist"),function(x) pull(x))

  d_list <- list()

  for(i in 1:length(bandwidth)){

    print(paste0("Applying kernel density estimator with bandwidth of ",bandwidth[i],"m ... "))

    d <- lapply(distances_list,function(x) quartic(x,bandwidth[i]))
    d <- mapply(`*`,d,attribute) #multiply by attribute of source lixel
    d <- mapply(`*`,d,source_lixels$n_events) #sum over number of events that occured on the source lixel

    d_list[[i]] <- data.frame(unlist(d))
  }

  d_cols <- as.data.frame(do.call(cbind,d_list))

  names(d_cols) <- paste0("kde_bw_",bandwidth)

  #sum densities over all lixels
  density <- cbind(LIXID,d_cols) %>%
    group_by(LIXID) %>%
    summarise_each(list(~sum(.)))

  network_kde <- left_join(lixel_list$target_lixel,density,by = "LIXID") %>%
    mutate(length = st_length(.)) %>%
    replace(., is.na(.), 0)

  return(list(network_kde=network_kde,neighbours = distances))

}
