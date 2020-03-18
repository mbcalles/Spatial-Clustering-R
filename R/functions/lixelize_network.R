

lixelize_network <- function(sf_network,max_lixel_length){
  print("Splitting input spatial lines by lixel length...")
  target_lixel <- split_lines(input_lines = sf_network,max_length = max_lixel_length,id = "ID")
  print("Create corresponding shortest distance network...")
  shortest_distance_network <- split_lines(input_lines = target_lixel,max_length = max_lixel_length/2,id = "ID")
  return(list(target_lixel=target_lixel,shortest_distance_network=shortest_distance_network))
}
