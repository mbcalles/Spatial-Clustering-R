
source("R/load_packages_functions.R")
#################### 1. Read in data ####################

### read census boundary geometry

# Returns data and geography as an sf-class data frame
#CRD Census Tract

crd_ct <- get_census(dataset='CA16', regions=list(CMA="59935"),
                     vectors = c("v_CA16_5807"),
                     level='DA', use_cache = FALSE,geo_format = "sf") %>%
  st_transform(crs = 26910)  #transform projection to NAD 83 UTM Zone10N

crd_city <- crd_ct %>%
  filter(`Region Name` == "Victoria"|
           `Region Name` == "Oak Bay"|
           `Region Name` == "Saanich"|
           `Region Name` == "Esquimalt"
           ) %>%
  group_by(`Region Name`) %>%
  summarise(`.groups`="drop") %>%
  mutate(`Region Name` = as.character(`Region Name`))


vic_city <- crd_city %>%
  filter(`Region Name` == "Victoria")

crd <- crd_city %>%
  st_union()

### read strava road shapefile as simple feature
e <- read_sf("input/exposure/STRAVA/Edges/victoria_osm_edges.shp") %>%
  st_transform(crs = 26910) %>% #transform projection to NAD 83 UTM Zone10N
  filter(CLAZZ!=1) %>% #remove ferry routes
  filter(st_intersects(x=.,y=st_as_sfc(st_bbox(crd)),sparse=FALSE)) %>% # subset road network by capital region bounding box
  distinct(.keep_all = TRUE) %>%  #remove duplicate edges
  mutate(LENGTH = st_length(.))


### read in a-spatial strava counts aggregated by road geometry and month
#extract all relevant .csv filenames
files <- list.files(path = paste0(getwd(),"/input/exposure/STRAVA/Edges"),
                    pattern = "*.csv", full.names = TRUE)
#extract filenames with monthly totals for road links
monthly_total <- files[str_detect(files,"_month_") & str_detect(files,"total")]

#read in each file representing monthly counts by road geometry id
m_total <- sapply(monthly_total, read.csv,simplify = FALSE) %>%
  bind_rows(.id = "id") %>% #bind dataframes together
  mutate(id = str_remove(id,
                         paste0(getwd(),"/input/exposure/STRAVA/Edges/victoria_201601_201709_ride_rollup_month_"))) %>%
  mutate(id = str_remove(id,
                         ".csv"))

### read bikemaps data

bm_all <- readLines('https://bikemaps.org/incidents.json') %>%  #Read in data from website
  st_read() %>% #convert to spatial object
  st_transform(crs = 26910) %>%
  filter(st_intersects(x=.,y=st_as_sfc(st_bbox(crd)),sparse=FALSE)) %>%  # subset bikemaps by captial region boundary
  filter(date >= ymd("2014-10-01"))


#################### 2. Aggregate Exposure/Crashes by Road Network Link ####################

###Aggregate monthly counts by edge id
total_by_edge_id <- m_total %>%
    group_by(edge_id) %>%
    summarise(n_months = n(),
              t_tactcnt = sum(tactcnt))

### Join total counts to geometry
e_total <- e %>%
  left_join(total_by_edge_id,by = c("ID"="edge_id")) %>%
  replace_na(list(t_tactcnt=0)) %>%   #replace na's with 0s (e.g. edge not in the strava data
  mutate(avg_monthly_tactn = t_tactcnt/21,
         avg_annual_daily_tactn = t_tactcnt/638) %>%
  mutate(CAADB = 40 + avg_annual_daily_tactn*50,
         CAADB_KM = CAADB*LENGTH) #(Roy et al 2019/Jestico et al 2016)


### Conflate bikemaps incidents with road geometry
bm_snp <- st_snap_points(bm_all,e_total,max_dist = 30)
bm_snp <- st_sf(bm_all %>%
                  st_drop_geometry() %>%
                  mutate(geom=bm_snp))
#join edge id to each point
bm_eID <- st_join(bm_snp,e_total["ID"],
                  join = st_is_within_distance, dist = 1)
#aggregate points by edge id
eID_sum_crashes <- bm_eID %>%
  filter(!is.na(ID)) %>%
  group_by(ID) %>%
  summarise(n_crashes = sum(p_type=="collision"),
            n_nearmiss = sum(p_type=="nearmiss")) %>%
  ungroup() %>%
  st_drop_geometry()

#join counts by edge id
e_total <- left_join(e_total,eID_sum_crashes ,
                       by = "ID")  %>%
  replace_na(list(n_crashes=0,n_nearmiss=0))

e_total <- e_total %>%
  mutate(
  all_incidents = n_crashes + n_nearmiss)

#################### 3. Create Lixelized Network ####################

network_sln <- SpatialLinesNetwork(e_total)
network_sln <- sln_clean_graph(network_sln)# Remove unconnected roads

### Lixelize cleaned network data

lixel_list_10m <- lixelize_network(
  sf_network = network_sln@sl,
  max_lixel_length = 10
)


#################### 4. Write Data ####################

st_write(vic_city, paste0(getwd(), "/input/processed/","city_of_victoria.gpkg"),)

st_write(crd, paste0(getwd(), "/input/processed/","study_area.gpkg"))

st_write(crd_city, paste0(getwd(), "/input/processed/","study_area_city.gpkg"))


st_write(e_total,paste0(getwd(), "/input/processed/","edge_ec_total_crash_strava.gpkg"))

st_write(bm_snp,paste0(getwd(), "/input/processed/","inc_bm_snp_30m.gpkg"))

save(lixel_list_10m, file = paste0(getwd(), "/input/processed/","lixel_list_edge_ec_total_crash_strava.RData"))

library(beepr)

beep(sound = 8)
