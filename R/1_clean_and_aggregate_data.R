
#################### 1. Read in data ####################

### read census boundary geometry

#CRD Census Tract

crd_ct <- get_census(dataset='CA16', regions=list(CMA="59935"),
                     vectors = c(bicycle="v_CA16_5807"),
                     level='CT', use_cache = FALSE, geo_format = 'sf') %>%
  st_transform(crs = 26910)  #transform projection to NAD 83 UTM Zone10N

crd_ct <- crd_ct %>%
  filter(`Region Name` == "Saanich" |
         `Region Name` == "Oak Bay" |
         `Region Name` == "Victoria" |
         `Region Name` == "Esquimalt" )

crd <- crd_ct %>%
  st_union()


### read icbc data

icbc <- read_csv(paste0(getwd(),"/input/crashes/ICBC 2014_2017/CRD_Dated.csv")) %>%
  mutate(date2 = ymd(date)) %>%
  unite(date_name,Month,Year, sep = " ",remove=FALSE) %>%
  select(-X1) %>%
  mutate(cid = row_number())

#243 crashes at 197 days and locations don't have locational information and are removed
icbc_sf <- icbc %>%
  filter(!(Longitude=="UNKNOWN" | Latitude =="UNKNOWN")) %>%
  mutate_at(vars(Latitude,Longitude),as.numeric) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326) %>%
  st_transform(crs = 26910)

icbc <- NULL

### read bikemaps data

bm <- readLines('https://bikemaps.org/incidents.json') %>%  #Read in data from website
  st_read() %>% #convert to spatial object
  st_transform(crs = 26910) %>%
  filter(st_intersects(x=.,y=st_as_sfc(st_bbox(crd)),sparse=FALSE)) %>% # subset bikemaps by captial region boundary
  mutate(dtm = ymd_hms(date),
         Month = month(dtm,label = TRUE,abbr = FALSE),
         Year = year(dtm)) %>%
  unite(date_name,Month,Year, sep = " ",remove = FALSE) %>%
  mutate(date_name = factor(date_name, levels = c(
    "January 2016",
    "February 2016",
    "March 2016",
    "April 2016",
    "May 2016",
    "June 2016",
    "July 2016",
    "August 2016",
    "September 2016",
    "October 2016",
    "November 2016",
    "December 2016",
    "January 2017",
    "February 2017",
    "March 2017",
    "April 2017",
    "May 2017",
    "June 2017",
    "July 2017",
    "August 2017",
    "September 2017"))) %>%
  filter(!is.na(date_name))



### read strava road shapefile as simple feature
e <- read_sf(paste0(getwd(),"/input/exposure/STRAVA/Edges/victoria_osm_edges.shp")) %>%
  st_transform(crs = 26910) %>% #transform projection to NAD 83 UTM Zone10N
  filter(CLAZZ!=1) %>% #remove ferry routes
  filter(st_intersects(x=.,y=st_as_sfc(st_bbox(crd)),sparse=FALSE)) %>% # subset road network by capital region bounding box
  distinct(.keep_all = TRUE) #remove duplicate edges


### read in a-spatial strava counts aggregated by road geometry and month
#extract all relevant .csv filenames
files <- list.files(path = paste0(getwd(),"/input/exposure/STRAVA/Edges"),
                    pattern = "*.csv", full.names = TRUE)
#extract filenames with monthly totals for road links
monthly_total <- files[str_detect(files,"_month_") & str_detect(files,"total")]
#monthly filenames with weekday totals for road links
monthly_weekday <- files[str_detect(files,"_month_") & str_detect(files,"weekday")]
#monthly filenames weekday totals for road links
monthly_weekend <- files[str_detect(files,"_month_") & str_detect(files,"weekend")]
#read in each file representing monthly counts by road geometry id
m_total <- sapply(monthly_total, read.csv,simplify = FALSE) %>%
  bind_rows(.id = "id") %>% #bind dataframes together
  mutate(id = str_remove(id,
                         paste0(getwd(),"/input/exposure/STRAVA/Edges/victoria_201601_201709_ride_rollup_month_"))) %>%
  mutate(id = str_remove(id,
                         ".csv")) %>%
  separate(id, into = c("year","month","agg")) %>%
  unite("year_month",year,month,remove = FALSE) %>%
  mutate(month_rank = case_when(
    year_month  == "2016_1" ~ 1,
    year_month  == "2016_2" ~ 2,
    year_month  == "2016_3" ~ 3,
    year_month  == "2016_4" ~ 4,
    year_month  == "2016_5" ~ 5,
    year_month  == "2016_6" ~ 6,
    year_month  == "2016_7" ~ 7,
    year_month  == "2016_8" ~ 8,
    year_month  == "2016_9" ~ 9,
    year_month  == "2016_10" ~ 10,
    year_month  == "2016_11" ~ 11,
    year_month  == "2016_12" ~ 12,
    year_month  == "2017_1" ~ 13,
    year_month  == "2017_2" ~ 14,
    year_month  == "2017_3" ~ 15,
    year_month  == "2017_4" ~ 16,
    year_month  == "2017_5" ~ 17,
    year_month  == "2017_6" ~ 18,
    year_month  == "2017_7" ~ 19,
    year_month  == "2017_8" ~ 20,
    year_month  == "2017_9" ~ 21))
#create dataframe linking month integer to name
month_join <- data.frame(month = as.character(seq(1,12,by = 1)),
                         month_name = c("January",
                                        "February",
                                        "March",
                                        "April",
                                        "May",
                                        "June",
                                        "July",
                                        "August",
                                        "September",
                                        "October",
                                        "November",
                                        "December"))
#create dataframe where each row refers to a specific road link and a month and year within the range of strava data
m_total_all_months <- expand.grid(edge_id = e$ID,
                                  month_rank = seq(1,21,by = 1)) %>%
  mutate(year_month = case_when(
    month_rank == 1  ~ "2016_1" ,
    month_rank == 2 ~ "2016_2",
    month_rank == 3  ~ "2016_3",
    month_rank == 4 ~ "2016_4" ,
    month_rank  == 5 ~ "2016_5",
    month_rank  == 6 ~ "2016_6",
    month_rank  == 7 ~ "2016_7",
    month_rank  == 8 ~ "2016_8",
    month_rank  == 9 ~ "2016_9",
    month_rank  == 10 ~ "2016_10",
    month_rank  == 11 ~ "2016_11",
    month_rank  == 12 ~ "2016_12",
    month_rank  == 13 ~ "2017_1" ,
    month_rank  == 14 ~ "2017_2" ,
    month_rank  == 15 ~ "2017_3" ,
    month_rank  == 16 ~ "2017_4" ,
    month_rank  == 17 ~ "2017_5" ,
    month_rank  == 18 ~ "2017_6" ,
    month_rank  == 19 ~ "2017_7" ,
    month_rank  == 20 ~ "2017_8" ,
    month_rank  == 21 ~ "2017_9" )) %>%
  separate(year_month,into = c("year","month"),remove = FALSE) %>%
  as_tibble() %>%
  arrange(edge_id,month_rank)
#Join the data frame with rows for each possible road link and month total with the data.
#NAs indicate that the road link and month was not available within the strava dataset, suggesting there was no
#strava riders that month for that road segment. Replace NAs with 0
m_total_all_months <- m_total_all_months %>%
  left_join(m_total) %>%
  replace_na(list(tactcnt=0))#replace na's with 0s (e.g. edge not in the strava data)
#NOTE to self, expand NA replacement to other counts, not just total counts for that month

#################### 2. Aggregate Exposure/Crashes by Road Network Link ####################

###Join monthly counts to road geometry
e_m_total <- left_join(e,m_total_all_months,
                       by = c("ID"="edge_id"))
#add more month and year columns
e_m_total <- e_m_total %>%
  left_join(month_join) %>%
  unite(col = date,month_name,year,sep = " ",remove = FALSE) %>%
  mutate(date_name = factor(date, levels = c(
    "January 2016",
    "February 2016",
    "March 2016",
    "April 2016",
    "May 2016",
    "June 2016",
    "July 2016",
    "August 2016",
    "September 2016",
    "October 2016",
    "November 2016",
    "December 2016",
    "January 2017",
    "February 2017",
    "March 2017",
    "April 2017",
    "May 2017",
    "June 2017",
    "July 2017",
    "August 2017",
    "September 2017")))

### Conflate ICBC crashes with road geometry
#snap points to nearest road link (max distance =30m)
icbc_snp <- st_snap_points(icbc_sf,e,max_dist = 30)
icbc_snp <- st_sf(icbc_sf %>%
                    st_drop_geometry() %>%
                    mutate(geom=icbc_snp))
#join edge id to each point
icbc_eID <- st_join(icbc_snp,e["ID"],
                    join = st_is_within_distance, dist = 0.001)
#aggregate points by edge id
eID_sum_crashes <- icbc_eID %>%
  filter(!is.na(ID)) %>%
  group_by(ID,date_name) %>%
  summarise(n_crashes = sum(CRASHES)) %>%
  ungroup() %>%
  st_drop_geometry()
#join counts by edge id
e_m_total <- left_join(e_m_total,eID_sum_crashes ,
                       by = c("ID","date_name"))  %>%
  replace_na(list(n_crashes=0))

### Conflate bikemaps incidents with road geometry
bm_snp <- st_snap_points(bm,e,max_dist = 30)
bm_snp <- st_sf(bm %>%
                  st_drop_geometry() %>%
                  mutate(geom=bm_snp))
#join edge id to each point
bm_eID <- st_join(bm_snp,e["ID"],
                  join = st_is_within_distance, dist = 0.001)
#aggregate points by edge id
eID_sum_crashes <- bm_eID %>%
  filter(!is.na(ID)) %>%
  group_by(ID,date_name) %>%
  summarise(n_bm_crashes = sum(p_type=="collision"),
            n_nearmiss = sum(p_type=="nearmiss")) %>%
  ungroup() %>%
  st_drop_geometry()
#join counts by edge id
e_m_total <- left_join(e_m_total,eID_sum_crashes ,
                       by = c("ID","date_name"))  %>%
  replace_na(list(n_bm_crashes=0,n_nearmiss=0))

e_m_total <- e_m_total %>% mutate(
  all_crashes = n_crashes + n_bm_crashes,
  all_incidents = n_crashes + n_bm_crashes + n_nearmiss)

e_o_total <- e_m_total %>%
  group_by(ID,OSM_NAME,KM,KMH) %>%
  summarise(total_tactn = sum(tactcnt),
            total_icbc_crashes = sum(n_crashes),
            total_bm_crashes = sum(n_bm_crashes),
            total_nearmiss = sum(n_nearmiss),
            total_all_crashes = sum(all_crashes),
            total_all_incidents = sum(all_incidents))

e <- NULL

e_o_total <- e_o_total %>%
  mutate(avg_monthly_tactn = total_tactn/21,
         avg_annual_daily_tactn = total_tactn/638) %>%
  mutate(CAADB = 40 + avg_annual_daily_tactn*50)  #(Roy et al 2019/Jestico et al 2016)

#################### 3. Create Lixelized Network ####################

network_sln <- SpatialLinesNetwork(e_o_total)
network_sln <- sln_clean_graph(network_sln)# Remove unconnected roads

### Lixelize cleaned network data

lixel_list_10m <- lixelize_network(
  sf_network = network_sln@sl,
  max_lixel_length = 10
)


#################### 5. Write to Shapefile ####################

st_write(crd, paste0(getwd(), "/input/processed/","study_area.gpkg"))

st_write(e_o_total,paste0(getwd(), "/input/processed/","edge_ec_201601_201709_total.gpkg"))

st_write(e_m_total,paste0(getwd(), "/input/processed/","edge_ec_201601_201709_monthly.gpkg"))

st_write(bm_snp,paste0(getwd(), "/input/processed/","inc_bm_201601_201709_snp_30m.gpkg"))

st_write(icbc_snp,paste0(getwd(), "/input/processed/","inc_icbc_201601_201709_snp_30m.gpkg"))


