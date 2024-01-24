###################################################################################
### Calculating the percentage of diurnal, nocturnal and crepuscular detections ###
###################################################################################

#Load relevant packages
library(tidyverse);library(lubridate);library(plyr);library(overlap);library(activity);library(sp)

#Load all necessary dataframes
caps = as.data.frame(read.csv("ECL captures summer spp projects_20210322.csv"))

#Remove Danum2020 surveys (Covariates not calculated for this survey)
caps = caps %>% 
  filter(!survey_id == "Danum_2020")

#Only select for humans
sort(unique(caps$Species))
human_caps <- caps %>% filter(caps$Species %in% c('Homo_sapiens', 'Homo_sapiens_hiker', 'Homo_sapiens_Poacher',
                                                  'Homo_sapiens_ranger', 'Homo_sapiens_researcher', 'Homo_sapiens_tourist'))

##########################################################################
#### Organizing camera coordinates to account for daylength variation ####
##########################################################################

#Load in coordinates of camera locations
ECL_metadata = as.data.frame(read.csv("ECL and Collaborator Camera Trap Metadata_20220802.csv"))
names(ECL_metadata)
sort(unique(ECL_metadata$Landscape))
ECL_metadata = ECL_metadata %>% 
  select(camera_id, survey_id, Landscape, Longitude, Latitude) %>% 
  filter(Landscape %in% c('Bukit_Barisan_Selatan_National_Park', 'Danum_Valley_Conservation_Area',
                          'Gunung_Leuser_National_Park', 'Khao_Yai_National_Park', 'Lambir_Hills_National_Park',
                          'Singapore', 'Ulu_Muda_Forest_Reserve', 'Pasoh_Forest_Reserve', 'Kerinci_Seblat_National_Park')) #Missing Khao Chong survey!

#Check if all cameras in caps are accounted for
setdiff(human_caps$camera_id, ECL_metadata$camera_id) # Also missing 2 cameras from lambir hills (LM198_CAM35 and LM201-CAM68)

#Load in metadata from Calebe
ECL_metadata_2 = as.data.frame(read.csv("Full_Metadata_by_Camera.csv"))
names(ECL_metadata_2)
ECL_metadata_2 = ECL_metadata_2 %>% 
  select(camera_id, survey_id, region, country, Y_lat, X_long) %>% 
  filter(survey_id == 'KhaoChong2018.ECL')

ECL_metadata_2 = ECL_metadata_2 %>% 
  mutate(survey_id = str_remove_all(survey_id, ".ECL")) %>% 
  mutate(camera_id = str_remove_all(camera_id, "KhaoChong2018.ECL.")) %>% 
  select(camera_id, survey_id, Y_lat, X_long)

#Remove landscape column so that we can rbind both metadata together
ECL_metadata = ECL_metadata %>% 
  select(-Landscape)

#Change column names of ECL_metadata_2
colnames(ECL_metadata_2) = c("camera_id", "survey_id", "Latitude", "Longitude")

#rbind both metadatas
ECL_metadata <- rbind(ECL_metadata, ECL_metadata_2)
setdiff(human_caps$camera_id, ECL_metadata$camera_id) # Missing Cameras: C13_sorted, SR5_sorted, SR9_sorted, LM198_CAM35 and LM201-CAM68

#Load in last remaining coordinates given by Jon and combined it with ECL_metadata
ECL_metadata_3 <- as.data.frame(read.csv("ECL_missing_coordinates_20221103.csv"))
ECL_metadata <- rbind(ECL_metadata, ECL_metadata_3)
human_caps <- human_caps %>% filter(!camera_id == 'SR4_sorted') #Don't have the coordinates for this camera!

#Remove survey_id from ECL_metadata as they are not organised the same
ECL_metadata <- ECL_metadata %>% 
  select(-survey_id)

#Merge captures and camera coordinates into one single dataset
human_caps <- merge(human_caps, ECL_metadata, by = 'camera_id')

#Keep environment clean!
rm(ECL_metadata, ECL_metadata_2, ECL_metadata_3)

#####################################################################
#### Accounting for daylength variation using camera coordinates ####
#####################################################################

#Here, we will be using camera coordinates and Suntime () function to take into account of daylength variation.
#Since SEA is near the equator, it would not vary that much but its good practice to do this!
#The radian time produced is a relative value of both sunset and sunrise hours.
#For instance, if a certain camera location has a sunrise time of 0630 hr, the Suntime () function will regard it as 1.57 radians or pie/2 and vice versa for sunset time (3pie/2)

#Since we are dealing with multiple surveys, there will be multiple time zones and therefore we need to split the dataset based on respective time zones

sort(unique(human_caps$survey_id))

#### Thailand surveys ####

thai <- human_caps %>% 
  filter(survey_id %in% c('KhaoYai2019', 'KhaoChong2018'))

#create a vector containing time in radians which will be used in circular analyses later on
time <- gettime(thai$Photo.Time, format = "%H:%M:%S", scale = c("radian")) 

#Create a vector containing date in the correct format
date <- as.POSIXct(thai$Photo.Date, tz= "Asia/Bangkok", format = '%d/%m/%Y')

# Create a SpatialPoints object with the location
coords <- data.frame(thai$Longitude,thai$Latitude)
coords2 <- sp::SpatialPoints(coords, proj4string = sp::CRS("+epsg=4087 +proj=longlat +datum=WGS84"))

#Correct for sunrise and sunset time based on lat and long
st <- sunTime(time, date, coords2)

#Merge it with main thailand dataset
thai$time.rad <- st

#Keep environment clean!
rm(coords,coords2,date,time,st)

#### Malaysian and Singapore surveys ####

my_sg <- human_caps %>% 
  filter(survey_id %in% c('Danum_Valley_2019a', 'Danum2018', 'Lambir2017', 'Pasoh_TEAM_2013', 'Pasoh_TEAM_2014',
                          'Pasoh_TEAM_2015', 'Pasoh_TEAM_2017', 'Singapore', 'Ulu_Muda_2015a', 'Ulu_Muda_2015b',
                          'Ulu_Muda_2015c', 'Ulu_Muda_2015d', 'Ulu_Muda_2016a', 'Ulu_Muda_2016b', 'Ulu_Muda_2016c'))

#create a vector containing time in radians which will be used in circular analyses later on
time <- gettime(my_sg$Photo.Time, format = "%H:%M:%S", scale = c("radian")) 

#Create a vector containing date in the correct format
date <- as.POSIXct(my_sg$Photo.Date, tz= "Asia/Singapore", format = '%d/%m/%Y')

# Create a SpatialPoints object with the location
coords <- data.frame(my_sg$Longitude,my_sg$Latitude)
coords2 <- sp::SpatialPoints(coords, proj4string = sp::CRS("+epsg=4087 +proj=longlat +datum=WGS84"))

#Correct for sunrise and sunset time based on lat and long
st <- sunTime(time, date, coords2)

#Merge it with main peninsular and singapore dataset
my_sg$time.rad <- st

#Keep environment clean!
rm(coords,coords2,date,time,st)

#### Sumatran surveys ####

sum <- human_caps %>% 
  filter(survey_id %in% c('BBS', 'Kerinci', 'Leuser'))

#create a vector containing time in radians which will be used in circular analyses later on
time <- gettime(sum$Photo.Time, format = "%H:%M:%S", scale = c("radian")) 

#Create a vector containing date in the correct format
date <- as.POSIXct(sum$Photo.Date, tz= "Asia/Jakarta", format = '%d/%m/%Y')

# Create a SpatialPoints object with the location
coords <- data.frame(sum$Longitude,sum$Latitude)
coords2 <- sp::SpatialPoints(coords, proj4string = sp::CRS("+epsg=4087 +proj=longlat +datum=WGS84"))

#Correct for sunrise and sunset time based on lat and long
st <- sunTime(time, date, coords2)

#Merge it with main sumatra dataset
sum$time.rad <- st

#Keep environment clean!
rm(coords,coords2,date,time,st)

#Merge all datasets together!
human_caps <- rbind(my_sg, thai)
human_caps <- rbind(human_caps, sum)

#Keep environment clean!
rm(my_sg, sum, thai)

############################################
#### Merge captures with FLII covariate ####
############################################

#Load in covariate dataset
covs <- as.data.frame(read.csv("ECL_metadata_cam_level_summer_spp_20210322.csv"))
names(covs)
head(covs)

#Filter out other covariates and only include FLII
covs <- covs %>% 
  select(camera_id, survey_id, forest_integrity)

#Check if all cameras are accounted for
setdiff(human_caps$camera_id, covs$camera_id) #All the cameras in captures are accounted for so all good!

#Merge both captures and FLII
human_caps <- merge(human_caps, covs, by = c('camera_id', 'survey_id'))

#Determine FLII sttaus by quantiles
human_caps$FLII_status <- "NA"
human_caps$FLII_status[human_caps$forest_integrity <= quantile(human_caps$forest_integrity,0.25)] = "Degraded"
human_caps$FLII_status[human_caps$forest_integrity >= quantile(human_caps$forest_integrity,0.75)] = "Intact"
human_caps <- human_caps %>% filter(FLII_status %in% c("Degraded", "Intact"))

#Keep environment clean!
rm(covs)

##Convert the radians back to normal time format

##First, let's define our diel niche:

# 1. Nocturnal hours (19:31 pm to 04:29 am)(5.109451 rad to 1.173734 rad)
# 2. Crepuscular hours (04:30 am to 07:30 am or 16:30 pm to 19:30 pm)(1.178097 rad to 1.963495 or 4.319690 to 5.105088)
# 3. Diurnal hours (07:31 am to 16:29 pm)(1.967859 rad to 4.315327 rad)

#Trim the dataset to only include columns that we need
human_caps <- human_caps %>% select("Species", "FLII_status", "time.rad")

#Include the respective time of day categories
human_caps$time_of_day = "NA"
human_caps$time_of_day[human_caps$time.rad >= 5.109451 | human_caps$time.rad <= 1.178096] = "night"
human_caps$time_of_day[human_caps$time.rad >= 1.967859 & human_caps$time.rad <= 4.319689] = "day"
human_caps$time_of_day[human_caps$time_of_day == "NA"] = "crepuscular"
sort(unique(human_caps$time_of_day))

#Calculate total day, crepuscular and night detections
df = human_caps %>% 
  group_by(FLII_status) %>% 
  dplyr::count(time_of_day) %>% 
  pivot_wider(names_from = "time_of_day", values_from = "n") 

#Save it!
write.csv(df, "SEA_Activity_raw_human_detections_extreme_20230602.csv", row.names = F)

perc <- df %>% 
  summarise(perc_change_day = (day[FLII_status == "Degraded"]/total_detections[FLII_status == "Degraded"]) - 
              (day[FLII_status == "Intact"]/total_detections[FLII_status == "Intact"]),
            perc_change_crep = (crepuscular[FLII_status == "Degraded"]/total_detections[FLII_status == "Degraded"]) - 
              (crepuscular[FLII_status == "Intact"]/total_detections[FLII_status == "Intact"]),
            perc_change_night = (night[FLII_status == "Degraded"]/total_detections[FLII_status == "Degraded"]) - 
              (night[FLII_status == "Intact"]/total_detections[FLII_status == "Intact"])*100)


#Save it!
write.csv(perc, "SEA_Activity_perc_change_humans_extreme_20230602.csv", row.names = F)





