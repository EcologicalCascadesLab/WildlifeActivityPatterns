#### Human activity pattern analysis across SEA ####


#By: Samuel Xin Tham Lee

#Load libraries
library(tidyverse)
library(overlap)
library(activity)
library(sp)
library(lubridate)

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

#Here, we will create a separate captures dataset that redefines FLII. 
#For this dataset, we will define FLII status as "degraded" when it is < Q1 and "intact" when it is > Q3
human_caps_extreme = human_caps

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
human_caps_extreme <- merge(human_caps_extreme, covs, by = c('camera_id', 'survey_id')) 

#Determine FLII status by median
human_caps$FLII_status = "Degraded"
human_caps$FLII_status[human_caps$forest_integrity > median(human_caps$forest_integrity)] = "Intact"

#Determine FLII sttaus by quantiles
human_caps_extreme$FLII_status <- "NA"
human_caps_extreme$FLII_status[human_caps_extreme$forest_integrity <= quantile(human_caps_extreme$forest_integrity,0.25)] = "Degraded"
human_caps_extreme$FLII_status[human_caps_extreme$forest_integrity >= quantile(human_caps_extreme$forest_integrity,0.75)] = "Intact"
human_caps_extreme <- human_caps_extreme %>% filter(FLII_status %in% c("Degraded", "Intact"))

#Keep environment clean!
rm(covs)

####### Loop to split captures per status ####

#Create a FLII status vector
status = c('Degraded', 'Intact')

#store captures data frames for each status in this list
human.caps = list() 

# test out the loop with specific variables
i=1

for(i in 1:length(status)){ #repeat for each status
  
  # specify status and subset caps for it
  a = status[i]
  b = human_caps[human_caps$FLII_status == a,]
  
  #Save it!
  human.caps[[i]] = b
  names(human.caps)[i]= a
  
}

## keep environment clean
rm(a,b,i)

####### Loop to split captures per status (for FLII split by quantiles) ####

#store captures data frames for each status in this list
human.caps.extreme = list() 

# test out the loop with specific variables
i=1

for(i in 1:length(status)){ #repeat for each status
  
  # specify status and subset caps for it
  a = status[i]
  b = human_caps_extreme[human_caps_extreme$FLII_status == a,]
  
  #Save it!
  human.caps.extreme[[i]] = b
  names(human.caps.extreme)[i]= a
  
}

## keep environment clean
rm(a,b,i)

####### Loop to Generate Actmods per status for each status #####

## store actmod objects in this list 
actmod.list = list()

## variables to practice the loop with 
i=1

for(i in 1:length(human.caps)){ # repeat for each status (i.e. highest level of the nested list)
  
  #subset for a specific status
  a = human.caps[[i]]
  b = names(human.caps)[i] #save the name
    
    if(nrow(a) <= 100){
      
      act = fitact(a$time.rad, sample = "model", adj = 1, reps = 10000) #COME HERE AND CHAGNE to 1000 for final run
      
    }else{
      
      act = fitact(a$time.rad, sample = "data", adj = 1, reps = 10000) #COME HERE AND CHAGNE to 1000 for final run
    } #end conditional 
    
  ## Save it
  actmod.list[[i]] = act
  names(actmod.list)[i] = b
  
} 

#keep environment clean
rm(i,a,b,act)

#Save it!
#saveRDS(actmod.list, 'SEA_Activity_temporal_shift_human_actmods_20230209.RDS')
#actmod.list = readRDS('SEA_Activity_temporal_shift_human_actmods_20230209.RDS')

####### Loop to Generate Actmods per status for each status (for FLII status split by quantiles) #####

## store actmod objects in this list 
actmod.list.extreme = list()

## variables to practice the loop with 
i=1

for(i in 1:length(human.caps.extreme)){ # repeat for each status (i.e. highest level of the nested list)
  
  #subset for a specific status
  a = human.caps.extreme[[i]]
  b = names(human.caps.extreme)[i] #save the name
  
  if(nrow(a) <= 100){
    
    act = fitact(a$time.rad, sample = "model", adj = 1, reps = 10000) #COME HERE AND CHAGNE to 1000 for final run
    
  }else{
    
    act = fitact(a$time.rad, sample = "data", adj = 1, reps = 10000) #COME HERE AND CHAGNE to 1000 for final run
  } #end conditional 
  
  ## Save it
  actmod.list.extreme[[i]] = act
  names(actmod.list.extreme)[i] = b
  
} 

#keep environment clean
rm(i,a,b,act)

#Save it!
saveRDS(actmod.list.extreme, 'SEA_Activity_human_actmods_extreme_20230529.RDS')


## Calculate p-value

human_actmods <- list(actmod.list, actmod.list.extreme) #Create a nested list 
names(human_actmods) <- c("actmod.list", "actmod.list.extreme")

#Create an empty list to save p-values later
p_value.list <- list()

i=1

for (i in 1:length(human_actmods)) { #repeat for each treatment 
  
  a = human_actmods[[i]] #Select a single treatment
  b = names(human_actmods)[i] #Save the name of the treatment
  
  #Compare activity patterns
  c = as.data.frame(compareCkern(a[[1]], a[[2]]))
  
  #Add a treatment and outputs column
  c$treatment = b
  c$outputs = rownames(c)
  
  #Save it 
  p_value.list[[i]] = c
  names(p_value.list)[i] = b
  
}

#Keep environment clean!
rm(a,b,c,i)

#Extract outputs and covert it into a dataframe
df = do.call(rbind, p_value.list)
rownames(df)= NULL

#Filter out the p-values and rename the p-value column
df <- df %>% 
  filter(outputs == "pNull") %>% 
  rename('pNull' = 'compareCkern(a[[1]], a[[2]])') %>% 
  select(-outputs)

##Calculate overlap

#Create an empty list to store results
overlap.list = list()

i=1;l=1

for (i in 1:length(human_actmods)) { #repeat for each treatment
  
  a = human_actmods[[i]]
  b = names(human_actmods)[i]
  
  c = min(length(a[[1]]@data), length(a[[2]]@data))
  
  if(c <= 75){ #calculate overlap estimate
    
    overlap.coef = overlapEst(a[[1]]@data, a[[2]]@data, adjust = 0.8, type = "Dhat1")
    
  }else{
    
    overlap.coef = overlapEst(a[[1]]@data, a[[2]]@data, adjust = 1, type = "Dhat4")
    
  }
  
  
  if(c <= 75){ 
    
    overlap.boot = bootstrap(a[[1]]@data, a[[2]]@data, nb = 10000, adjust = 0.8, type = "Dhat1")
    
    
  }else{
    
    overlap.boot = bootstrap(a[[1]]@data, a[[2]]@data, nb = 10000, adjust = 1, type = "Dhat4")
    
    
  }
  
  overlap.CI = bootCI(overlap.coef, overlap.boot, conf = 0.95) #calculate 95% CI of overlap estimate
  
  d = as.data.frame(overlap.CI)
  
  d = d["norm0",]
  
  e = as.data.frame(overlap.coef)
  
  f = cbind(e,d)
  
  f$treatment = b
  
  overlap.list[[i]] = f
  names(overlap.list)[i] = b
  
}


#Keep environment clean!
rm(a,b,c,d,e,f,i,l,overlap.boot,overlap.CI,overlap.coef)

#Merge it with df
a = do.call(rbind, overlap.list)
df <- merge(df, a, by = 'treatment')
rm(a)

write.csv(df, "SEA_Activity_human_overlaps_p-values_20230601.csv", row.names = F)


#### Extract and determine the change in peak activity ####

#Create dataframe to fill
colnamestoextract = c("status", "peak_activity_radians")
empty = as.data.frame(matrix(NA,nrow = 2,ncol = 2))
colnames(empty) = colnamestoextract

# loop to pull peak activity value

i=1

for(i in 1:length(actmod.list)){
  
  a = actmod.list[[i]] #Subset a term from the list called act
  
  b = names(actmod.list)[i] #Save the names of the term
  
  empty$status[i] = b #Save the names of the term
  
  c = a@pdf #Extract the probability density function values
  
  d = as.data.frame(c) #Convert it into a dataframe
  
  e = filter(d, y == max(y)) # Remove other rows and now only contains maximum density values 
  
  empty$peak_activity_radians[i] = e$x # X will be the peak activity where most of our species detections (density) are found
  
}

#Keep environment clean
rm(a,b,c,d,e,i,colnamestoextract)

# Convert peak_activity in radians to seconds
empty <- empty %>% mutate(peak_activity_seconds = (peak_activity_radians*86400)/(2*pi))

#Convert seconds into time format
str(empty)
empty$peak_activity_seconds <- as.integer(empty$peak_activity_seconds)
empty <- empty %>% mutate(peak_activity = seconds_to_period(peak_activity_seconds))

#Select columns that we need
empty <- select(empty, status, peak_activity)

#Save it!
write.csv(empty, 'SEA_Activity_temporal_shift_peak_activity_human_20230117.csv', row.names = F)

#### Extract and determine the change in peak activity (for FLII split by quantiles) ####

#Create dataframe to fill
colnamestoextract = c("status", "peak_activity_radians")
empty = as.data.frame(matrix(NA,nrow = 2,ncol = 2))
colnames(empty) = colnamestoextract

# loop to pull peak activity value

i=1

for(i in 1:length(actmod.list.extreme)){
  
  a = actmod.list.extreme[[i]] #Subset a term from the list called act
  
  b = names(actmod.list.extreme)[i] #Save the names of the term
  
  empty$status[i] = b #Save the names of the term
  
  c = a@pdf #Extract the probability density function values
  
  d = as.data.frame(c) #Convert it into a dataframe
  
  e = filter(d, y == max(y)) # Remove other rows and now only contains maximum density values 
  
  empty$peak_activity_radians[i] = e$x # X will be the peak activity where most of our species detections (density) are found
  
}

#Keep environment clean
rm(a,b,c,d,e,i,colnamestoextract)

# Convert peak_activity in radians to seconds
empty <- empty %>% mutate(peak_activity_seconds = (peak_activity_radians*86400)/(2*pi))

#Convert seconds into time format
str(empty)
empty$peak_activity_seconds <- as.integer(empty$peak_activity_seconds)
empty <- empty %>% mutate(peak_activity = seconds_to_period(peak_activity_seconds))

#Select columns that we need
empty <- select(empty, status, peak_activity)

#Save it!
write.csv(empty, 'SEA_Activity_peak_activity_human_extreme_20230529.csv', row.names = F)

#Keep environment clean!
rm(empty, human.caps.extreme)

#### Loop to prepare actmods for visualisation ####

pdf_cov = list()

i=1

for (i in 1:length(actmod.list)) {#repeat for each condition
  
  #Select a single condition
  a = actmod.list[[i]]
  b = names(actmod.list)[i]
  
  # #Create an empty dataframe store multiple activity distributions here
  # #by matching the format of the pdf 
  # c = data.frame(a[[1]]@pdf)
  # c = c[0,]
  
  # extract the probability density function 
    c = data.frame(a@pdf)
    c$status = b
    
    #save this!
    #c = rbind(c, f)
  
  #last save
  pdf_cov[[i]] = c
  names(pdf_cov)[i] = b
}

#Keep environment clean!
rm(a,b,c,i)

#Visualise activity using ggplot

pdf <- do.call(rbind,pdf_cov)

plot =
  ggplot(pdf, aes(x=x, y=y, fill = status))+
  #geom_line(size = 1)+
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.6)+
  #geom_rug(data = df[df$FLII_status == "Degraded",], aes(x = Time.rad, y=0), position = "jitter", sides = "b", color = "#4169e1", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
  #geom_rug(data = df[df$FLII_status == "Intact",], aes(x = Time.rad, y=0), position = "jitter", sides = "b", color = "#FC6600", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
  scale_x_continuous(breaks = c(0, 3.14, 6.3), 
                     labels = c("Midnight", "Noon", "Midnight"))+
  coord_cartesian(ylim = c(0, 1.0))+ ## Need to limit y axis to positive values due to geom_rug's jitter
  annotate("rect", xmin = 0, xmax = 1.17, ymin = 0, ymax= c(0,1.0), color = NA, alpha = .4)+
  annotate("rect", xmin = 1.18, xmax = 1.96, ymin = 0, ymax= c(0,1.0), color = NA, alpha = .1)+
  annotate("rect", xmin = 4.32, xmax = 5.1, ymin = 0, ymax= c(0,1.0), color = NA, alpha = .1)+
  annotate("rect", xmin = 5.11, xmax = 6.3, ymin = 0, ymax= c(0,1.0), color = NA, alpha = .4)+
  #geom_vline(xintercept = 1.57, colour = '#ffc808', alpha = 0.5, size = 1) +
  #geom_vline(xintercept = 4.71, colour = '#ffc808', alpha = 0.5, size = 1) +
  theme_test()+
  #labs(x = "\nTime of day", y = "Activity\n")+
  scale_fill_manual(name = pdf$status, values = c("#FC6600", "#4169e1"))+
  theme(axis.title = element_blank(), axis.text = element_text(size = 18),
        legend.title = element_blank(), legend.text = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")

plot

#### Loop to prepare actmods (FLII split by quantiles) for visualisation ####

pdf_cov_extreme = list()

i=1

for (i in 1:length(actmod.list.extreme)) {#repeat for each condition
  
  #Select a single condition
  a = actmod.list.extreme[[i]]
  b = names(actmod.list.extreme)[i]
  
  # #Create an empty dataframe store multiple activity distributions here
  # #by matching the format of the pdf 
  # c = data.frame(a[[1]]@pdf)
  # c = c[0,]
  
  # extract the probability density function 
  c = data.frame(a@pdf)
  c$status = b
  
  #save this!
  #c = rbind(c, f)
  
  #last save
  pdf_cov_extreme[[i]] = c
  names(pdf_cov_extreme)[i] = b
}

#Keep environment clean!
rm(a,b,c,i)

#Visualise activity using ggplot

pdf_extreme <- do.call(rbind,pdf_cov_extreme)

plot =
  ggplot(pdf_extreme, aes(x=x, y=y, fill = status))+
  #geom_line(size = 1)+
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.6)+
  #geom_rug(data = df[df$FLII_status == "Degraded",], aes(x = Time.rad, y=0), position = "jitter", sides = "b", color = "#4169e1", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
  #geom_rug(data = df[df$FLII_status == "Intact",], aes(x = Time.rad, y=0), position = "jitter", sides = "b", color = "#FC6600", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
  scale_x_continuous(breaks = c(0, 3.14, 6.3), 
                     labels = c("Midnight", "Noon", "Midnight"))+
  coord_cartesian(ylim = c(0, 1.0))+ ## Need to limit y axis to positive values due to geom_rug's jitter
  annotate("rect", xmin = 0, xmax = 1.17, ymin = 0, ymax= c(0,1.0), color = NA, alpha = .4)+
  annotate("rect", xmin = 1.18, xmax = 1.96, ymin = 0, ymax= c(0,1.0), color = NA, alpha = .1)+
  annotate("rect", xmin = 4.32, xmax = 5.1, ymin = 0, ymax= c(0,1.0), color = NA, alpha = .1)+
  annotate("rect", xmin = 5.11, xmax = 6.3, ymin = 0, ymax= c(0,1.0), color = NA, alpha = .4)+
  #geom_vline(xintercept = 1.57, colour = '#ffc808', alpha = 0.5, size = 1) +
  #geom_vline(xintercept = 4.71, colour = '#ffc808', alpha = 0.5, size = 1) +
  theme_test()+
  #labs(x = "\nTime of day", y = "Activity\n")+
  scale_fill_manual(name = pdf_extreme$status, values = c("#FC6600", "#4169e1"))+
  theme(axis.title = element_blank(), axis.text = element_text(size = 18),
        legend.title = element_blank(), legend.text = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")

plot


