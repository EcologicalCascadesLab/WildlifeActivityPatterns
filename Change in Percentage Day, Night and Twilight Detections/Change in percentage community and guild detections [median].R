###################################################################################
### Calculating the percentage of diurnal, nocturnal and crepuscular detections ###
###################################################################################

#Load relevant packages
library(tidyverse);library(lubridate);library(plyr)

#Set working directory
setwd("C:/Users/AMD/Dropbox/Sam Lee Honours/SEA Activity Data Analysis Sam Honours/SEA Activity Temporal Shift Analysis Sam Honours/SEA Activity percentage change in detections community-level Sam Honours")

#Load the dataset
caps <- as.data.frame(read.csv("SEA_Activity_dataset_for_community_guild-level_analyses_20230306.csv"))

#Check the datasets
head(caps);names(caps);anyNA(caps)

##Convert the radians back to normal time format

##First, let's define our diel niche:

# 1. Nocturnal hours (19:31 pm to 04:29 am)(5.109451 rad to 1.173734 rad)
# 2. Crepuscular hours (04:30 am to 07:30 am or 16:30 pm to 19:30 pm)(1.178097 rad to 1.963495 or 4.319690 to 5.105088)
# 3. Diurnal hours (07:31 am to 16:29 pm)(1.967859 rad to 4.315327 rad)

#Trim the dataset to only include columns that we need
caps <- caps %>% select("Species", "FLII_status", "time.rad", "trophic_guild")

#Include the respective time of day categories
caps$time_of_day = "NA"
caps$time_of_day[caps$time.rad >= 5.109451 | caps$time.rad <= 1.178096] = "night"
caps$time_of_day[caps$time.rad >= 1.967859 & caps$time.rad <= 4.319689] = "day"
caps$time_of_day[caps$time_of_day == "NA"] = "crepuscular"
sort(unique(caps$time_of_day))

#####################
## Community-level ##
#####################

##Create three separate datasets that corresponds to each species composition: all species, pigs and macaques excluded and only pigs and macaques included.

#Create all_species dataset
all_species = caps

#Create no_pigs_macaques dataset
no_pigs_macaques = caps %>% filter(!Species %in% c("Sus_barbatus", "Sus_scrofa", "Macaca_fascicularis", "Macaca_nemestrina"))

#Create no_pigs_macaques dataset
only_pigs_macaques = caps %>% filter(Species %in% c("Sus_barbatus", "Sus_scrofa", "Macaca_fascicularis", "Macaca_nemestrina"))

#Keep environment clean!
rm(caps)

## Here, we will create a list so that we can loop it!
act_list = list(all_species, no_pigs_macaques, only_pigs_macaques)
names(act_list) <- c("all_species", "no_pigs_macaques", "only_pigs_macaques")
rm(all_species, no_pigs_macaques, only_pigs_macaques) #Keep environment clean!

## Loop to calculate the raw total detections of each species composition

total_sp_list <- list()

i=1

for (i in 1:length(act_list)) {#repeat for each species dataset
  
  #Select for a single dataset
  a = act_list[[i]]
  b = names(act_list)[i]
  
  caps = a %>% dplyr::count(FLII_status, name = "total_detections") #Calculate total detections in both intact and degraded forests
  
  c = a %>% 
    group_by(FLII_status) %>% 
    dplyr::count(time_of_day) %>% 
    pivot_wider(names_from = "time_of_day", values_from = "n") #Calculate total day, crepuscular and night detections
  
  caps = merge(caps, c, by = "FLII_status") #Merge both datasets together
  
  caps$species_composition = b #Add species composition column to clarify which community is being analysed
  
  #Save it
  total_sp_list[[i]] = caps
  names(total_sp_list)[i] = b
}

#Keep environment clean!
rm(a,b,c,caps,i)

#Rbind them together
df = do.call(rbind, total_sp_list)

#Save it!
write.csv(df, "SEA_Activity_raw_total_detections_20230531.csv", row.names = F)

perc <- df %>% 
  group_by(species_composition) %>% 
  summarise(perc_change_day = (day[FLII_status == "Degraded"]/total_detections[FLII_status == "Degraded"]) - 
              (day[FLII_status == "Intact"]/total_detections[FLII_status == "Intact"]),
            perc_change_crep = (crepuscular[FLII_status == "Degraded"]/total_detections[FLII_status == "Degraded"]) - 
              (crepuscular[FLII_status == "Intact"]/total_detections[FLII_status == "Intact"]),
            perc_change_night = (night[FLII_status == "Degraded"]/total_detections[FLII_status == "Degraded"]) - 
              (night[FLII_status == "Intact"]/total_detections[FLII_status == "Intact"])*100)


perc$species_composition <- c("all_species", "no_pigs_macaques", "only_pigs_macaques")

write.csv(perc, "SEA_Activity_changes_in_perc_detections_20230531.csv", row.names = F)

#################
## Guild-level ##
#################

#Loop to split captures into guilds and save it within a nested list
guilds <- sort(unique(caps$trophic_guild))

guild_caps <- list()

i=1

for (i in 1:length(guilds)) {#repeat for each guild
  
  a = guilds[[i]]
  
  c = caps[caps$trophic_guild == a,]
  
  guild_caps[[i]] = c
  names(guild_caps)[i] = a
}

#Keep environment clean!
rm(a,b,c,i)

## Loop to calculate the raw total detections of each guild

total_guild <- list()

i=1

for (i in 1:length(guild_caps)) {#repeat for each species dataset
  
  #Select for a single dataset
  a = guild_caps[[i]]
  b = names(guild_caps)[i]
  
  df = a %>% dplyr::count(FLII_status, name = "total_detections") #Calculate total detections in both intact and degraded forests
  
  c = a %>% 
    group_by(FLII_status) %>% 
    dplyr::count(time_of_day) %>% 
    pivot_wider(names_from = "time_of_day", values_from = "n") #Calculate total day, crepuscular and night detections
  
  df = merge(df, c, by = "FLII_status") #Merge both datasets together
  
  df$guild = b #Add guild column to clarify which guild is being analysed
  
  #Save it
  total_guild[[i]] = df
  names(total_guild)[i] = b
}

#Keep environment clean!
rm(a,b,c,df,i)

#Rbind them together
df = do.call(rbind, total_guild)

#save it
write.csv(df, "SEA_Activity_raw_guild_detections_20230605.csv", row.names = F)

#Calculate percentage
perc_guild <- df %>% 
  group_by(guild) %>% 
  summarise(perc_change_day = (day[FLII_status == "Degraded"]/total_detections[FLII_status == "Degraded"]) - 
              (day[FLII_status == "Intact"]/total_detections[FLII_status == "Intact"]),
            perc_change_crep = (crepuscular[FLII_status == "Degraded"]/total_detections[FLII_status == "Degraded"]) - 
              (crepuscular[FLII_status == "Intact"]/total_detections[FLII_status == "Intact"]),
            perc_change_night = (night[FLII_status == "Degraded"]/total_detections[FLII_status == "Degraded"]) - 
              (night[FLII_status == "Intact"]/total_detections[FLII_status == "Intact"])*100)

perc_guild$guild <- guilds

#save this!
write.csv(perc_guild, "SEA_Activity_guild_perc_change_20230605.csv", row.names = F)

#Keep environment clean!
rm(perc_guild,guild_caps,total_guild,df)
