#### Long-tailed macaque activity pattern analysis across SEA ####


#By: Samuel Xin Tham Lee

#Load libraries
library(tidyverse)
library(overlap)
library(activity)
library(sp)
library(lubridate)

#Load all necessary dataframes
caps = as.data.frame(read.csv("SEA_Activity_dataset_for_community_guild-level_analyses_20230306.csv"))

#Only select for long-tailed macaques
macaque_caps <- caps %>% filter(caps$Species == "Macaca_fascicularis")

#Keep environment clean!
rm(caps)


####### Loop to split captures per status ####

#Create a FLII status vector
status = c('Degraded', 'Intact')

#store captures data frames for each status in this list
macaque.caps = list() 

# test out the loop with specific variables
i=1

for(i in 1:length(status)){ #repeat for each status
  
  # specify status and subset caps for it
  a = status[i]
  b = macaque_caps[macaque_caps$FLII_status == a,]
  
  #Save it!
  macaque.caps[[i]] = b
  names(macaque.caps)[i]= a
  
}

## keep environment clean
rm(a,b,i)

####### Loop to Generate Actmods per status for each status #####

## store actmod objects in this list 
actmod.list = list()

## variables to practice the loop with 
i=1

for(i in 1:length(macaque.caps)){ # repeat for each status (i.e. highest level of the nested list)
  
  #subset for a specific status
  a = macaque.caps[[i]]
  b = names(macaque.caps)[i] #save the name
    
    if(nrow(a) <= 100){
      
      act = fitact(a$time.rad, sample = "model", adj = 1, reps = 10000) #COME HERE AND CHAGNE to 1000 for final run
      
    }else{
      
      act = fitact(a$time.rad, sample = "data", adj = 1, reps = 10000) #COME HERE AND CHAGNE to 1000 for final run
    } #end conditional 
    
  ## Save it
  actmod.list[[i]] = act
  names(actmod.list)[i] = b
  
} # end loop per species

#keep environment clean
rm(i,a,b,act)

#Save it!
saveRDS(actmod.list, 'SEA_Activity_temporal_shift_long-tailed_macaque_actmods_20230209.RDS')
actmod.list = readRDS('SEA_Activity_temporal_shift_long-tailed_macaque_20230209.RDS')

#### Extract and determine the change in peak activity at guild-level ####

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
write.csv(empty, 'SEA_Activity_temporal_shift_peak_activity_long-tailed_macaque_20230117.csv', row.names = F)

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
  coord_cartesian(ylim = c(0, 0.8))+ ## Need to limit y axis to positive values due to geom_rug's jitter
  annotate("rect", xmin = 0, xmax = 1.17, ymin = 0, ymax= c(0,0.8), color = "#333333", alpha = .4)+
  annotate("rect", xmin = 1.18, xmax = 1.96, ymin = 0, ymax= c(0,0.8), color = "grey", alpha = .1)+
  annotate("rect", xmin = 4.32, xmax = 5.1, ymin = 0, ymax= c(0,0.8), color = "grey", alpha = .1)+
  annotate("rect", xmin = 5.11, xmax = 6.3, ymin = 0, ymax= c(0,0.8), color = "#333333", alpha = .4)+
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
