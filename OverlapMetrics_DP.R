# Load required packages
library(tidyverse)

###### Overlap Metrics #####

## Area Overlap
## for binary data
## measures proportion of an area where two species co-occur
area_overlapfn <- function(prey, pred, area){
  total_area <- sum(area, na.rm = T)
  sum(area[pred > 0 & prey > 0], na.rm = T)/total_area
}

## Range Overlap
## for binary data
## measures the proportion of one species range where the other co-occurs
range_overlapfn<-function(prey, pred, area){
  area_prey <- sum(area[prey > 0], na.rm = T)
  sum(area[pred > 0 & prey > 0], na.rm = T)/area_prey
}

## Schoener's D
## density or probability of occurrence data
## measures how equally predator and prey share available resources
schoeners_overlapfn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  1 - 0.5 * (sum(abs(p_prey-p_pred), na.rm = T))
}

## Bhattacharyya's coefficient
## density or probability of occurrence data
## measures whether two species use space independently
bhatta_coeffn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum(sqrt(p_prey*p_pred), na.rm = T)
}

###### Applying overlap metrics to actual data ######

# load data and calculate overlap metrics for each day
sdm2015 <- readRDS("~/Documents/UCSC/Dissertation/Sanctuary Futures/SanctuaryFutures/data/Processed Data/Year/humpanch_2015.rds") # load a specific year of SDM data
daily2015zoom <- sdm2015 %>% # make new dataframe called dailyzoom2015 which is a copy of sdm2015 and do the following below
  mutate(Humpback_Core = ifelse(Humpback_HS >= 0.28,1,0),Anchovy_Core = ifelse(Anchovy_HS >= 0.5,1,0),Area=1) %>% # create a binary data column for range and area overlap of each species
  filter(Lat>=34 & Lat < 35,Lon <= -119 & Lon > -120) %>% # retain values only within this lat-lon range
  group_by(Date) %>% # organize by date
  summarise(AO=area_overlapfn(prey=Anchovy_Core,pred=Humpback_Core,area=Area),RO=range_overlapfn(prey=Anchovy_Core,pred=Humpback_Core,area=Area), # apply the overlap metrics and output as new dataframe
            Schoener=schoeners_overlapfn(prey=Anchovy_HS,pred=Humpback_HS), Bhatty=bhatta_coeffn(prey=Anchovy_HS,pred=Humpback_HS)) %>% 
  pivot_longer(cols=c(AO,RO,Schoener,Bhatty),names_to="Metric_Name",values_to="Overlap_Metric") # change the dataframe into long format for easier plotting
