# Mus Musculus (Hous mouse) Analysis - BOLD Database 

#In this script we extract Mus Musculus data from the bold database and analyze the metadata and DNA barcode sequences.

#Load Libraries 
library(tidyverse)

#Extract mus mysculus data from BOLD database (API call) - Extract 02/20/2022
#mus_musculus <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Mus&musculus&format=tsv")


#Write file to current directory
getwd()
#write_tsv(mus_musculus, "mus_musculus.tsv")

#Read file in from directory
#For data from BOLD, each row is an individual biological, so each row is a consistent type of sample.
mus_musculus <- read_tsv("/Users/alkabenawra/Documents/OneDrive/Job Search 2022/Technical Question Practice/mus_musculus.tsv")

#Explore dataset
summary(mus_musculus)
str(mus_musculus)
class(mus_musculus) #it's a tibble
names(mus_musculus)


#Subset dataframe
mus_musculus <- mus_musculus[, c(1, 8, 10, 12, 14, 16, 18, 20, 22, 25, 26, 32, 47, 48, 55:59, 70, 72)]

#Exploratory Data Analysis ----

#Lets look at range of latituddes (the North pole is at 90 degrees (North), and the South pole is at -90 )
hist(mus_musculus$lat)
hist(mus_musculus[[13]])

class(mus_musculus[13])
class(mus_musculus[[13]])

#Lets's look at longitude range (-180 to 180)
hist(mus_musculus$lon)
summary(mus_musculus$lon)


#Find mean latitude and longitude
mean(mus_musculus$lon, na.rm = TRUE)
mean(mus_musculus$lat, na.rm = TRUE)

#What sample is collect furthest north?
mus_musculus$processid[which.max(mus_musculus$lat)]

#What sample is collected furthest south?
mus_musculus$processid[which.min(mus_musculus$lat)]

#What is the west most sample (-)?
mus_musculus$processid[which.min(mus_musculus$lon)]

#What is the east most sample (+)?
mus_musculus$processid[which.max(mus_musculus$lon)]

#Find max longitude and latitude
max(mus_musculus$lon, na.rm = TRUE)
min(mus_musculus$lat, na.rm = TRUE)

#How many records are collected north of the 10 degrees N? 10
sum(mus_musculus$lat >10, na.rm = TRUE)

#How many taxonomic species names are in the dataset?
length(unique(mus_musculus$species_name))

#How many BINs are in the sample?
length(unique(mus_musculus$bin_uri))

#What is the ratio of BINs to species? What does the result mean?
length(unique(mus_musculus$bin_uri))/length(unique(mus_musculus$species_name))

#How many records have a BIN number?
?grep
length(grep(":", mus_musculus$bin_uri))


#Filter out records with a bin
unique(mus_musculus$bin_uri)
mus_musculus2 <- mus_musculus %>%
  filter(is.na(bin_uri) == FALSE)

#Can filter out with grep too
mus_musculus2 <- mus_musculus[grep(":", c(mus_musculus$bin_uri)), ]

#Verified bins with na were remove
unique(mus_musculus2$bin_uri)
sum(is.na(mus_musculus2$bin_uri))
sum(is.na(mus_musculus$bin_uri)) #in original dataset there was 65 records with NA as bin


#Create a new dataframe called hundred_mus, which contains rows 1 to 100 from mus_musculus.
hundred_mus <- mus_musculus[c(1:100), ]

#Create a new dataframe called TestingIndexing, from original dataframe mus_musculus, which contains rows 50 to 100 and columns 10 to 15.
TestingIndexing <- mus_musculus[c(50:100), c(10:15)]

#Using dataframe Hundred_mus for the next tasks...How many unique collectors collected these samples?
length(unique(hundred_mus$collectors))

#Does every specimen have a unique processid?
nrow(hundred_mus)
length(unique(hundred_mus$processid))

#Find records with non-unique process ids
subset(hundred_mus,  duplicated(hundred_mus$processid))

#Count the number of records mined by the university of Toronto
length(mus_musculus[grep("Toronto", mus_musculus$collectors), ])


#unique countries 
length(unique(hundred_mus$country))


#Which country has the most samples among this 100 Daphnia? 
which.max(table(hundred_mus$country))


#How many records are missing geographic coordinate data?
sum(is.na(hundred_mus$lat))
sum(is.na(hundred_mus$lon))


#How many records have geographic coordinate data?
sum(!is.na(hundred_mus$lat))
sum(!is.na(hundred_mus$lon))

#How many records have a BIN assignment?
sum(!is.na(hundred_mus$bin_uri))

