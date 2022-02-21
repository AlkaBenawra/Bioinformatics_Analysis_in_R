# Mus Musculus (Hous mouse) Analysis - BOLD Database 

#In this script we extract Mus Musculus data from the bold database and analyze the metadata and DNA barcode sequences.

#Load Libraries 
library(tidyverse)
library(vegan)

#Extract mus mysculus data from BOLD database (API call) - Extract 02/20/2022
#mus_musculus <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Mus&musculus&format=tsv")


#Write file to current directory
getwd()
#write_tsv(mus_musculus, "mus_musculus.tsv")

#Read file in from directory
#For data from BOLD, each row is an individual biological, so each row is a consistent type of sample.
mus_musculus <- read_tsv("/Users/alkabenawra/Documents/OneDrive/Job Search 2022/Technical Question Practice/mus_musculus.tsv")

#Explore dataset ---
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

#Which countries have the most barcode data?
count.by.country <- mus_musculus %>%
  group_by(country) %>%
  summarize(count = length(processid)) %>%
  arrange(desc(count))

#Another way to get counts
mus_musculus %>%
  group_by(country) %>%
  summarize(count = n()) 


###### INTRODUCTION TO BIODIVERSITY ANALYSIS USING VEGAN----

#For this part, we will use the BINs as a molecular proxy for species identities to explore the completeness of sampling in the genus Daphnia for the DNA barcoding campaign.Barcode Index Numbers (BINs) are a unique identifier assigned on the basis of clustering patterns in the DNA barcode sequence data. The BIN algorithm was calibrated against traditional species identifications for selected, well-studied animal groups. BINs are very helpful for biodiversity analysis, especially of invertebrates, where obtaining morphological species-level identifications is very difficult or impossible, due to the high diversity and many undescribed species.

#Grouping the data by BIN and counting the number of records in each BIN.
mus.by.bin <- mus_musculus %>%
  group_by(bin_uri) %>%
  summarise(count = n())

unique(mus.by.bin$bin_uri)

#Get data in comm data object format
mus.comm.object <- mus.by.bin %>%
  pivot_wider(names_from = bin_uri, values_from = count)

#Display rarefied species richness for community ecologists
rare_curve <- rarecurve(mus.comm.object)

#The curve shows the rate of species discovery is slowing for the house mouse, but if we continued sampling we would add more species. 

#As we add countries, do we add a lot of new BINs?
#Count of the number of specimens per BIN per country.
mus.by.country.bin <- mus_musculus %>%
  group_by(country, bin_uri) %>%
  filter(!is.na(country) & !is.na(bin_uri)) %>%
  summarise(count = n())

#Create comm object
mus.country.bin.comm.object <- mus.by.country.bin %>%
  pivot_wider(names_from = bin_uri, values_from = count)

#Turn BIN NA's to zeros
mus.country.bin.comm.object[is.na(mus.country.bin.comm.object)] <- 0
#specaccum(mus.country.bin.comm.object)

str(mus.country.bin.comm.object)
summary(mus.country.bin.comm.object)
#We need our data all to be in numerical format to have a correct comm object type to pass to the function specaccum.

#remove column with country names and put country names as row names
rownames(mus.country.bin.comm.object) 

#Remove row column 
mus.country.bin.comm.object <- mus.country.bin.comm.object %>%
  remove_rownames() %>%
  column_to_rownames(var = "country")

#Run a species accumulation curve analysis. Resampling is performed to see how BINs accumulate as sites (in this case countries) are added.
mus.accum <- specaccum(mus.country.bin.comm.object)

#Plot the model (The number of species for a certain number of sampled sites or individuals.)
plot(mus.accum)

#What is that plot telling us? As countries are sampled, additional species are added to the dataset. So, not all species of Mus musculus are found everywhere.
  


######PART 3 - Analysis-

# Create a subset as a new object retaining only records found in the northern (lat >0) AND western hemispheres (lon < 0)
mus.subset <- mus_musculus %>%
  filter(lat >0 & lon <0 & !is.na(lat) & !is.na(lon))


#First, checking dimensions. We should NOT have lost columns. column number should be the same in original and subset datasets. Yes, both have 80 variables (columns).
dim(mus_musculus)
dim(mus.subset)

#Next, checking distribution of latitude and longitude values. Yes, the min and max for each is as expected. Latitudes should be positive and between 0 and 90, and longitudes should be negative and between 0 and -180.
summary(mus.subset$lat)
summary(mus.subset$lon)


#How many unique BINs (bin_uri column) are found in United States in this dataset?
mus_musculus %>%
  group_by(country) %>%
  filter(country =="United States") %>%
  summarize(count = n())

table(mus_musculus$country)

#How many unique bins are there for ther united states
mus_musculus %>%
  filter(country == "United States") %>%
  summarise(count = length(unique(bin_uri)))


#How many records in this dataset are missing country information?
sum(is.na(mus_musculus$country))

#How many are missing a BIN assignment? 
sum(is.na(mus_musculus$bin_uri))

#We can also answer such questions by looking at an overall summary for all countries if we wanted to. This will just output the results to the screen. (Remember, we could assign the output to a data object, if we wanted to.)
mus_musculus %>%
  group_by(country) %>%
  summarize(count.unique.bins = length(unique(bin_uri))) %>%
  arrange(desc(count.unique.bins))


#Summarize the count of records by country in this subset.
table(mus_musculus$country)

#Sort the countries in decreasing order of records
sort(table(mus_musculus$country), decreasing = TRUE)

#Get the top three countries with the highest number of records
sort(table(mus_musculus$country), decreasing = TRUE)[1:3]



