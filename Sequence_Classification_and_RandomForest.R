#Sequence Classification and Random Forest

#Load libraries
library(tidyverse)
#install.packages("randomForest")
library(randomForest)
library(BiocManager)
library(Biostrings)

#We will use sequence features to identify sequences to one of two genes and to identify sequences to taxonomic groups.


#Ran the below line at March 3rd, 2022, at 8:15pm. Building upon our "Cotesia" test dataset (which is from the wasp family Braconidae) by updating our download for September 2019 and also downloading a second genus within a different family of wasps (family Ichneumonidae). Both of the genera being analyzed belong to the order Hymenoptera (which consists of wasps, bees, ants, and allies).

#wasps <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Cotesia|Campoletis&format=tsv")

#Write dataset to local director
#write_tsv(wasps, "wasps.tsv")   

#Read in dataset 
wasps <- read_tsv("wasps.tsv")


#Explore dataset
summary(wasps)
str(wasps)
names(wasps)
class(wasps)


#Show the number of records for each marker code in the dataset
table(wasps$markercode) #COI-5P and 28S have the most records in the dataset.

#Look at sequence in order to consider making a subset.

sort(table(wasps$markercode), decreasing = TRUE)

#histogram and summary of sequence lengths for COI. We can see a wide range of sequence lengths.
hist(nchar(wasps$nucleotides[wasps$markercode == "COI-5P"], keepNA = FALSE))

summary(nchar(wasps$nucleotides[wasps$markercode == "COI-5P"], keepNA = FALSE))


#histogram and summary of sequence lengths for  28S. We can see a wide range of sequence lengths.
hist(nchar(wasps$nucleotides[wasps$markercode == "28S"], keepNA = FALSE))

summary(nchar(wasps$nucleotides[wasps$markercode == "28S"], keepNA = FALSE))


#Make a new dataset with COI and 28S records only. Filter out records with NA for nucleotide sequence. 

wasps_filtered <- wasps %>%
  filter(!is.na(nucleotides)) %>%
  filter(markercode == "COI-5P"| markercode == "28S") 

#Total number of records with seuqence
length(wasps_filtered$nucleotides)

#Unique markercodes in filtered dataset
unique(wasps_filtered$markercode)
table(wasps_filtered$markercode)

#Checking there are no NA's in new dataset
sum(is.na(wasps_filtered$nucleotides)) # 0 records with NA in nucleotide column

#We are next going to create another dataset with a filter for sequence length as well. We will filter the dataset to retain records having a nucleotide  sequence (i.e. not NA) and with the marker being either COI-5P (this means the 5 primer end of the  COI gene) or 28S. Note the usage of a single "|" to indicate "OR", element-wise. We will also impose sequence length filters to constrain the range of variability within genes. Variability can be due to different primer sets being used, or to poor quality sequences. 28S, being a non-coding rRNA gene, would be expected to have length variability due to indels (insertions and deletions). Constraining the length variability would give us more consistent k-mer counts. I used the inter-quartile range from the above summaries as a guide. Specific choices for such filtering steps would vary depending upon the goals of the analysis. And, exploring results across a range of parameters can also be helpful. I named the dataframe wasps.length to remind us that this dataset is filtered for length. We may choose to use this length-filtered dataset for some analyses below. Note the usage of "&" to mean "and" and "|" to mean "or".

summary(nchar(wasps_filtered$nucleotides[wasps_filtered$markercode == "COI-5P"]))

summary(nchar(wasps_filtered$nucleotides[wasps_filtered$markercode == "28S"]))


wasps.length <- wasps %>%
  filter(!is.na(nucleotides)) %>%
  filter((markercode == "COI-5P" & 
           nchar(nucleotides) > 635 &
           nchar(nucleotides) < 658) |
           (markercode == "28S" & 
           nchar(nucleotides) > 353 & 
           nchar(nucleotides) < 536))

#Check new dataset to make sure everything was filtered as expected
unique(wasps.length$markercode)
sum(is.na(wasps.length$nucleotides))

summary(nchar(wasps.length$nucleotides[wasps.length$markercode == "COI-5P"]))

summary(nchar(wasps.length$nucleotides[wasps.length$markercode == "28S"]))


#Calculating sequence features ----

#Convert nucleotide sequence to DNAstringset so we can analyze it using biostrings package
wasps_filtered$nucleotides <- DNAStringSet(wasps_filtered$nucleotides)

#Calculating the nucleotide frequencies (similar to what we did in class 6 script) and appending onto our wasps1 dataframe using cbind().
wasps_filtered_updated <- cbind(wasps_filtered, 
as.data.frame(letterFrequency(DNAStringSet(wasps_filtered$nucleotides), letters = c("A", "T", "G", "C", "N"))))

#Adding A, T, G, and C proportions  as these may be helpful features less influenced by sequence length variation than are counts.

wasps_filtered_updated$Aprop <- wasps_filtered_updated$A / sum(wasps_filtered_updated$A + wasps_filtered_updated$T + wasps_filtered_updated$C + wasps_filtered_updated$G)

wasps_filtered_updated$Tprop <- wasps_filtered_updated$T / sum(wasps_filtered_updated$A + wasps_filtered_updated$T + wasps_filtered_updated$C + wasps_filtered_updated$G)

wasps_filtered_updated$Gprop <- wasps_filtered_updated$G / sum(wasps_filtered_updated$A + wasps_filtered_updated$T + wasps_filtered_updated$C + wasps_filtered_updated$G)

wasps_filtered_updated$Cprop <- wasps_filtered_updated$C / sum(wasps_filtered_updated$A + wasps_filtered_updated$T + wasps_filtered_updated$C + wasps_filtered_updated$G)


#Adding dinucleotide frequency (k-mers of length 2)
wasps_filtered_updated <- cbind(wasps_filtered_updated, as.data.frame(dinucleotideFrequency(DNAStringSet(wasps_filtered_updated$nucleotides))))

#Adding trinucleotide frequency (k-mers of length 3)
wasps_filtered_updated <- cbind(wasps_filtered_updated, as.data.frame(trinucleotideFrequency(DNAStringSet(wasps_filtered_updated$nucleotides))))

#If you want to add k-mers of a different lenght use oligonucleotideFrequency
#wasps_filtered_updated <- cbind(wasps_filtered_updated, as.data.frame(oligonucleotideFrequency(DNAStringSet(wasps_filtered_updated$nucleotides), x =  4)))


#Training Classification Model 1: Different Genes ----

#First, we are going to give the random forest algorithm a relative "easy" classification problem. We will ask whether the COI gene and 28S genes can be accurately classified from simple sequence features. As we saw in an earlier script, there was a large difference in nucleotide composition between these two genes in Cotesia, and therefore we SHOULD be able to classify sequences to one of these two genes on the basis of their sequence feature


#We will confine this test to the Cotesia genus, which has both genes represented.
wasps.Cotesia <- wasps_filtered_updated %>%
  filter(genus_name == "Cotesia")

#See count by markercode
table(wasps.Cotesia$markercode)

#As our maximum sample size for these two genes is 82 for the 28S gene in wasps.Cotesia, we will below sample 20 individuals (about 25% of the total for the 28S gene) from each gene to serve as the validation data set. These samples will be entirely set aside during the model training process. These samples should not be used until validating the model later. sample_n() is a very useful function for sampling a specific number of rows. Here, we are sampling 20 rows from each group (i.e. each marker code). See the resulting object to confirm.
set.seed(217)

#testing dataset
validation.data <- wasps.Cotesia %>%
  group_by(markercode) %>%
  sample_n(20)
#sample_n samples 20 rows from each group

#Now, we are creating a training dataset that does NOT overlap with the validation set. To do this, we will first remove processids that were among the samples randomly selected for the validation dataset. We can do this by asking for processid's that are not (!) in the validation set using %in%. Second, we will pick 62 individuals of each gene to serve in the training set. Note that this is a small sample size; this script is an EXAMPLE.

training.data <- wasps.Cotesia %>%
  filter(!(processid %in% validation.data$processid)) %>%
  group_by(markercode) %>%
  sample_n(62)

names(training.data)

#Next, we are building a classifier to separate the COI and 28S genes in these datasets, using the A, T, and G proportion as predictors. The response variable is markercode; we are trying to predict which gene a sequence belongs to, on the basis of simple sequence features alone.

gene_classifer <- randomForest(x = training.data[ , 81:83], y = as.factor(training.data$markercode), ntree = 1000, importance = TRUE)

#Perfect performance!
gene_classifer



#We can also dive into the specifics, such as by looking at the following. See the documentation for the randomForest() function for more information.
gene_classifer$importance
gene_classifer$oob.times

#Make model again with just two features
gene_classifier2 <- randomForest(x = training.data[ , 81:82], y = as.factor(training.data$markercode), ntree = 100, importance = TRUE)

#Still has perfect performance
gene_classifier2


#NOTE: This was a simple EXAMPLE, where we are only trying to distinguish TWO genes, which have very different sequence properties, as we saw in our earlier Biostrings script.

# Training Classification model #2: Genus ----

#Create classifer that can assign COI sequence to the correct wasp genus

#Restrict our analysis to the COI gene and sequences of similar length.

#Checking class of nucleotides
class(wasps_filtered$nucleotides)

#filtering wasps.filter down to just COI-5P
wasp.COI <- wasps_filtered_updated %>%
  filter(markercode == "COI-5P")

#Checking amount of data by genus
table(wasp.COI$genus_name)
#1368 records for Campoletis and 4873 records for 4873

#Making separate validation dataset, with 100 records per genus.
set.seed(123)
validation.data.coi <- wasp.COI %>%
  group_by(genus_name) %>%
  sample_n(100)


#Now, creating training dataset for COI, with 1100 records per genus, non-overlapping with the validation dataset. We could use more data, but I will sample down a bit for this exercise.
set.seed(43)
training.data.coi <- wasp.COI %>%
  filter(!processid %in% validation.data.coi$processid) %>%
  group_by(genus_name) %>%
  sample_n(1100)

#Now training the classifier. Note we would want rows to be included in multiple times, so we want the number of trees to be larger than the number of data points. For this example, you can reduce the number of trees if you have a slow computer.

genus_classifier <- randomForest(x = training.data.coi[ , 81:83], y = as.factor(training.data.coi$genus_name), ntree = 10000, importance = TRUE)

genus_classifier 
#I am actually very surprised how well this worked. Let's try just 2 features.

genus_classifier2 <- randomForest(x = training.data.coi[ , 81:82], y = as.factor(training.data.coi$genus_name), ntree = 10000, importance = TRUE)

genus_classifier2

#For fun, let's try more features, including both nucleotide proportions and dinucleotide frequencies.
genus_classifier3 <- randomForest(x = training.data.coi[ , 81:105], y = as.factor(training.data.coi$genus_name), ntree = 10000, importance = TRUE)

genus_classifier3

#Let's assess the robustness of our models by using cross-validation and see how much of a reduction in performance we get by reducing features. Below, we are using 10-fold cross-validation. The data are divided into 10 equal portions. Each time 9 portions are used for training the model and one portion to test the accurate the of the model.

cross.validate <- rfcv(trainx = training.data.coi [ , 81:83],trainy = as.factor(training.data.coi$genus_name), cv.fold = 10)

cross.validate

with(cross.validate, plot(n.var, error.cv, log="x", type="o", lwd=2))
#We can see that error rates are very low for this example. Error rates are even quite low with just a single predictor.

#Now let's return to our independent validation data to test the ability of our model to make predictions. We are using our predictive model genus_classifier to predict the groups (i.e. to predict genus name) for our data on the basis of three features (nucleotide proportions).

prediction <- predict(genus_classifier, validation.data.coi[ , 81:83])


######PART 6 - CODING CHALLENGE----

#Working in your groups, come up with your own classification problem. You may choose an example that is similar to the one above, or you could choose your own using a dataset of your choosing. Remember that R has a wide variety of built-in data sets that you can also explore.


# We want to see how sepal and petal characteristics of plants can be used to predict the species type using the iris dataset 

#Type this to see list of available datasets
data("iris")

summary(iris)
str(iris)
names(iris)
class(iris)
dim(iris)

#Number of records per a species
table(iris$Species)

#Check what unique species are in dataset 
unique(iris$Species)

#Checking sepal length variability in different plant species
summary(iris$Sepal.Length[iris$Species == "setosa"])

summary(iris$Sepal.Length[iris$Species == "versicolor"])

summary(iris$Sepal.Length[iris$Species == "virginica"])

#Create validation data set 
#I will use 10 records from each group
validation.data.iris <- na.omit(iris)

#No columns have NA for any records
dim(iris)
dim(validation.data.iris)

set.seed(123)
validation.data.iris <- validation.data.iris %>%
  group_by(Species) %>%
  sample_n(10)

#Now, creating training dataset for , with 40 records per species, non-overlapping with the validation dataset. We could use more data, but I will sample down a bit for this exercise.

#Make sure to filter out data used in validation set
row.names(validation.data.iris)

training.data.iris <- iris %>%
  filter(!row.names(iris) %in% validation.data.iris) %>%
  group_by(Species) %>%
  sample_n(40)


#filter(!row.names(iris) %in% validation.data$processid) %>%; you are filter out samples with rownames in the validation data set


#Next, we are building a classifier to separate the plant species names in these datasets, using sepal length, sepal width, petal length, and petal with predictors. The response variable is spccies; we are trying to predict which species the characteristics belongs to, on the basis of simple plant features alone.

iris_species_classifier <- randomForest(x = training.data.iris[ , 1:4], y = as.factor(training.data.iris$Species), ntree = 100, importance = TRUE)

iris_species_classifier


##Now, let's check out observed data against our predictions for the validation datset.

iris_classifier_predictions <- predict(iris_species_classifier, validation.data.iris[ , 1:4])


table(observed = validation.data.iris$Species, predicted = iris_classifier_predictions)
#The model had 100% accuracy

#So, we have very strong prediction performance for our "unseen" data for this model!! That is what we want. We don't want the model only to be fitted to the training data but to have found a generality that can be applied to unseen data.

#Although this case was more challenging than separating COI and 28S, this is still a relatively "easy" classification case. For this example, I chose two genera that are in different families of wasps. What would happen if we are trying to separate closely related species or genera? Also, it is important to note that there can be errors (e.g. due to misidentification or contamination) in public sequence databases. Often, the solution to overcome such errors is to use a lot of data, because it can be difficult to detect all such errors.

