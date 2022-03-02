#Biostrings Package

#Load Library
library(tidyverse)
library(ape)
library(RSQLite)

#Bioconductor packages

#Code for installation for R version 3.5 or greater. You can install whichever Bioconductor packages you need.
#Instructions from: https://bioconductor.org/install/
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install(c("Biostrings", "DECIPHER", "muscle"))

library(DECIPHER)
library(muscle)
library(Biostrings)

#Analyzing DNA sequence dataset from a group of parasitoid wasps (genus Cotesia), containing sequences from a single marker gene (COI and 28S).

#viewing the quick overview vignette
vignette("BiostringsQuickOverview", package="Biostrings")

#Viewing the Pairwise Alignments vignette
vignette("PairwiseAlignments", package="Biostrings")

#Viewing Multiple Alignment vignette
vignette("MultipleAlignments", package="Biostrings")

#Load data file for Cotesia, a genus of parasitoid wasps, BOLD API call (03/02/2022).
#Cotesia <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Cotesia&format=tsv")

#write file to local directory ]
#getwd()
#write_tsv(Cotesia, "Cotesia.tvs")

#Read cortesia file in 
Cotesia <- read_tsv("Cotesia.tvs")

#Explore data set
summary(Cotesia)
str(Cotesia)
names(Cotesia)
class(Cotesia) #It is a tibble

#Genes in dataset
unique(Cotesia$markercode)

#The genes in the dataset include COI-5P, ITS, 28S, Wnt1, 16S, CYTB, ND3, ND4L, COII, ND4, ND6, ND1, ND5-0, COXIII, COI-3P.

#Show how many sequences are associated with each markercode in the dataset. 
sort(table(Cotesia$markercode), decreasing = TRUE )


count.by.marker. <- Cotesia %>%
  group_by(markercode) %>%
  filter(is.na(nucleotides) == FALSE) %>%
  summarize(n = length(processid)) %>%
  arrange(desc(n)) %>%
  print()

count.by.marker. <- Cotesia %>%
  group_by(markercode) %>%
  filter(is.na(nucleotides) == FALSE) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>%
  print()

#We can see that most of the data are COI-5P (the 5' end of the cytochrome c oxidase subunit I mitochondrial gene), which is the standard marker gene for DNA barcoding of animals. There are 82 records having markercode 28s (the nuclear large ribosomal RNA subunit gene). The other genes all have a small sample size. Note that this is to be expected for BOLD data, due to the focus on DNA barcoding. Other markers are generally regarded as "complementary" or "additional" to the main barcode gene region. If you want to analyze a larger variety, it would be wise to use NCBI databases instead as a resource. We will cover accessing NCBI resources in a future lesson.


#Filter Cotesia dataset for the gene COI-5P and only include nucleotides with ATCG
Cotesia_COI5P <- Cotesia %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect(nucleotides, "[ATCG]"))

         
#Verify dataset only has COI-5P gene
unique(Cotesia_COI5P$markercode)


#Create dataset with only 238
Cotesia_28S <- Cotesia %>%
  filter(markercode == "28S" & str_detect(nucleotides, "[ATCG]"))

#Verify dataset only has 28S data
unique(Cotesia_28S$markercode)

#Remove original dataset to clear up workspace
rm(Cotesia)


#Analyze DNA Sequences ----

#Look at datatype for DNA squence
class(Cotesia_COI5P$nucleotides) #character

#To use biostrings package we need to convert the nucleotide sequence to a DNAStringSet. A DNAStringSet contains character (string) data but belongs to a stricter object class.

#Change nucleotide class to DNAStringSet
DNAStringSet(Cotesia_COI5P$nucleotides)
Cotesia_COI5P$nucleotides <- as.vector(x)
class(Cotesia_COI5P$nucleotides)

Cotesia_28S$nucleotides <- as.vector(DNAStringSet(Cotesia_28S$nucleotides))


#This below line will calculate the nucleotide counts of the specified letter, in this case A.
Cotesia.COI.AFreq <- as.data.frame(letterFrequency(DNAStringSet(Cotesia_COI5P$nucleotides), letters = "A"))

#This next line will calculate the nucleotide frequency of the specified letters, A or T. So, this calculates the frequency at which an A or T is observed across positions. Have a look at the new object.
Cotesia.COI.ATFreq <- as.data.frame(letterFrequency(DNAStringSet(Cotesia_COI5P$nucleotides), letters = "AT"))

#This calculates the frequency of each nucleotide separately. Again, be sure to have a look in the viewer to follow what this is doing.
Cotesia.COI.NucFreq <- as.data.frame(letterFrequency(DNAStringSet(Cotesia_COI5P$nucleotides), letters = c("A", "C", "G", "T", "N")))

Cotesia.28S.NucFreq <- as.data.frame(letterFrequency(DNAStringSet(Cotesia_28S$nucleotides), letters = c("A", "C", "G", "T", "N")))

#We could calculate the AT frequency (as a proportion of the total sequence, excluding Ns). Here, we are using the function mutate() to add a column onto the end of Cotesia.COI.NucFreq
Cotesia.COI.NucFreq <- Cotesia.COI.NucFreq %>%
  mutate(ATproportion = ((A + T) / (A + T + G + C)))

#Let's look at the histogram of AT proportions. (This is an exploratory data analysis. I often use base functions list hist() and plot() for exploratory analysis. For reports or for publication, I would recommend to use other tools, such as through package ggplot.)
hist(Cotesia.COI.NucFreq$ATproportion)

#Mean AT proportion
mean(Cotesia.COI.NucFreq$ATproportion)


#Get ATproportion for 28S
names(Cotesia.28S.NucFreq)
Cotesia.28S.NucFreq <- Cotesia.28S.NucFreq %>%
  mutate(ATproportion = ((A + T) / (A+T+C+G+N)))

#viewing histogram. We will come back to the unusual value in a future class
hist(Cotesia.28S.NucFreq$ATproportion)
mean(Cotesia.28S.NucFreq$ATproportion)

#Are the mean AT frequencies significantly different between COI and 28S? Answer: Yes.
t.test(Cotesia.28S.NucFreq$ATproportion, Cotesia.COI.NucFreq$ATproportion)

#Calculating k-mers. Below, we are calculating the counts of short oligonucleotide sequences of length 4. Counts of all possible 4-nucleotide strings are tabulated. k-mers can be a very useful sequence feature. See papers I posted to the CourseLink on uses of k-mers for sequence classification if you are interested in more depth. We could choose a different length for our k-mers.
oligo4 <- oligonucleotideFrequency(DNAStringSet(Cotesia_COI5P$nucleotides), 4)
length(oligo4)

#Have a look at oligo4 in the viewer to understand better what this is doing. Then, for comparison, let's try a different number and then have a look at the object. Example:
oligo3 <- oligonucleotideFrequency(DNAStringSet(Cotesia_COI5P$nucleotides), 3)

#But... now, let's return to that unusual 28S value...

#What is the minimum value of ATproportion in 28S?
min(Cotesia.28S.NucFreq$ATproportion)

#Which sample exhibits the minimum value? this will give row number.
which.min(Cotesia.28S.NucFreq$ATproportion)

#Or, we could have used indexing to obtain the processid (sample identifier) of the sample with the minimum AT proportion in one step.
Cotesia_28S$processid[which.min(Cotesia.28S.NucFreq$ATproportion)]

#Checking that there was just one very low value. What is the next min value after leaving out row 65? Yes, there was just a single extremely low outlier, as the next-lowest sequence is closer to the mean AT proportion.
min(Cotesia.28S.NucFreq$ATproportion[-65])

#Converting that specific sequence to a character vector
to_blast <- as.character(Cotesia_28S$nucleotides[65])

#checking class
class(to_blast)

#printing to screen
to_blast

#We will review BLAST another day and as part of your Assignment #2. In brief, we can use the NCBI implementation of the BLAST algorithm to compare our sequence against publicly available sequences. For now, we can use the web tool. In a future lesson, we will go over how to access NCBI tools from R. Link to BLAST sequence against the nucleotide sequence database:
#https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

#Our unusual Cotesia sequence is possibly a misidentification. BLASTEd to one 100% match (likely itself) but then mainly to other Hymenoptera of other taxonomic identifications. This sequence could be dropped or explored further prior to analysis, depending upon the goals of a study. In such a case, I would recommend to code in the data drop into the script, with an explanation.

######PART 6 - CODING CHALLENGES----

#This section contains commenting only. I would recommend to try to solve these tasks yourself and then you can compare your answers against the posted EXAMPLE answers. As usual in R, there is more than one way to solve these tasks.

#I would suggest to use the broom buttom to wipe clean your workspace prior to commencing this section.

#Read in the Cotesia dataset again.
Cotesia <- read_tsv("Cotesia.tvs")

#Then check data parsing and data summary.
summary(Cotesia)
str(Cotesia)
class(Cotesia)

#How many unique countries are represented in this dataset? Answer: 75
length(unique(Cotesia$country))

#Which country is the best represented in this dataset? Canada, has 2030 records
sort(table(Cotesia$country), decreasing = TRUE)[1]

which.max(table(Cotesia$country))
table(Cotesia$country)[11]


#Find the nucleotide composition of the Wnt1 gene and calculate the AT frequency.
Cotesia.Wnt1 <- Cotesia %>%
  filter(markercode == "Wnt1" & !is.na(nucleotides))

unique(Cotesia.Wnt1$markercode)


#Build a histogram of the AT frequency values for the Wnt1 gene. What shape is the distribution? Is there evidence of outliers?

#First, converting class for nucleotide sequences and Finding nucleotide frequencies
Cotesia.Wnt1.NucFreq <- as.data.frame(letterFrequency(DNAStringSet(Cotesia.Wnt1$nucleotides), letters = c("A", "T", "G", "C", "N")))

#Adding a column to the same dataframe specifically with the AT frequency.
Cotesia.Wnt1.NucFreq <- Cotesia.Wnt1.NucFreq %>%
  mutate(ATproportion = ((A +T) / (A +T + G + C + N)))

#Histogram. Intereseting: quite a narrow range of AT frequencies in this gene.
hist(Cotesia.Wnt1.NucFreq$ATproportion)

