#R_ Assignment

### Tanner M. Cook
### Ph.D. Student at Iowa State University
#### Copyright privileges are owned by Tanner Cook and Iowa State University. For use please contact Tanner Cook (tmcook@iastate.edu)


## This file contains code and descriptions from the R assignment in EEOB 546X.

## Load relevant libraries, setting directory, read files 
library("tidyverse")
library("tidyr")
library("ggplot2")

install.packages("reshape")
library("reshape")

setwd("C:/Users/Tanner/Desktop/Tanner/Iowa_State_Work/PhD_2020_2021/Bioinformatics_for_Biological_data/EEOB 546 R_Assignment")
Fang_genotypes= read_tsv("fang_et_al_genotypes.txt")
SNP_position= read_tsv("snp_position.txt")


## Inspecting data
head(Fang_genotypes)
head (SNP_position)

tail(Fang_genotypes)
tail(SNP_position)

dim(Fang_genotypes)
dim (SNP_position)

nrow(Fang_genotypes)
nrow(SNP_position)

str(Fang_genotypes)
str(SNP_position)


view(Fang_genotypes)
view(SNP_position)

### From this I learned that Fang_et_al_genotypes has a header row with Sample ID, Location information, and Group information 
  ### along with SNP ID labels. By transforming the data we are able to produce a table  ### that has three headers.
  
# Maize

## Isolate ZMMLR, ZMMIL,ZMMMR containing Rows and then Transpose Fang_et_al_Genotypes 
ZMM_Grouped_Fang_geno = filter(Fang_genotypes, grepl("ZMM", Group))
View(ZMM_Grouped_Fang_geno)

Transposed_ZMM_Grouped_Fang_geno= t(ZMM_Grouped_Fang_geno)

## Remove the three rows of headers in Genotype data except for sample ID
Headerless_Transposed_ZMM_Grouped_Fang_geno= Transposed_ZMM_Grouped_Fang_geno[-c(1,2,3),]
view(Headerless_Transposed_ZMM_Grouped_Fang_geno)


## Join SNP and Genotype Data

Headerless_Transposed_ZMM_Grouped_Fang_geno.df <- as.data.frame(Headerless_Transposed_ZMM_Grouped_Fang_geno)

row_nm_tocol_Headerless= rownames_to_column(Headerless_Transposed_ZMM_Grouped_Fang_geno.df, 'SNP_ID')
### Since the SNP ID's were set as a row name, we had to convert it to a column, but to do so 
### we had to make it a data frame

Join_SNP_Geno= right_join(SNP_position,row_nm_tocol_Headerless, by ='SNP_ID')

##Remove columns to ensure SNP_ID", the second column is "Chromosome", the third column is "Position", 
##and subsequent columns are genotype data from either maize or teosinte individuals.
  ### Assumption that genotype data means nucleotides only

Cleaned_Join_SNP_Geno = subset(Join_SNP_Geno, select = -c(cdv_marker_id, alt_pos:count_gene))


## Chromosome Files for MAize
###Since these already have ? marks for missing data, we do not have to impute here

Chr1_Maize_Cleaned_Join_SNP_Geno= filter(Cleaned_Join_SNP_Geno, Chromosome == 1)
Chr2_Maize_Cleaned_Join_SNP_Geno= filter(Cleaned_Join_SNP_Geno, Chromosome == 2)
Chr3_Maize_Cleaned_Join_SNP_Geno= filter(Cleaned_Join_SNP_Geno, Chromosome == 3)
Chr4_Maize_Cleaned_Join_SNP_Geno= filter(Cleaned_Join_SNP_Geno, Chromosome == 4)
Chr5_Maize_Cleaned_Join_SNP_Geno= filter(Cleaned_Join_SNP_Geno, Chromosome == 5)
Chr6_Maize_Cleaned_Join_SNP_Geno= filter(Cleaned_Join_SNP_Geno, Chromosome == 6)
Chr7_Maize_Cleaned_Join_SNP_Geno= filter(Cleaned_Join_SNP_Geno, Chromosome == 7)
Chr8_Maize_Cleaned_Join_SNP_Geno= filter(Cleaned_Join_SNP_Geno, Chromosome == 8)
Chr9_Maize_Cleaned_Join_SNP_Geno= filter(Cleaned_Join_SNP_Geno, Chromosome == 9)
Chr10_Maize_Cleaned_Join_SNP_Geno= filter(Cleaned_Join_SNP_Geno, Chromosome == 10)




# Teosinte

## Isolate ZMPBA, ZMPIL,ZMPJA containing Rows and then Transpose Fang_et_al_Genotype

ZMP_Grouped_Fang_geno = filter(Fang_genotypes, grepl("ZMP", Group))
View(ZMP_Grouped_Fang_geno)

Transposed_ZMP_Grouped_Fang_geno= t(ZMP_Grouped_Fang_geno)

## Remove the three rows of headers in Genotype data except for sample ID
Headerless_Transposed_ZMP_Grouped_Fang_geno= Transposed_ZMP_Grouped_Fang_geno[-c(1,2,3),]
view(Headerless_Transposed_ZMP_Grouped_Fang_geno)


##Join SNP and Genotype Data

Headerless_Transposed_ZMP_Grouped_Fang_geno.df <- as.data.frame(Headerless_Transposed_ZMP_Grouped_Fang_geno)

row_ZMP_nm_tocol_Headerless= rownames_to_column(Headerless_Transposed_ZMP_Grouped_Fang_geno.df, 'SNP_ID')
### Since the SNP ID's were set as a row name, we hd to convert it to a column, but to do so 
### we had to make it a data frame

Join_ZMP_SNP_Geno= right_join(SNP_position,row_ZMP_nm_tocol_Headerless, by ='SNP_ID')

##Remove columns to ensure SNP_ID", the second column is "Chromosome", the third column is "Position", 
##and subsequent columns are genotype data from either maize or teosinte individuals.
### Assumption that genotype data means nucleotides only

Cleaned_ZMP_Join_SNP_Geno = subset(Join_ZMP_SNP_Geno, select = -c(cdv_marker_id, alt_pos:count_gene))


## Chromosome Files for Teosinte
###Since these already have ? marks for missing data, we do not have to impute here

Chr1_Teosinte_Cleaned_Join_SNP_Geno= filter(Cleaned_ZMP_Join_SNP_Geno, Chromosome == 1)
Chr2_Teosinte_Cleaned_Join_SNP_Geno= filter(Cleaned_ZMP_Join_SNP_Geno, Chromosome == 2)
Chr3_Teosinte_Cleaned_Join_SNP_Geno= filter(Cleaned_ZMP_Join_SNP_Geno, Chromosome == 3)
Chr4_Teosinte_Cleaned_Join_SNP_Geno= filter(Cleaned_ZMP_Join_SNP_Geno, Chromosome == 4)
Chr5_Teosinte_Cleaned_Join_SNP_Geno= filter(Cleaned_ZMP_Join_SNP_Geno, Chromosome == 5)
Chr6_Teosinte_Cleaned_Join_SNP_Geno= filter(Cleaned_ZMP_Join_SNP_Geno, Chromosome == 6)
Chr7_Teosinte_Cleaned_Join_SNP_Geno= filter(Cleaned_ZMP_Join_SNP_Geno, Chromosome == 7)
Chr8_Teosinte_Cleaned_Join_SNP_Geno= filter(Cleaned_ZMP_Join_SNP_Geno, Chromosome == 8)
Chr9_Teosinte_Cleaned_Join_SNP_Geno= filter(Cleaned_ZMP_Join_SNP_Geno, Chromosome == 9)
Chr10_Teosinte_Cleaned_Join_SNP_Geno= filter(Cleaned_ZMP_Join_SNP_Geno, Chromosome == 10)

## Replace ? with -
P = as.data.frame(lapply(Cleaned_Join_SNP_Geno, gsub, pattern="?", replacement= "-", fixed= TRUE))
view(P)

R = as.data.frame(lapply(Cleaned_ZMP_Join_SNP_Geno, gsub, pattern="?", replacement= "-", fixed= TRUE))
view(R)

#Chromosome Files for Maize- ? swapped for -

Chr1_dash_maize_Cleaned_Join_SNP_Geno= filter(P, Chromosome == 1)
Chr2_dash_maize_Cleaned_Join_SNP_Geno= filter(P, Chromosome == 2)
Chr3_dash_maize_Cleaned_Join_SNP_Geno= filter(P, Chromosome == 3)
Chr4_dash_maize_Cleaned_Join_SNP_Geno= filter(P, Chromosome == 4)
Chr5_dash_maize_Cleaned_Join_SNP_Geno= filter(P, Chromosome == 5)
Chr6_dash_maize_Cleaned_Join_SNP_Geno= filter(P, Chromosome == 6)
Chr7_dash_maize_Cleaned_Join_SNP_Geno= filter(P, Chromosome == 7)
Chr8_dash_maize_Cleaned_Join_SNP_Geno= filter(P, Chromosome == 8)
Chr9_dash_maize_Cleaned_Join_SNP_Geno= filter(P, Chromosome == 9)
Chr10_dash_maize_Cleaned_Join_SNP_Geno= filter(P, Chromosome == 10)

#Chromosome files for Teosinte- ? swapped for -
Chr1_dash_Teosinte_Cleaned_Join_SNP_Geno= filter(R, Chromosome == 1)
Chr2_dash_Teosinte_Cleaned_Join_SNP_Geno= filter(R, Chromosome == 2)
Chr3_dash_Teosinte_Cleaned_Join_SNP_Geno= filter(R, Chromosome == 3)
Chr4_dash_Teosinte_Cleaned_Join_SNP_Geno= filter(R, Chromosome == 4)
Chr5_dash_Teosinte_Cleaned_Join_SNP_Geno= filter(R, Chromosome == 5)
Chr6_dash_Teosinte_Cleaned_Join_SNP_Geno= filter(R, Chromosome == 6)
Chr7_dash_Teosinte_Cleaned_Join_SNP_Geno= filter(R, Chromosome == 7)
Chr8_dash_Teosinte_Cleaned_Join_SNP_Geno= filter(R, Chromosome == 8)
Chr9_dash_Teosinte_Cleaned_Join_SNP_Geno= filter(R, Chromosome == 9)
Chr10_dash_Teosinte_Cleaned_Join_SNP_Geno= filter(R, Chromosome == 10)


#Visualization

##SNP's Per Chromosome Maize

ZMM_SNPS_Per_Chrom <- table(Cleaned_Join_SNP_Geno$Chromosome)
ZMM_SNPS_Per_Chrom<- melt(ZMM_SNPS_Per_Chrom)
ggplot(ZMM_SNPS_Per_Chrom, aes(x=as.factor(Var.1), y=value))+ geom_bar(stat='identity')
##SNP's Per Chromosome Teosinte

ZMP_SNPS_Per_Chrom <- table(Cleaned_ZMP_Join_SNP_Geno$Chromosome)
ZMP_SNPS_Per_Chrom<- melt(ZMP_SNPS_Per_Chrom)
ggplot(ZMP_SNPS_Per_Chrom, aes(x=as.factor(Var.1), y=value))+ geom_bar(stat='identity')
### Confirmation that they should have the same number

## Distribution of SNP'S on Chromosomes-Maize

ggplot(data = Cleaned_Join_SNP_Geno) + geom_point(mapping = aes(x=Position, y=Chromosome))

## Distribution of SNP'S on Chromosomes-Teosinte

ggplot(data = Cleaned_ZMP_Join_SNP_Geno) + geom_point(mapping = aes(x=Position, y=Chromosome))

## New Column _Maize
 ### Provides us data on the types of differences present in nucleotide data
 Cleaned_Join_SNP_Geno$ZygosityAT <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "A/T", "HZ", "Hm")
 Cleaned_Join_SNP_Geno$ZygosityAG <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "A/G", "HZ", "Hm")
 Cleaned_Join_SNP_Geno$ZygosityAC <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "A/C", "HZ", "Hm")
 Cleaned_Join_SNP_Geno$ZygosityTA <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "T/A", "HZ", "Hm")
 Cleaned_Join_SNP_Geno$ZygosityGA <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "G/A", "HZ", "Hm")
 Cleaned_Join_SNP_Geno$ZygosityCA <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "C/A", "HZ", "Hm")
 Cleaned_Join_SNP_Geno$ZygosityGC <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "G/C", "HZ", "Hm")
 Cleaned_Join_SNP_Geno$ZygosityGT <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "G/T", "HZ", "Hm")
 Cleaned_Join_SNP_Geno$ZygosityCG <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "C/G", "HZ", "Hm")
 Cleaned_Join_SNP_Geno$ZygosityTG <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "T/G", "HZ", "Hm")
 Cleaned_Join_SNP_Geno$ZygosityTC <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "T/C", "HZ", "Hm")
 Cleaned_Join_SNP_Geno$ZygosityCT <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "C/T", "HZ", "Hm")

 
 Cleaned_Join_SNP_Geno$Zygositymissing1 <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "?/*", Cleaned_Join_SNP_Geno$Chromosome, "Missing")
 Cleaned_Join_SNP_Geno$Zygositymissing2 <- ifelse(Cleaned_Join_SNP_Geno[c(4:1573)]== "*/?", Cleaned_Join_SNP_Geno$Chromosome, "Missing")
 
 
 Count_HZ= count(Cleaned_Join_SNP_Geno, vars="HZ")
 Count_HM= count(Cleaned_Join_SNP_Geno, vars="HM")
 
 #Update-
 ggplot(data = Cleaned_Join_SNP_Geno) + geom_point(mapping=aes(y="HZ", x="HM"))
 ggplot(data = Cleaned_Join_SNP_Geno) + geom_point(mapping=aes(y=SNP_ID, x="HM"))
 ggplot(data = Cleaned_Join_SNP_Geno) + geom_point(mapping=aes(y=SNP_ID, x="HZ"))
 
 ## New Column _Teosinte
 ### Provides us data on the types of differences present in nucleotide data

 Cleaned_ZMP_Join_SNP_Geno$ZygosityAT <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "A/T", "HZ", "Hm")
 Cleaned_ZMP_Join_SNP_Geno$ZygosityAG <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "A/G", "HZ", "Hm")
 Cleaned_ZMP_Join_SNP_Geno$ZygosityAC <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "A/C", "HZ", "Hm")
 Cleaned_ZMP_Join_SNP_Geno$ZygosityTA <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "T/A", "HZ", "Hm")
 Cleaned_ZMP_Join_SNP_Geno$ZygosityGA <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "G/A", "HZ", "Hm")
 Cleaned_ZMP_Join_SNP_Geno$ZygosityCA <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "C/A", "HZ", "Hm")
 Cleaned_ZMP_Join_SNP_Geno$ZygosityGC <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "G/C", "HZ", "Hm")
 Cleaned_ZMP_Join_SNP_Geno$ZygosityGT <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "G/T", "HZ", "Hm")
 Cleaned_ZMP_Join_SNP_Geno$ZygosityCG <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "C/G", "HZ", "Hm")
 Cleaned_ZMP_Join_SNP_Geno$ZygosityTG <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "T/G", "HZ", "Hm")
 Cleaned_ZMP_Join_SNP_Geno$ZygosityTC <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "T/C", "HZ", "Hm")
 Cleaned_ZMP_Join_SNP_Geno$ZygosityCT <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "C/T", "HZ", "Hm")
 
 Cleaned_ZMP_Join_SNP_Geno$Zygositymissing1 <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "?/*", Cleaned_Join_SNP_Geno$Chromosome, "Missing")
 Cleaned_ZMP_Join_SNP_Geno$Zygositymissing2 <- ifelse(Cleaned_ZMP_Join_SNP_Geno[c(4:978)]== "*/?", Cleaned_Join_SNP_Geno$Chromosome, "Missing")
 
 #Update
 ggplot(data = Cleaned_ZMP_Join_SNP_Geno) + geom_point(mapping=aes(y="HZ", x="HM"))
 ggplot(data = Cleaned_ZMP_Join_SNP_Geno) + geom_point(mapping=aes(y=SNP_ID, x="HM"))
 ggplot(data = Cleaned_ZMP_Join_SNP_Geno) + geom_point(mapping=aes(y=SNP_ID, x="HZ"))
 
 
  
  ## My own Visualization 
### Amplicons per chromosome-Maize

barplot(tapply(Join_SNP_Geno$count_amplicons, format(Join_SNP_Geno$Chromosome), FUN=sum))

### Amplicons per chromosome-Teosinte

barplot(tapply(Join_ZMP_SNP_Geno$count_amplicons, format(Join_SNP_Geno$Chromosome), FUN=sum))

write.table(Chr1_Maize_Cleaned_Join_SNP_Geno, "Chr1_Maize.txt", sep="\t")
write.table(Chr2_Maize_Cleaned_Join_SNP_Geno, "Chr2_Maize.txt", sep="\t")
write.table(Chr3_Maize_Cleaned_Join_SNP_Geno, "Chr3_Maize.txt", sep="\t")
write.table(Chr4_Maize_Cleaned_Join_SNP_Geno, "Chr4_Maize.txt", sep="\t")
write.table(Chr5_Maize_Cleaned_Join_SNP_Geno, "Chr5_Maize.txt", sep="\t")
write.table(Chr6_Maize_Cleaned_Join_SNP_Geno, "Chr6_Maize.txt", sep="\t")
write.table(Chr7_Maize_Cleaned_Join_SNP_Geno, "Chr7_Maize.txt", sep="\t")
write.table(Chr8_Maize_Cleaned_Join_SNP_Geno, "Chr8_Maize.txt", sep="\t")
write.table(Chr9_Maize_Cleaned_Join_SNP_Geno, "Chr9_Maize.txt", sep="\t")
write.table(Chr10_Maize_Cleaned_Join_SNP_Geno, "Chr10_Maize.txt", sep="\t")
write.table(Chr1_Teosinte_Cleaned_Join_SNP_Geno, "Chr1_Teosinte.txt", sep="\t")
write.table(Chr2_Teosinte_Cleaned_Join_SNP_Geno, "Chr2_Maize.txt", sep="\t")
write.table(Chr3_Teosinte_Cleaned_Join_SNP_Geno, "Chr3_Teosinte.txt", sep="\t")
write.table(Chr4_Teosinte_Cleaned_Join_SNP_Geno, "Chr4_Teosinte.txt", sep="\t")
write.table(Chr5_Teosinte_Cleaned_Join_SNP_Geno, "Chr5_Teosinte.txt", sep="\t")
write.table(Chr6_Teosinte_Cleaned_Join_SNP_Geno, "Chr6_Teosinte.txt", sep="\t")
write.table(Chr7_Teosinte_Cleaned_Join_SNP_Geno, "Chr7_Teosinte.txt", sep="\t")
write.table(Chr8_Teosinte_Cleaned_Join_SNP_Geno, "Chr8_Teosinte.txt", sep="\t")
write.table(Chr9_Teosinte_Cleaned_Join_SNP_Geno, "Chr9_Teosinte.txt", sep="\t")
write.table(Chr10_Teosinte_Cleaned_Join_SNP_Geno, "Chr10_Teosinte.txt", sep="\t")
write.table(Chr1_dash_maize_Cleaned_Join_SNP_Geno, "Chr1_dash_Maize.txt", sep="\t")
write.table(Chr2_dash_maize_Cleaned_Join_SNP_Geno, "Chr2_dash_Maize.txt", sep="\t")
write.table(Chr3_dash_maize_Cleaned_Join_SNP_Geno, "Chr3_dash_Maize.txt", sep="\t")
write.table(Chr4_dash_maize_Cleaned_Join_SNP_Geno, "Chr4_dash_Maize.txt", sep="\t")
write.table(Chr5_dash_maize_Cleaned_Join_SNP_Geno, "Chr5_dash_Maize.txt", sep="\t")
write.table(Chr6_dash_maize_Cleaned_Join_SNP_Geno, "Chr6_dash_Maize.txt", sep="\t")
write.table(Chr7_dash_maize_Cleaned_Join_SNP_Geno, "Chr7_dash_Maize.txt", sep="\t")
write.table(Chr8_dash_maize_Cleaned_Join_SNP_Geno, "Chr8_dash_Maize.txt", sep="\t")
write.table(Chr9_dash_maize_Cleaned_Join_SNP_Geno, "Chr9_dash_Maize.txt", sep="\t")
write.table(Chr10_dash_maize_Cleaned_Join_SNP_Geno, "Chr10_dash_Maize.txt", sep="\t")
write.table(Chr1_dash_Teosinte_Cleaned_Join_SNP_Geno, "Chr1_dash_Teosinte.txt", sep="\t")
write.table(Chr2_dash_Teosinte_Cleaned_Join_SNP_Geno, "Chr2_dash_Teosinte.txt", sep="\t")
write.table(Chr3_dash_Teosinte_Cleaned_Join_SNP_Geno, "Chr3_dash_Teosinte.txt", sep="\t")
write.table(Chr4_dash_Teosinte_Cleaned_Join_SNP_Geno, "Chr4_dash_Teosinte.txt", sep="\t")
write.table(Chr5_dash_Teosinte_Cleaned_Join_SNP_Geno, "Chr5_dash_Teosinte.txt", sep="\t")
write.table(Chr6_dash_Teosinte_Cleaned_Join_SNP_Geno, "Chr6_dash_Teosinte.txt", sep="\t")
write.table(Chr7_dash_Teosinte_Cleaned_Join_SNP_Geno, "Chr7_dash_Teosinte.txt", sep="\t")
write.table(Chr8_dash_Teosinte_Cleaned_Join_SNP_Geno, "Chr8_dash_Teosinte.txt", sep="\t")
write.table(Chr9_dash_Teosinte_Cleaned_Join_SNP_Geno, "Chr9_dash_Teosinte.txt", sep="\t")
write.table(Chr10_dash_Teosinte_Cleaned_Join_SNP_Geno, "Chr10_dash_Teosinte.txt", sep="\t")



  

