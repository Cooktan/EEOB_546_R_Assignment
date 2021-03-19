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

#Chromosome files for Teosinte
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

## New Column
New_Teosinte_Col= if (Cleaned_Join_SNP_Geno[1:987]== "A/A" | "C/C" | "G/G" | "T/T")
New_Teosinte_Col=Cleaned_Join_SNP_Geno$Zygosity
Cleaned_Join_SNP_Geno$Zygosity= c(1:989){Cleaned_Join_SNP_Geno[c(3:989)== "A/A" | "C/C" | "G/G" | "T/T"]== "HM"

  New_Column_Maize = as.data.frame(lapply(Cleaned_Join_SNP_Geno, gsub, pattern="A/A"| "G/G" | "T/T"| "C/C", replacement= "HM", fixed= TRUE))
  view(New_Column_Maize)
  
  Cleaned_Join_SNP_Geno$Zygosity <-ifelse(Cleaned_Join_SNP_Geno[c(3:983)]=="A/A" | "C/C" | "G/G" | "T/T", "1","0") 
  
  ## My own Visualization 
### Amplicons per chromosome-Maize

barplot(tapply(Join_SNP_Geno$count_amplicons, format(Join_SNP_Geno$Chromosome), FUN=sum))

### Amplicons per chromosome-Teosinte

barplot(tapply(Join_ZMP_SNP_Geno$count_amplicons, format(Join_SNP_Geno$Chromosome), FUN=sum))



