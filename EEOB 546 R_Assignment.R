library("tidyverse")
library("tidyr")

library("ggplot2")

setwd("C:/Users/Tanner/Desktop/Tanner/Iowa_State_Work/PhD_2020_2021/Bioinformatics_for_Biological_data/EEOB 546 R_Assignment")
Fang_genotypes= read_tsv("fang_et_al_genotypes.txt")
SNP_position= read_tsv("snp_position.txt")

view(Fang_genotypes)
view(SNP_position)
getwd()
