# lzhang


library(readxl)
library(tidyr)
library(dplyr)
library(data.table)

####### read
df_id <- read_xlsx("FRY_OTF_LTS_phenotypes.xlsx", sheet = 1)
dim(df_id) # 810 4
head(df_id)
# `Master Genopheno line name` `sample name in SMAP outputs`                                                       OTF                LTS  
# <chr>                        <chr>                                                                               <chr>              <chr>
# 1 ANTARCTICA                   MOLP2_pl_7_E11_GenSPI_truetest_ANTARCTICA_1.fq.gz_.extendedFrags.fastq_Q30_sort.bam 54.748319090000003 NA   
# 2 BIKINI                       NA                                                                                  NA                 NA   
# 3 CARA                         NA                                                                                  NA                 NA 

df_id2 <- read_xlsx("FRY_OTF_LTS_phenotypes.xlsx", sheet = 2)
dim(df_id2) # 765 5
head(df_id2)
# `sample name in SMAP outputs`                                            OTF                LTS              `GenSPI file name` Master Genopheno lin…¹
# <chr>                                                                    <chr>              <chr>            <chr>              <chr>                 
# MOLP2_pl_1_A10_GenSPI_2015_99_1.fq.gz_.extendedFrags.fastq_Q30_sort.bam  25.282499999999999 25.93            T15_133_GenSPI.fq… T15_133               
# MOLP2_pl_1_A11_GenSPI_2015_141_1.fq.gz_.extendedFrags.fastq_Q30_sort.bam 48.8125            39.4925          T15_175_GenSPI.fq… T15_175               
# MOLP2_pl_1_A12_GenSPI_2015_55_1.fq.gz_.extendedFrags.fastq_Q30_sort.bam  53.39              39.82            T15_088_GenSPI.fq… T15_088 

df_hap <- fread("data/FRY_SMAP_outputs/haplotypes_c20_f5_m0_discrete_calls_filtered.tsv", header = T)
df_hap[1:3,1:3]
# Locus                  Haplotypes       MOLP2_pl_1_A10_GenSPI_2015_99_1.fq.gz_.extendedFrags.fastq_Q30_sort.bam
# chr01:1489825-1489947/+          0                                                                    <NA>
# chr01:1489825-1489947/+         10                                                                    <NA>
# chr01:1489825-1489947/+        100                                                                    <NA>

df_gbs <- fread("data/MASTER_GENOPHENO/MASTER_GENOPHENO", header=T)
dim(df_gbs) # 810 46410
df_gbs[1:4,1:7]
# popOL    pop       line  crispBLUP_OTF LTS chr01_262875 chr01_262876
# truetest test ANTARCTICA      54.74832  NA            0            1
# truetest test     BIKINI      40.21322  NA           -1            0
# no       test       CARA      39.02520  NA            0           NA

sum(!is.na(match(tail(names(df_hap), -2), df_id$`sample name in SMAP outputs`)))
# 327
sum(!is.na(match(tail(names(df_hap), -2), df_id2$`sample name in SMAP outputs`)))
# 387

######## subset haplotype lines existing in id file
idx <- which(names(df_hap) %in% df_id2$`sample name in SMAP outputs`)
df_hap_sub <- df_hap[,..idx]
df_hap_sub[1:3,1:3]

idx <- match(names(df_hap_sub), df_id2$`sample name in SMAP outputs`)
names(df_hap_sub) <- df_id2$`Master Genopheno line name`[idx]
df_hap_sub[1:3,1:3]
df_hap_sub_t <- t(df_hap_sub) %>% data.frame()
df_hap_sub_t[1:3,1:3]
df_hap_sub_t <- data.frame(line = rownames(df_hap_sub_t), df_hap_sub_t)
df_hap_sub <- df_hap_sub_t
for(i in 2:ncol(df_hap_sub_t)) {
  df_hap_sub[[i]] <- as.integer(df_hap_sub[[i]])
}
df_hap_sub_t[89:92, 49:56]
df_hap_sub[89:92, 49:56]

sapply(df_hap_sub[, -c(1)], max, na.rm = T) %>% table()
# 0   1   2   3   4 
# 82 299 470 536 613 

####### replace haplotypes 1-3 to 1 and 4 to 2
for (i in 1:ncol(df_hap_sub)) {
  df_hap_sub[[i]][which(df_hap_sub[[i]] %in% 1:3)] <- 1
  df_hap_sub[[i]][which(df_hap_sub[[i]] == 4)] <- 2
  df_hap_sub[[i]][which(is.na(df_hap_sub[[i]]))] <- 5
}
df_hap_sub[89:92, 49:56]
sapply(df_hap_sub[, -c(1)], max) %>% table()

###### order gbs to have lines in haplotype file at the very end
idx <- match(df_hap_sub$line, df_gbs$line)
idx <- na.omit(idx)
df_gbs_common <- df_gbs[idx,] %>% 
  select(-starts_with("pop"))
df_gbs_notco <- df_gbs %>% 
  select(-starts_with("pop")) %>% 
  filter(!line %in% df_gbs_common$line)

# sanity check
sum(!is.na(match(df_gbs_common$line, df_gbs_notco$line)))==0
nrow(df_gbs_common) + nrow(df_gbs_notco) == nrow(df_gbs)

df_gbs <- rbind(df_gbs_notco, df_gbs_common)

# change -1 0 1 to 0 1 2
df <- lapply(4:ncol(df_gbs), function(i) {
  col <- df_gbs[[i]]
  col[which(col==1)] <- 2
  col[which(col==0)] <- 1
  col[which(col==-1)] <- 0
  col[is.na(col)] <- 5
  return(col)
}) %>% do.call(what = cbind)

df <- data.frame(df)
names(df) <- names(df_gbs)[4:ncol(df_gbs)]
df_gbs <- cbind(df_gbs[,c(1:3)], df)

fwrite(df_gbs[,-c(2:3)], "output/gbs", quote = F, sep = " ", row.names = F, col.names = F)

##### remove lines in haplotype that are not in gbs
df_hap_sub <- df_hap_sub %>% 
  filter(line %in% df_gbs_common$line)

fwrite(df_hap_sub, "output/haplotype", quote = F, sep = " ", row.names = F, col.names = F)
