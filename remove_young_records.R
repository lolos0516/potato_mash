# lzhang
# remove the youngest generation records for XV
library(dplyr)

full <- read.table("blupf90/vce/renum_gbs/renf90.dat", header = F)
table(full$V5)
# 0 2015 2016 2017 
# 109  238   86  296
296/nrow(full) # 0.41

red <- full %>% 
  filter(V5!=2017)

write.table(red, "blupf90/vce/renum_gbs/renf90.dat.reduce", quote = F,
            sep = " ", row.names = F, col.names = F)
