# lzhang
# after renuming gbs manually create an Xref ID for haplotype

library(dplyr)

gbs_xrefid <- read.table("output/gbs_XrefID", header = F, col.names = c("id", "id_orig"))
haplotype <- read.table("output/haplotype", header = F)
a22 <- read.table("output/gbs_a22", header = F)
#haplotype[1:3,1:3]
names(haplotype)[1] <- "id_orig"
names(a22)[1] <- "id_orig"

xrefid <- haplotype %>% 
  select(id_orig) %>% 
  left_join(gbs_xrefid) %>% 
  select(id, id_orig)

xrefid %>% filter(is.na(id))
# id  id_orig
# NA  T17_667
# NA  T17_599
# NA T5156_11
write.table(xrefid, "output/haplotype_XrefID", quote = F, col.names = F, row.names = F, sep = " ")

xrefid <- a22 %>% 
  select(id_orig) %>% 
  left_join(gbs_xrefid) %>% 
  select(id, id_orig)

head(xrefid)
write.table(xrefid, "output/gbs_a22_XrefID", quote = F, col.names = F, row.names = F, sep = " ")
