# lzhang
# manually create Hi = Ai + Gi - A22i
library(dplyr)
library(reshape2)

####### function
convert_to_matrix <- function(df, xrefid_df) {
  rank <- max(df$i)
  
  out <- matrix(rep(NaN, rank*rank), rank, rank)
  
  for(row in 1:nrow(df)) {
    out[df$i[row], df$j[row]] <- df$gii[row]
  }
 rownames(out) <- xrefid_df$id
 colnames(out) <- xrefid_df$id
  return(out)
}

calc_hi <- function(ai_mat, gima22i) {
  # make 
  # |0, 0       |
  # |0, gi-a22i |
  place_holder <- matrix(rep(NaN, nrow(ai_mat)*nrow(ai_mat)), nrow(ai_mat), nrow(ai_mat))
  place_holder[upper.tri(place_holder)] <- 0
  diag(place_holder) <- 0
  colnames(place_holder) <- colnames(ai_mat)
  rownames(place_holder) <- rownames(ai_mat)
  
  for(i in 1:nrow(gima22i)) {
    for(j in i:ncol(gima22i)) {
      row_id <- rownames(gima22i)[i]
      col_id <- colnames(gima22i)[j]
      place_holder[row_id,col_id] <- gima22i[row_id,col_id] 
    }
  }
  
  # hi = ai + above
  hi_mat <- ai_mat + place_holder
  colnames(hi_mat) <- as.character(1:ncol(hi_mat))
  rownames(hi_mat) <- as.character(1:nrow(hi_mat))
  
  return(hi_mat)
}

convert_to_table <- function(hi_mat) {
  
}

##### global par
tau <- 1
omega <- 0

###### read
# fake a22i0
a22i <- read.table("blupf90/vce/pregs_gbs_a22/Gi_ascii", header=F, 
                    col.names=c("i", "j", "gii"))
xrefid_a22 <- read.table("output/gbs_a22_clean_XrefID", header = F,
                     col.names = c("id", "id_orig")) %>% 
  select(id)
# fake gi
gi <- read.table("blupf90/vce/pregs_hap/Gi_ascii", header =F, 
                    col.names=c("i", "j", "gii"))
xrefid_gi <- read.table("output/haplotype_clean_XrefID", header = F,
                        col.names = c("id", "id_orig")) %>% 
  select(id)

# fake ai
ai <- read.table("blupf90/vce/pregs_gbs/Gi_ascii", header =F, 
                 col.names=c("i", "j", "gii"))
xrefid_ai <- read.table("output/gbs_clean_XrefID", header = F,
                        col.names = c("id", "id_orig")) %>% 
  select(id)

# fake a22
a22 <- read.table("blupf90/vce/pregs_gbs_a22/G", header=F, 
                   col.names=c("i", "j", "gii"))
# fake g
g <- read.table("blupf90/vce/pregs_hap/G", header =F, 
                 col.names=c("i", "j", "gii"))

# sanity check
max(a22i$i);max(gi$i);max(ai$i) # 702, 702, 806
match(xrefid_gi$id, xrefid_a22$id) %>% is.na() %>% sum() # 0
match(xrefid_ai$id, xrefid_gi$id) %>% is.na() %>% sum() # 104

# matrix
a22i_mat <- convert_to_matrix(a22i, xrefid_a22)
gi_mat <- convert_to_matrix(gi, xrefid_gi)
ai_mat <- convert_to_matrix(ai, xrefid_ai)
a22_mat <- convert_to_matrix(a22, xrefid_a22)
g_mat <- convert_to_matrix(g, xrefid_gi)

# sumstat. Their means are very similar. Don't need to tune.
mean(diag(g_mat));sd(diag(g_mat)) # 0.5656 0.0047
mean(diag(a22_mat));sd(diag(a22_mat)) # 0.5094 0.1040
mean(g_mat[upper.tri(g_mat)]);sd(g_mat[upper.tri(g_mat)]) # -0.000757 0.0385
mean(a22_mat[upper.tri(a22_mat)]);sd(a22_mat[upper.tri(a22_mat)]) # -0.0005795 0.0284

# Scale G based on A22.
# The variable x can be:
# 0: no scaling
# 1: mean(diag(G))=1, mean(offdiag(G))=0. This implies that the estimated variance components and mean 
# refer to the genotyped population Legarra, 2016
# 2: mean(diag(G))=mean(diag(A22)), mean(offdiag(G))=mean(offdiag(A22)) Christensen et al., 2012 
# This is the default
# 3: meann(G)=mean(A22) Vitezica et al. (2011)
# 4: rescale G using Fst adjustment. As in Powell et al. (2010) and Vitezica et al. (2011)
# 9: arbitrary parameters: specify two additional numbers a and b in a +  bG as OPTION tunedG 9 a b.
tune_g <- function(gi_mat, a22i_mat) {
  gi_mat-a22i_mat
}
gima22i <- tune_g(gi_mat, a22i_mat)

hi_mat <- calc_hi(ai_mat, gima22i)

# sanity check
hi_mat[105:107,105:107]
ai_mat[105:107,105:107]
hi_mat[1:3,1:3]*6
ai_mat[1:3,1:3]
summary(diag(hi_mat))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -87.679   6.114   6.774   7.021   7.359 100.000 
hist(diag(hi_mat), breaks=100)
hist(diag(hi_mat), breaks=100, ylim = c(0,10))
filter(hi, value<0 & Var1==Var2) # 1 row 118, -87
filter(ai, gii<0 &i==j)
filter(a22i, gii<0 &i==j)
filter(gi, gii<0 &i==j)

filter(ai, i==118 & j==118) # 6.5698
id <- xrefid_ai$id[118] # 40
gima22i[as.character(id), as.character(id)] # -94
a22i_mat[as.character(id), as.character(id)] # 100
gi_mat[as.character(id), as.character(id)] # 5.75
a22_mat[as.character(id), as.character(id)] # 0.01
g_mat[as.character(id), as.character(id)] # 0.584
summary(diag(g_mat-a22_mat))
summary(diag(gi_mat-a22i_mat))

summary(hi_mat[upper.tri(hi_mat)])
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -38.24235  -0.01053   0.11456   0.11231   0.24277  11.24201 

# output
hi <- melt(hi_mat, na.rm = T, as.is = T) 
head(hi)
# optional change negative diagonal to positive
hi$value[which(hi$value<0 & hi$Var1==hi$Var2)] <- 7

write.table(hi, "blupf90/vce/Hinv_manual", quote = F, sep = " ", row.names = F, col.names = F)
