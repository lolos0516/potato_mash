# lzhang
# manually create Hi = Ai + Gi - A22i
library(dplyr)
library(reshape2)
library(ggcorrplot)
library(stringr)

####### function
convert_to_matrix <- function(df) {
  rank <- length(unique(df$i))
  
  out <- matrix(rep(NaN, rank*rank), rank, rank)
  colnames(out) <- df$i[df$j==df$j[1]]
  rownames(out) <- colnames(out)
  
  for(row in 1:nrow(df)) {
    out[df$i[row], df$j[row]] <- df$gij[row]
  }
  return(out)
}

calc_hi <- function(ai_mat, gima22i) {
  # make 
  # |0, 0       |
  # |0, gi-a22i |
  # due to the id order in Ai and others are different, need to index them
  hi_mat <- ai_mat
  
  # hi = ai + gi - a22i
  for(i in rownames(gima22i)) {
    for(j in colnames(gima22i)) {
      hi_mat[i,j] <- hi_mat[i,j]+gima22i[i,j]
    }
  }
  
  return(hi_mat)
}

convert_mat_to_long <- function(hi_mat, xrefid) {
  hi <- melt(hi_mat, na.rm = T, as.is = T) 
  
  # change id to renum id
  hi <- left_join(hi, xrefid, by = c("Var1" = "id_orig")) %>% 
    left_join(xrefid, by = c("Var2" = "id_orig"))
# print(head(hi))
  hi <- hi %>% 
    filter(id.x <= id.y) %>% # only keep upper tri
    select(id.x, id.y, value)
  names(hi) <- c("Var1", "Var2", "value")
  return(hi)
}

##### global par
tau <- 1
omega <- 0

###### read
# fake a22i
a22i <- read.table("blupf90/vce/pregs_gbs_a22/Gi_Orig.txt", header=F, col.names=c("i", "j", "gij"))
# fake gi
gi <- read.table("blupf90/vce/pregs_hap/Gi_Orig.txt", header =F, col.names=c("i", "j", "gij"))
# fake ai
ai <- read.table("blupf90/vce/pregs_gbs/Gi_Orig.txt", header =F, col.names=c("i", "j", "gij"))
# fake a22
a22 <- read.table("blupf90/vce/pregs_gbs_a22/G_Orig.txt", header=F, col.names=c("i", "j", "gij"))
# fake g
g <- read.table("blupf90/vce/pregs_hap/G_Orig.txt", header =F, col.names=c("i", "j", "gij"))
# xrefid
ped <- read.table("blupf90/vce/renum_gbs/renadd02.ped", header = F)
names(ped)[1] <- "id"; names(id)[10] <- "id_orig"
# for the convenience of idx 
gi_pre <- read.table("blupf90/vce/pregs_gbs/Gi", header = F, col.names=c("i", "j", "gij"))

# in BLUPF90, Ainv ID order is different from G/Gi/A22i/A22 orders. Therefore, the new Hinv need to 
# follow the same order.

# sanity check
max(a22i$i);max(gi$i);max(ai$i) # 702, 702, 806
all(a22i$i==gi$i);all(gi$i==g$i)

# matrix
a22i_mat <- convert_to_matrix(a22i)
gi_mat <- convert_to_matrix(gi)
ai_mat <- convert_to_matrix(ai)
a22_mat <- convert_to_matrix(a22)
g_mat <- convert_to_matrix(g)

# sumstat. Their means are very similar. Don't need to tune.
mean(diag(g_mat));sd(diag(g_mat)) # 0.5656 0.0047
mean(diag(a22_mat));sd(diag(a22_mat)) # 0.5094 0.1040
mean(g_mat[upper.tri(g_mat)]);sd(g_mat[upper.tri(g_mat)]) # -0.000757 0.0385
mean(a22_mat[upper.tri(a22_mat)]);sd(a22_mat[upper.tri(a22_mat)]) # -0.0005795 0.0284

cor(cbind(diag(g_mat), diag(a22_mat))) # 0.331
cor(cbind(g_mat[upper.tri(g_mat)], a22_mat[upper.tri(a22_mat)])) # 0.771

mean(diag(g_mat)-diag(a22_mat));sd(diag(g_mat)-diag(a22_mat)) # 0.056 0.099
hist(diag(g_mat)-diag(a22_mat), breaks = 100)

mean(g_mat[upper.tri(g_mat)]-a22_mat[upper.tri(a22_mat)]) #  -0.0001778824
sd(g_mat[upper.tri(g_mat)]-a22_mat[upper.tri(a22_mat)]) # 0.02455699
hist(g_mat[upper.tri(g_mat)]-a22_mat[upper.tri(a22_mat)], breaks = 100) 

# too slow
# ggplot(g %>% mutate(g_a22=g$gij-a22$gij), aes(i, j, fill = g_a22)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "#075AFF",
#                        mid = "#FFFFCC",
#                        high = "#FF0000")

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
tune_g <- function(gi_mat, a22i_mat, tau=1, omega=1) {
  tau*gi_mat-omega*a22i_mat
}
gima22i <- tune_g(gi_mat, a22i_mat)

hi_mat <- calc_hi(ai_mat, gima22i)
hi <- convert_mat_to_long(hi_mat, ped[,c(1,10)]) 

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
# optional change negative diagonal to positive
hi$value[which(hi$value<0 & hi$Var1==hi$Var2)] <- mean(diag(hi_mat))

# sort to gi sequence
hi <- hi %>% 
  mutate(ij = paste0(Var1, "_", Var2))
gi_pre <- gi_pre %>% 
  mutate(ij = paste0(i, "_", j))
idx <- match(gi_pre$ij, hi$ij)
hi <- hi[idx,]

dfstr = data.frame(Var1=hi$Var1, #str_pad(hi$Var1, 10, "left"),
                   Var2=hi$Var1, #str_pad(hi$Var2, 10, "left"),
                   value=sprintf(paste0('%.', 14-1*(hi$value<0),'f'), hi$value))
write.table(dfstr,
            "blupf90/vce/Hinv_sorted", quote = F, sep = " ", row.names = F, col.names = F)

dfstr = data.frame(Var1=str_pad(hi$Var1, 10, "left"),
                   Var2=str_pad(hi$Var2, 10, "left"),
                   value=sprintf(paste0('%.', 21-1*(hi$value<0),'f'), hi$value))
write.table(dfstr, "blupf90/vce/Hinv_manual", quote = F, sep = " ", row.names = F, col.names = F)

