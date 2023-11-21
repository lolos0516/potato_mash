library(dplyr)

gbs_solution <- read.table("blupf90/vce/vce_gbs/solutions", header = F, skip = 1,
                           col.names = c("trait", "effect", "level",  "solution", " s.e."))
blend_solution <- read.table("blupf90/vce/vce_blend/solutions", header = F, skip = 1,
                           col.names = c("trait", "effect", "level",  "solution", " s.e."))
hap_solution <- read.table("blupf90/vce/vce_hap/solutions", header = F, skip = 1,
                           col.names = c("trait", "effect", "level",  "solution", " s.e."))
dat <- read.table("blupf90/vce/renum_gbs/renf90.dat", header = F, 
                  col.names = c("trait", "year", "id", "id_orig", "year_orig"))
young <- dat$id[dat$year_orig==2017]

gbs_sub <- filter(gbs_solution, level %in% young)
hap_sub <- filter(hap_solution, level %in% young)
blend_sub <- filter(blend_solution, level %in% young)

cor(cbind(gbs_sub$solution, hap_sub$solution, blend_sub$solution))
# [,1]     [,2]      [,3]
# [1,] 1.0000000 0.285907 0.2806476
# [2,] 0.2859070 1.000000 0.7016790
# [3,] 0.2806476 0.701679 1.0000000

cor(cbind(gbs_solution$solution, hap_solution$solution, blend_solution$solution))
# [,1]      [,2]      [,3]
# [1,] 1.0000000 0.2445246 0.8162987
# [2,] 0.2445246 1.0000000 0.5701508
# [3,] 0.8162987 0.5701508 1.0000000
