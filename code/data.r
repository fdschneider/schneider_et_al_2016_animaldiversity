###############################################################################
# Code for the article 'Animal diversity and ecosystem functioning in dynamic food webs'
# 
# Copyright (C) 2016 Christian Guill & Florian D. Schneider
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################


##################################
### read in simulation results ###
##################################

# list of filenames
items <- paste("data/pdef_2.4/pdef_2.4_", c(paste0("0", as.character(1:9)), as.character(10:15)), sep ="")

# loop: open files and append to object webstats
webstats <- read.table(items[1], skip =1, header = T, sep = "")
for(i in 2:length(items)) webstats <- rbind(webstats, read.table(items[i], skip =1, header = T, sep = ""))

##################################
###### derived parameters ########
##################################

S_b_ini <- 30 # initial number of basal species
webstats$S_ini <- webstats$S_c+S_b_ini        # initial number of species
webstats$S_c_fin <- round(webstats$S_c*webstats$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats$S_c_fin[webstats$S_c_fin == "NaN"] <- 0                     # sets NaN = 0
webstats$S_b_fin <- round(S_b_ini*webstats$P_b, digits = 0)		# final number of basal species
webstats$S_fin <- webstats$S_c_fin+webstats$S_b_fin				# final number of species
webstats$FP <- webstats$FD_fin/webstats$FD_ini				# functional persistence
webstats$iniC <- webstats$iniL/webstats$S_ini^2     # initial Connectance
webstats$finC <- webstats$finL/webstats$S_fin^2    # final Connectance

webstats$meanB_b <- webstats$meanB_b*webstats$S_b_fin
webstats$meanB_c <- webstats$meanB_c*webstats$S_c_fin

webstats$Mav <- webstats$meanB/(webstats$findens_c+webstats$findens_b) #  average bodymass
webstats$Mav_c <- webstats$geom_mean_m_c #webstats$meanB_c/webstats$findens_c   # average bodymass of predators
webstats$Mav_b <- webstats$geom_mean_m_p #webstats$meanB_b/webstats$findens_b  # average bodymass of basalspecies

webstats$meanTLav_ini[which(webstats$meanTLav_ini == "NaN")] <- 0   # replaces NaN with 0
webstats$meanTLav_fin[which(webstats$meanTLav_fin == "NaN")] <- 0   # replaces NaN with 0
webstats$meanTLmin_ini[which(webstats$meanTLmin_ini == "NaN")] <- 0   # replaces NaN with 0
webstats$meanTLmin_fin[which(webstats$meanTLmin_fin == "NaN")] <- 0   # replaces NaN with 0 
webstats$meanTLeff[which(webstats$meanTLeff == "NaN")] <- 0    # replaces NaN with 0 


##################################
##### drawing valid datasets  ####
##################################

# from all generated webs, draw valid ones where at least one predator survived.  
#set.seed(11) 
#webstats <- webstats[sort(sample(webstats$WebID[webstats$fintop_b == 0  & webstats$meanB_c >= 0 & webstats$meanB_b >= 0], 7500, replace = FALSE)),]   # selects 40.000 random food webs
length(which(webstats$fintop_b <= 0))
webstats <- subset(webstats, fintop_b == 0  & meanB_b > 0)
#webstats <- subset(webstats, meanB_b > 0 & meanB_c > 0)

write.table(webstats, "data/webstats.txt", sep = ",")



###########################################
###  Read data for sensitivity analysis ###
###########################################

items <- c(
  "data/pdef_2.4/pdef_2.4_H50_01",
  "data/pdef_2.4/pdef_2.4_H50_02",
  "data/pdef_2.4/pdef_2.4_H50_03"
)

# loop: open files and append to object webstats
webstats_H50 <- read.table(items[1], skip =1, header = T, sep = "")
for(i in 2:length(items)) webstats_H50 <- rbind(webstats_H50, read.table(items[i], skip =1, header = T, sep = ""))
webstats_H50$S_c_fin <- round(webstats_H50$S_c*webstats_H50$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats_H50 <- subset(webstats_H50, fintop_b == 0  & meanB_b > 0)

webstats_H50$S_c_fin <- round(webstats_H50$S_c*webstats_H50$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats_H50$S_c_fin[webstats_H50$S_c_fin == "NaN"] <- 0                     # sets NaN = 0
webstats_H50$S_b_fin <- round(S_b_ini*webstats_H50$P_b, digits = 0)		# final number of basal species

webstats_H50$meanB_b <- webstats_H50$meanB_b*webstats_H50$S_b_fin
webstats_H50$meanB_c <- webstats_H50$meanB_c*webstats_H50$S_c_fin

webstats_H50$H50 <- 50

write.table(webstats_H50, "data/webstats_H50.txt", sep = ",")


###########################################

items <- c(
  "data/pdef_2.4/pdef_2.4_Sb10_01",
  "data/pdef_2.4/pdef_2.4_Sb10_03",
  "data/pdef_2.4/pdef_2.4_Sb10_02",
  "data/pdef_2.4/pdef_2.4_Sb50_01",
  "data/pdef_2.4/pdef_2.4_Sb50_02",
  "data/pdef_2.4/pdef_2.4_Sb50_02"
)

# loop: open files and append to object webstats
webstats_Sb10 <- read.table(items[1], skip =1, header = T, sep = "")
for(i in 2:3) webstats_Sb10 <- rbind(webstats_Sb10, read.table(items[i], skip =1, header = T, sep = ""))

webstats_Sb10$S_c_fin <- round(webstats_Sb10$S_c*webstats_Sb10$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats_Sb10$S_c_fin[webstats_Sb10$S_c_fin == "NaN"] <- 0                     # sets NaN = 0
webstats_Sb10 <- subset(webstats_Sb10, fintop_b == 0 & S_c_fin != 0 & meanB_b > 0)


webstats_Sb10$S_c_fin <- round(webstats_Sb10$S_c*webstats_Sb10$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats_Sb10$S_c_fin[webstats_Sb10$S_c_fin == "NaN"] <- 0                     # sets NaN = 0
webstats_Sb10$S_b_fin <- round(S_b_ini*webstats_Sb10$P_b, digits = 0)		# final number of basal species

webstats_Sb10$meanB_b <- webstats_Sb10$meanB_b*webstats_Sb10$S_b_fin
webstats_Sb10$meanB_c <- webstats_Sb10$meanB_c*webstats_Sb10$S_c_fin

write.table(webstats_Sb10, "data/webstats_Sb10.txt", sep = ",")

webstats_Sb50 <- read.table(items[4], skip =1, header = T, sep = "")
for(i in 5:6) webstats_Sb50 <- rbind(webstats_Sb50, read.table(items[i], skip =1, header = T, sep = ""))
webstats_Sb50$S_c_fin <- round(webstats_Sb50$S_c*webstats_Sb50$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats_Sb50$S_c_fin[webstats_Sb50$S_c_fin == "NaN"] <- 0
webstats_Sb50 <- subset(webstats_Sb50, fintop_b == 0 & S_c_fin != 0 & meanB_b > 0)

webstats_Sb50$S_c_fin <- round(webstats_Sb50$S_c*webstats_Sb50$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats_Sb50$S_c_fin[webstats_Sb50$S_c_fin == "NaN"] <- 0                     # sets NaN = 0
webstats_Sb50$S_b_fin <- round(S_b_ini*webstats_Sb50$P_b, digits = 0)		# final number of basal species

webstats_Sb50$meanB_b <- webstats_Sb50$meanB_b*webstats_Sb50$S_b_fin
webstats_Sb50$meanB_c <- webstats_Sb50$meanB_c*webstats_Sb50$S_c_fin

write.table(webstats_Sb50, "data/webstats_Sb50.txt", sep = ",")

###########################################

items <- c(
  "data/pdef_2.4/pdef_2.4_D01_01",
  "data/pdef_2.4/pdef_2.4_D01_02",
  "data/pdef_2.4/pdef_2.4_D01_03",
  "data/pdef_2.4/pdef_2.4_D05_01",
  "data/pdef_2.4/pdef_2.4_D05_02",
  "data/pdef_2.4/pdef_2.4_D05_03"
)

# loop: open files and append to object webstats
webstats_D01 <- read.table(items[1], skip =1, header = T, sep = "")
for(i in 2:3) webstats_D01 <- rbind(webstats_D01, read.table(items[i], skip =1, header = T, sep = ""))
webstats_D01$S_c_fin <- round(webstats_D01$S_c*webstats_D01$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats_D01$S_c_fin[webstats_D01$S_c_fin == "NaN"] <- 0                     # sets NaN = 0
webstats_D01 <- subset(webstats_D01, fintop_b == 0 & S_c_fin != 0 & meanB_b > 0)

webstats_D01$S_c_fin <- round(webstats_D01$S_c*webstats_D01$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats_D01$S_c_fin[webstats_D01$S_c_fin == "NaN"] <- 0                     # sets NaN = 0
webstats_D01$S_b_fin <- round(S_b_ini*webstats_D01$P_b, digits = 0)		# final number of basal species

webstats_D01$meanB_b <- webstats_D01$meanB_b*webstats_D01$S_b_fin
webstats_D01$meanB_c <- webstats_D01$meanB_c*webstats_D01$S_c_fin

write.table(webstats_D01, "data/webstats_D01.txt", sep = ",")

webstats_D05 <- read.table(items[4], skip =1, header = T, sep = "")
for(i in 5:6) webstats_D05 <- rbind(webstats_D05, read.table(items[i], skip =1, header = T, sep = ""))

webstats_D05$S_c_fin <- round(webstats_D05$S_c*webstats_D05$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats_D05$S_c_fin[webstats_D05$S_c_fin == "NaN"] <- 0
webstats_D05 <- subset(webstats_D05, fintop_b == 0 & S_c_fin != 0 & meanB_b > 0)

webstats_D05$S_c_fin <- round(webstats_D05$S_c*webstats_D05$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats_D05$S_c_fin[webstats_D05$S_c_fin == "NaN"] <- 0                     # sets NaN = 0
webstats_D05$S_b_fin <- round(S_b_ini*webstats_D05$P_b, digits = 0)		# final number of basal species

webstats_D05$meanB_b <- webstats_D05$meanB_b*webstats_D05$S_b_fin
webstats_D05$meanB_c <- webstats_D05$meanB_c*webstats_D05$S_c_fin

write.table(webstats_D05, "data/webstats_D05.txt", sep = ",")

