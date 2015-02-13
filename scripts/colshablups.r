#####workflow#####
#reshape data so each genotype has associated phenotypes
#start line 94
#run variable slope and intercept model for each trait
#get blups for each genotype in each condition
#map traits using increasingly complex QTL mapping techniques
#start with biomass traits 


#merge data frames
setwd("~/git.repos/NAM/input/")
envelope <- read.table("envelopenumtoplantnumrep3_cleaned.csv", 
	                    header=TRUE, sep = ",", na.strings = "-", stringsAsFactors = FALSE)
head(envelope)
str(envelope)
dim(envelope)

positions <- read.table("plantID_positions_set3_cleaned.csv", 
	                    header=TRUE, sep = ",", na.strings = "-", stringsAsFactors = FALSE)
head(positions)
str(positions)
colnames(positions)[3] <- paste("flatbarcode")
colnames(positions)[4] <- paste("flat_position")
dim(positions)

weights <- read.table("weights-rep3final_cleaned.csv", 
	                    header=TRUE, sep = ",", na.strings = "-", stringsAsFactors = FALSE)
head(weights)
str(weights)
dim(weights)

?merge
env_weight <- merge(envelope, weights, all.y = TRUE)
dim(env_weight)
head(env_weight)


env_weight_pos <- merge(positions, env_weight, all.y = TRUE)
head(env_weight_pos)
dim(env_weight_pos)
str(env_weight_pos)

env_weight_pos$trt <- sub("(3)(\\w)(\\d+)","\\2", env_weight_pos$flatbarcode)
head(env_weight_pos)
tail(env_weight_pos)

#remove last two columns of NA values

dim(env_weight_pos)
env_weight_pos <- env_weight_pos[-c(1744:1745),]

#move trt to front of dataframe
colnames(env_weight_pos)
env_weight_pos <- env_weight_pos[,c(13,1:12)]
head(env_weight_pos)


flowers <- read.table("colsha_traits.csv", 
	                    header=TRUE, sep = ",", na.strings = "-", stringsAsFactors = FALSE)
head(flowers)
head(env_weight_pos)
dim(flowers)

env_wght_pos_flr <- merge(env_weight_pos, flowers, all.x = TRUE)
dim(env_wght_pos_flr)
head(env_wght_pos_flr)
colnames(env_wght_pos_flr)
#remove duplicated metadata
env_wght_pos_flr <- env_wght_pos_flr[,-c(14:17)]
dim(env_wght_pos_flr)
head(env_wght_pos_flr)
colnames(env_wght_pos_flr)


setwd("~/git.repos/NAM/output/")
write.table(env_wght_pos_flr, file="col_sha_NAM_GWAS_phenotypes.csv", sep = ",")

#slice data into GWAS, Parents, Col/Sha
library(plyr)
?arrange
env_wght_pos_flr_srt <- env_wght_pos_flr[order(env_wght_pos_flr$genotype),]
head(env_wght_pos_flr_srt)
env_wght_pos_flr_srt$genotype

#colsha phenotypes only
colsha_ril_phenotypes <- env_wght_pos_flr_srt[1:1430,]
head(colsha_ril_phenotypes)
tail(colsha_ril_phenotypes)
#the last three digit column value is the true number of these rils
colsha_ril_phenotypes$original.number
colsha_ril_phenotypes$genotype <- sub("(13RV)(00)(\\d)","RIL_\\3", colsha_ril_phenotypes$original.number)
colsha_ril_phenotypes$genotype <- sub("(13RV)(0)(\\d)","RIL_\\3", colsha_ril_phenotypes$genotype)
colsha_ril_phenotypes$genotype <- sub("(13RV)(\\d+)","RIL_\\2", colsha_ril_phenotypes$genotype)
colsha_ril_phenotypes$genotype

GWAS_phenotypes <- env_wght_pos_flr_srt[1431:length(env_wght_pos_flr_srt),]

setwd("~/git.repos/NAM/output/")

write.table(colsha_ril_phenotypes, file="col_sha_ril_phenotypes.csv", sep = ",")
write.table(GWAS_phenotypes, file=".csv", sep = ",")


#############
#############
#############
#blups
#############
#############
#############
?read.table
?naomit
colsha  <- read.table("col_sha_ril_phenotypes.csv", 
                      header=TRUE, sep = ",", na.strings = c("-","NA"),
                      stringsAsFactors = FALSE)
head(colsha)

#add shelf 
colsha$shelf <- sub("(3)(\\w)(\\d)(\\d+)","\\3", colsha$flatbarcode)
colsha$shelf
str(colsha)
dim(colsha)
colsha <- colsha[,c(1,21,2:20)]
head(colsha)


#for shortened name
str(colsha)
colsha[,1:10] <- lapply(colsha[,1:10], as.factor)
str(colsha)
dim(colsha)

colsha[,11:21] <- lapply(colsha[,11:21], as.numeric)
head(colsha)
tail(colsha)
str(colsha)

#set H for the reference light level
colsha$trt <- relevel(colsha$trt, ref = "H")

#wet dry weight correlations
colsha <- na.omit(colsha)
names(colsha)
str(colsha)
plot(colsha$biomass_dry_leaf, colsha$biomass_wet_leaf)
?cor
cor(colsha$biomass_dry_leaf, colsha$biomass_wet_leaf)
#s[1] 0.870805

plot(colsha$biomass_dry_reproductive, colsha$biomass_wet_reproductive)
?cor
cor(colsha$biomass_dry_reproductive, colsha$biomass_wet_reproductive)
#[1] 0.9372184
########
########
########
#plot phenotypes
########
########
########

dfplot <- function(data.frame)
{
  df <- data.frame
  ln <- length(names(data.frame))
  for(i in 1:ln){
    mname <- substitute(df[,i])
      if(is.factor(df[,i])){
        plot(df[,i],main=names(df)[i])}
        else{hist(df[,i],main=names(df)[i])}
  }
}

dfplot(colsha)
colnames(colsha)
#compare transformed and untransformed data in analysis
# colsha$biomass_dry_leaf_log <- -log(colsha$biomass_dry_leaf)


#######
#######
#Stats Models
#######
#######
library(lme4)
#library(lmerTest)
names(colsha)

drybiomass_model1 <- lmer(biomass_dry_leaf ~ trt + (1|shelf) +(0+trt|genotype), 
                          data = colsha, REML = FALSE)
summary(drybiomass_model1)

drybiomass_model2 <- lmer(biomass_dry_leaf ~ (1|shelf)+ (1|genotype),
                          data = colsha, REML = FALSE)
anova(drybiomass_model1,drybiomass_model2)


drybiomass_model3 <- lmer(biomass_dry_leaf ~ (1|trt) + (1|shelf),
                          data = colsha, REML = FALSE)
anova(drybiomass_model1,drybiomass_model3)


drybiomass_model4 <- lmer(biomass_dry_leaf ~ trt + (0+trt|genotype),
                          data = colsha, REML = FALSE)
anova(drybiomass_model1,drybiomass_model4)

#floweringtime
boltdays_model1 <- lmer(boltdays ~ trt + (1|shelf) + (0+trt|genotype), 
                          data = colsha, REML = FALSE)
summary(boltdays_model1)
boltdays_model2 <- lmer(boltdays ~ (1|shelf)+ (1|genotype),
                          data = colsha, REML = FALSE)
anova(boltdays_model1,boltdays_model2)


boltdays_model3 <- lmer(boltdays ~ trt + (1|shelf),
                          data = colsha, REML = FALSE)
anova(boltdays_model1,boltdays_model3)


boltdays_model4 <- lmer(boltdays ~ trt + (0+trt|genotype),
                          data = colsha, REML = FALSE)
anova(boltdays_model1,boltdays_model4)




#overwhelming support for genotype effects
drybiomass_model1_log <- lmer(biomass_dry_leaf_log ~ (0+trt|genotype),
                              data = colsha, REML = FALSE)
summary(drybiomass_model1_log)
drybiomass_model2_log <- lmer(biomass_dry_leaf_log ~ (1|genotype),
                              data = colsha, REML = FALSE)
anova(drybiomass_model1_log,drybiomass_model2_log)

names(colsha)
varlist <- names(colsha)[11:ncol(colsha)]
varlist

models <- lapply(varlist, function(x) {
    lmer(substitute(i ~ trt + (1|shelf) + (1+trt|genotype), list(i = as.name(x))), 
                      data = colsha)
})

#take a look
str(models)

models[[1]]
models[[2]]

#name the model list
names(models) <- varlist

#print a few traits to take a look
coef(models[[1]])$genotype
coef(models[[2]])$genotype

varlist

#extract the blups!!!!
#this works, but could be cleaned up
#it first extracts random intercept (mean) of each RIL
#it then extracts the random slope of the crowded treatment
rils <- data.frame(RILs = "")
for (trait in varlist) {
  rils <- merge(rils, data.frame(RILs=rownames(coef(models[[trait]])$genotype),                          
                                  placeholder=coef(models[[trait]])$genotype[,1]),
                all.y = T)
  colnames(rils)[length(colnames(rils))] <- trait
  rils <- merge(rils, data.frame(RILs=rownames(coef(models[[trait]])$genotype),                          
                                  placeholder=coef(models[[trait]])$genotype[,2]),
                all.y = T)
  colnames(rils)[length(colnames(rils))] <- paste(trait, "_shade", sep = "")
}
head(rils)
head(rils)
tail(rils)

rils$RILs <- as.character(rils$RILs)
str(rils)

head(rils)
tail(rils)

#get deviance value per average line for each trait
#make sure to only run this set of code once or you will have to remake the rils object as above
#can map this or map the estimates directly, RQTL does not mind and you get similar results
# rils$biomass_dry_leaf_shade <- rils$biomass_dry_leaf + rils$biomass_dry_leaf_shade
# rils$biomass_wet_leaf_shade <- rils$biomass_wet_leaf + rils$biomass_wet_leaf_shade
# rils$biomass_wet_reproductive_shade <- rils$biomass_wet_reproductive + rils$biomass_wet_reproductive_shade
# rils$biomass_dry_reproductive_shade <- rils$biomass_dry_reproductive + rils$biomass_dry_reproductive_shade
# rils$height1_shade <- rils$height1 + rils$height1_shade
# rils$height2_shade <- rils$height2 + rils$height2_shade
# rils$boltdays_shade <- rils$boltdays + rils$boltdays_shade
# rils$h1h2_shade <- rils$h1h2 + rils$h1h2_shade
# rils$h2h3_shade <- rils$h2h3 + rils$h2h3_shade
# rils$h1h3_shade <- rils$h1h3 + rils$h1h3_shade


setwd("~/git.repos/NAM/output/")                
write.table(rils,file="colsha_blups.csv",sep=",",row.names=FALSE) 

head(rils)

#format for RQTL
rils.t <- as.data.frame(t(rils))
head(rils.t)
colnames(rils.t)
dim(rils.t)
rils.t[24,] <- rils.t[1,]
rownames(rils.t)[24] <- "id"
rownames(rils.t)
rils.t <- rils.t[-1,]
head(rils.t)
tail(rils.t)

setwd("~/git.repos/NAM/output/")   
write.table(rils.t,file="colsha_blups_rqtl_full.csv",sep=",",row.names=TRUE, col.names = FALSE) 


#RQTL mapping of traits
library(qtl)
?read.cross
colsha_traits <- read.cross("csvsr", genfile ="colsha_map.csv", 
                         phefile="colsha_blups_rqtl_full.csv", genotypes=c("AA","BB"), na.strings = "NA")
head(colsha_traits)

class(colsha_traits)[1] <- "riself"
colsha_traits <- jittermap(colsha_traits)
colsha_traits

colsha_traits <- est.rf(colsha_traits)
plot.rf(colsha_traits) 
colsha_traits
newmap <- est.map(colsha_traits,verbose=T,error.prob=.01)

plot.map(colsha_traits,newmap) #some compression in this colsha_traits set
colsha_traits

colsha_traits <- replace.map(colsha_traits,newmap) #use new map
plot(colsha_traits) 
colsha_traits

colsha_traits <- sim.geno(colsha_traits,step=1,n.draws=64) 
#uses imputation
  #creates n.draw number of populations where the missing 
  #genotypes are filled in based on recombination frequencies

?calc.genoprob()
#just calculates the probabilities at different locations
colsha_traits <- calc.genoprob(colsha_traits, step = 1)
colsha_traits <- calc.genoprob(colsha_traits, error.prob=0.01)


summary(colsha_traits)

######biomass
scanone.perm.imp.1 <- scanone(colsha_traits, method = "imp", pheno.col = 1, n.perm = 1000) 
summary(scanone.perm.imp.1) 

perm95.1 <- summary(scanone.perm.imp.1)[1]
perm95.1

scanone.imp.1 <- scanone(colsha_traits, pheno.col = 1, method = "imp")
plot(scanone.imp.1, bandcol = "gray90")
abline(h=perm95.1,lty=2) 


######biomass_shade
scanone.perm.imp.2 <- scanone(colsha_traits, method = "imp", pheno.col = 2, n.perm = 1000) 
summary(scanone.perm.imp.2) 

perm95.2 <- summary(scanone.perm.imp.2)[1]
perm95.2

scanone.imp.2 <- scanone(colsha_traits, pheno.col = 2, method = "imp")
plot(scanone.imp.2, bandcol = "gray90")
abline(h=perm95.2,lty=2) 

######biomass_wet
scanone.perm.imp.3 <- scanone(colsha_traits, method = "imp", pheno.col = 3, n.perm = 1000) 
summary(scanone.perm.imp.3) 

perm95.3 <- summary(scanone.perm.imp.3)[1]
perm95.3

scanone.imp.3 <- scanone(colsha_traits, pheno.col = 3, method = "imp")
plot(scanone.imp.3, bandcol = "gray90")
abline(h=perm95.3,lty=2) 


######repro_biomass_dry
scanone.perm.imp.5 <- scanone(colsha_traits, method = "imp", pheno.col = 5, n.perm = 1000) 
summary(scanone.perm.imp.5) 

perm95.5 <- summary(scanone.perm.imp.5)[1]
perm95.5

scanone.imp.5 <- scanone(colsha_traits, pheno.col = 5, method = "imp")
plot(scanone.imp.5, bandcol = "gray90")
abline(h=perm95.5,lty=2) 



######repro_biomass_wet
scanone.perm.imp.7 <- scanone(colsha_traits, method = "imp", pheno.col = 7, n.perm = 1000) 
summary(scanone.perm.imp.7) 

perm95.7 <- summary(scanone.perm.imp.7)[1]
perm95.7

scanone.imp.7 <- scanone(colsha_traits, pheno.col = 7, method = "imp")
plot(scanone.imp.7, bandcol = "gray90")
abline(h=perm95.7,lty=2) 




######height_2
scanone.perm.imp.11 <- scanone(colsha_traits, method = "imp", pheno.col = 11, n.perm = 1000) 
summary(scanone.perm.imp.11) 

perm95.11 <- summary(scanone.perm.imp.2)[1]
perm95.11

scanone.imp.11 <- scanone(colsha_traits, pheno.col = 11, method = "imp")
plot(scanone.imp.11, bandcol = "gray90")
abline(h=perm95.11,lty=2)



######height_shade_2
scanone.perm.imp.12 <- scanone(colsha_traits, method = "imp", pheno.col = 12, n.perm = 1000) 
summary(scanone.perm.imp.12) 

perm95.12 <- summary(scanone.perm.imp.2)[1]
perm95.12

scanone.imp.12 <- scanone(colsha_traits, pheno.col = 12, method = "imp")
plot(scanone.imp.12, bandcol = "gray90")
abline(h=perm95.12,lty=2)


######height_3
scanone.perm.imp.13 <- scanone(colsha_traits, method = "imp", pheno.col = 13, n.perm = 1000) 
summary(scanone.perm.imp.13) 

perm95.13 <- summary(scanone.perm.imp.2)[1]
perm95.13

scanone.imp.13 <- scanone(colsha_traits, pheno.col = 13, method = "imp")
plot(scanone.imp.13, bandcol = "gray90")
abline(h=perm95.13,lty=2)


######height_shade_3
scanone.perm.imp.14 <- scanone(colsha_traits, method = "imp", pheno.col = 14, n.perm = 1000) 
summary(scanone.perm.imp.14) 

perm95.14 <- summary(scanone.perm.imp.2)[1]
perm95.14

scanone.imp.14 <- scanone(colsha_traits, pheno.col = 14, method = "imp")
plot(scanone.imp.14, bandcol = "gray90")
abline(h=perm95.14,lty=2)


######bolt
scanone.perm.imp.15 <- scanone(colsha_traits, method = "imp", pheno.col = 15, n.perm = 1000) 
summary(scanone.perm.imp.15) 

perm95.15 <- summary(scanone.perm.imp.2)[1]
perm95.15

scanone.imp.15 <- scanone(colsha_traits, pheno.col = 15, method = "imp")
plot(scanone.imp.15, bandcol = "gray90")
abline(h=perm95.14,lty=2)

######bolt_shade
scanone.perm.imp.15 <- scanone(colsha_traits, method = "imp", pheno.col = 16, n.perm = 1000) 
summary(scanone.perm.imp.16) 

perm95.16 <- summary(scanone.perm.imp.2)[1]
perm95.16

scanone.imp.16 <- scanone(colsha_traits, pheno.col = 16, method = "imp")
plot(scanone.imp.16, bandcol = "gray90")
abline(h=perm95.14,lty=2)

######h1h2
scanone.perm.imp.17 <- scanone(colsha_traits, method = "imp", pheno.col = 17, n.perm = 1000) 
summary(scanone.perm.imp.17) 

perm95.17 <- summary(scanone.perm.imp.2)[1]
perm95.17

scanone.imp.17 <- scanone(colsha_traits, pheno.col = 17, method = "imp")
plot(scanone.imp.17, bandcol = "gray90")
abline(h=perm95.14,lty=2)


######h2h3
scanone.perm.imp.19 <- scanone(colsha_traits, method = "imp", pheno.col = 19, n.perm = 1000) 
summary(scanone.perm.imp.19) 

perm95.19 <- summary(scanone.perm.imp.2)[1]
perm95.19

scanone.imp.19 <- scanone(colsha_traits, pheno.col = 19, method = "imp")
plot(scanone.imp.19, bandcol = "gray90")
abline(h=perm95.14,lty=2)


######h1h3
scanone.perm.imp.21 <- scanone(colsha_traits, method = "imp", pheno.col = 21, n.perm = 1000) 
summary(scanone.perm.imp.21) 

perm95.21 <- summary(scanone.perm.imp.2)[1]
perm95.21

scanone.imp.21 <- scanone(colsha_traits, pheno.col = 21, method = "imp")
plot(scanone.imp.21, bandcol = "gray90")
abline(h=perm95.21,lty=2)

#### Get the scanone output for Tiffany
str(scanone.imp.21)

######### QUICK DATA FOR TIFFANY ######
scanone.mr.7 <- scanone(colsha_traits, pheno.col = 7, method = "mr")
scanone.mr.4 <- scanone(colsha_traits, pheno.col = 4, method = "mr")
scanone.mr.5 <- scanone(colsha_traits, pheno.col = 5, method = "mr")


traits <- as.data.frame(scanone.mr.7$chr)
traits$pos <- scanone.mr.7$pos
traits$trait1_lod <- scanone.mr.7$lod
traits$trait2_lod <- scanone.mr.4$lod
traits$trait3_lod <- scanone.mr.5$lod
str(traits) 
str(colsha_traits)
write.table(traits, "traits_for_tiffany.csv", row.names = FALSE, col.names = TRUE, sep = ",")








