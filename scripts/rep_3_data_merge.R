library(reshape2)
library(ggplot2)

# change working directory for your machine -----------------
setwd("/Users/Cody_2/Box Sync/NAM_CAM/data")

?read.csv
imageData <- read.csv("Rep 3 Image Analysis.csv", header = TRUE, na.strings = "NA")
genoData <- read.csv("plant_id_position_col_sha_rep3_clean.csv", header = TRUE, na.strings = "NA")
WWRData <- read.csv("Rep 3 WWR Schedule.csv", header = TRUE)

# exploring imported files -----------------
head(imageData)
str(imageData)
tail(imageData)

hist(imageData$Area)
hist(imageData$Perimeter)

head(genoData)
str(genoData)

# add column names  -----------------
colnames(WWRData) <- c("Flat_Number", "20131108", "20131113")

# made copies of dataframes to make sure loops work
imgData <- imageData
geno <- genoData
WWR <- melt(WWRData)
head(WWR)

# the following loops are written by James and are hardcoded for the position they
# are refering to
for (i in 1:length(WWR$Flat_Number)) {
  j <- substr(WWR$Flat_Number[i], 3, 3)
  WWR$Flat_ID[i] <- paste(j, WWR$value[i], WWR$variable[i], sep = "_")
}

for (i in 1:length(imgData$Tray_Number)) {
  j <- substr(imgData$Image_Name[i], 18, 18)
  imgData$Flat_ID[i] <- paste(j, imgData$Tray_Number[i], imgData$Date[i], sep = "_")
}

head(imgData)
tail(imgData)
dim(imgData)

head(WWR)
tail(WWR)
dim(WWR)

# Merge the WWR and ImgData with the new flat ID column in each
image <- merge(WWR, imgData, by = "Flat_ID")
dim(image)
head(image)

# comment
for (i in 1:length(image$Flat_ID)) {
  j <- as.numeric(substr(image$Flat_Number[i], 3, 3))
  z <- as.numeric(substr(image$Flat_Number[i], 5, 6))
  image$ID[i] <- paste(j, z, image$Plant_Position[i], sep = "_")
}
head(image)
tail(image)

head(geno)
str(geno)

for (i in 1:length(geno$Flat.number)) {
  k <- as.numeric(substr(geno$Flat.number[i], 3, 3))
  j <- as.numeric(substr(geno$Flat.number[i], 5, 6))
  geno$ID[i] <- paste(k, j, geno$Position.in.flat[i], sep = "_")
}

dim(geno)
dim(image)


geno_image <- merge(geno, image, by = "ID")
dim(geno_image)
head(geno_image, 50)

geno_image[grep("c1|c2|c3", geno_image$Image_Name), "Treatment"] <- "Shade"
geno_image[grep("c4|c5|c6", geno_image$Image_Name), "Treatment"] <- "Sun"
head(geno_image)

names(geno_image)
geno_image <- geno_image[,-c(19)]
names(geno_image)


# infile the other phenotype data
setwd("~/git.repos/NAM/output/")
env_wght_pos_flr  <- read.table("col_sha_NAM_GWAS_phenotypes.csv", sep = ",", header = TRUE)
head(env_wght_pos_flr)
tail(env_wght_pos_flr)

hist(env_wght_pos_flr$h2h3)
hist(env_wght_pos_flr$height1)

str(env_wght_pos_flr)

env_wght_pos_flr$height1

is.numeric(env_wght_pos_flr$biomass_wet_reproductive)

head(geno_image)
head(env_wght_pos_flr)
dim(geno_image)
unique(geno_image$Plant.Barcode)

# change column order of geno_image
names(geno_image)
geno_image2 <- geno_image[,c(3, 7, 8, 10, 13, 14, 15, 17, 18, 19, 20)]
head(geno_image2)

dim(geno_image2)
dim(env_wght_pos_flr)

# there are two dates in the data set and some of the column data is duplicated
# subset these and merge the two subsetted geno_image2 dataframes with env_wght_pos_flr
geno_image_fd <- subset(geno_image2, Date == "20131113")
head(geno_image_fd)
date_sd <- subset(geno_image2, Date == "20131108")
head(date_sd)

dim(env_wght_pos_flr)
dim(date_sd)
combined_1 <- merge(date_sd, env_wght_pos_flr, by = c("Plant.Barcode", "original.number"), all.x = TRUE )
dim(combined_1)
head(combined_1)
tail(combined_1)
names(combined_1)
comb_1_clean <- combined_1[,c(1, 2, 4, 6, 7, 10, 18, 8, 9, 19:29)]
head(comb_1_clean)
# rename columns
colnames(comb_1_clean)[8] <- "Area_20131108"
colnames(comb_1_clean)[9] <- "Perimeter_20131108"
dim(comb_1_clean)

# merge the other subsetted dataframe with our combined dataframe
combined_2 <- merge(geno_image_fd, comb_1_clean, by = c("Plant.Barcode", "original.number"), all.y = TRUE )
dim(combined_2)
head(combined_2)
names(combined_2)
combined_2 <- combined_2[, -c(3, 5, 6, 11, 12, 13, 14, 15)]
head(combined_2)
names(combined_2)

combined_2 <- combined_2[,c(1:4, 7, 8, 5, 6, 9:21)]
head(combined_2)
dim(combined)
head(combined)
colnames(combined_2)[7] <- "Area_20131113"
colnames(combined_2)[8] <- "Perimeter_20131113"
head(combined_2)
tail(combined_2)
colnames(combined_2)[3] <- "Flat_Number"
colnames(combined_2)[4] <- "Plant_Position"
colnames(combined_2)[5] <- "Treatment"
combined_2$Treatment
combined_2$Plant.Barcode


# clean up and reorder to make stats modeling easier
combined_3 <- combined_2
combined_3$treatment <- sub("(\\d)(\\w)(\\d)(\\d)(\\d)(\\d)","\\2", combined_2$Plant.Barcode)
names(combined_3)
combined_3 <- combined_3[,c(1:4,22,6:21)]
head(combined_3)
combined_3$treatment[combined_3$treatment == "H"] <- "Sun"
combined_3$treatment[combined_3$treatment == "L"] <- "Shade"
head(combined_3)

setwd("/Users/Cody_2/Box Sync/NAM_CAM/data")
write.table(combined_3, "col_sha_biomass_area_flr_combined.csv", sep = ",", row.names = FALSE)



