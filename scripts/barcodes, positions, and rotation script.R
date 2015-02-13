setwd(dir="/Users/christine/Desktop/")

##Not sure why, but I called each shelf "rep1, rep2". It's confusing, but just don't think of rep as replicate and you'll be fine
##change the first number to the round number.  For example, the next round (july 2014) could be number round 4.  
##round 1 was full nam set, some troubles and incomplete round
##round 2 was full nam set, complete round, fully awesome
##round 3 was Colxsha set, 6 replicates, fully awesome

##so this would be rep1<-rep(#4H40)
##first digit, round #
##second digit, trt (H=sun, L=shade ("low light")
##third digit, shelf number
##fourth digit, spacer "0"
##fifth and sixth digit, flat number
rep1<-rep("3H40", 24)
numbers<-paste(0,1:9, sep="")
rep1flat<-c(numbers, 10:24)
rep1final<-paste(rep1, rep1flat, sep="")


rep2<-rep("3H50", 24)
numbers<-paste(0,1:9, sep="")
rep1flat<-c(numbers, 10:24)
rep2final<-paste(rep2, rep1flat, sep="")
rep1final

rep3<-rep("3H60", 24)
numbers<-paste(0,1:9, sep="")
rep1flat<-c(numbers, 10:24)
rep3final<-paste(rep3, rep1flat, sep="")
rep1final

rep4<-rep("3L10", 24)
numbers<-paste(0,1:9, sep="")
rep1flat<-c(numbers, 10:24)
rep4final<-paste(rep4, rep1flat, sep="")
rep4final


rep5<-rep("3L20", 24)
numbers<-paste(0,1:9, sep="")
rep1flat<-c(numbers, 10:24)
rep5final<-paste(rep5, rep1flat, sep="")
rep5final

rep6<-rep("3L30", 24)
numbers<-paste(0,1:9, sep="")
rep1flat<-c(numbers, 10:24)
rep6final<-paste(rep6, rep1flat, sep="")

###This makes your complete list of flat barcodes.
rep1total<-c(rep4final, rep5final, rep6final, rep1final, rep2final, rep3final)
rep1total

###TO print flat barcodes
write.csv(rep1total, file = "rep1total.csv")

##This makes your complete list of flats and plant positions, so that you can later assign a genotype to each position
rep3totalforplantbarcode<-cbind(c(rep(rep4final,each=16),rep(rep5final, each=16),rep(rep6final,each=16),rep(rep1final,each=16),
                            rep(rep2final,each=16),rep(rep3final,each=16)),rep(1:16,144))
write.csv(rep3totalforplantbarcode, "~/Desktop/allflatposrep3.csv")


##This part creates the randomized layout of the genotypes in the flats, not sure why I recreated the positions, but whatever
##assign a number 1-n (where n is the number of genotypes) and sample
##paste to complete list of flats/positions from above
plantpos<-c(rep1final, rep2final, rep3final)
plantposs<-rep(plantpos, 16)
plantpossort<-sort(plantposs)
inflat<-rep(1:16, 72)
plantposortogether<-cbind(plantpossort,inflat)
genotype<-sample(1:1152)
plantpositions<-cbind(plantposortogether,genotype)
plantpositions
summary(plantpositions)

### This is the megamaster layout file, that tells you where each plant is.  I print this out and use it to arrange the tubes after seeds are added  
##You'll want a separate sheet with what genotype number (1-n) goes with which genotype (RILname)
write.csv(plantpositions, file = "plantpositions.csv")



##just make a randomized set for the special round, with 384 genotypes
numbers<-rep(1:384)
shelf1<-sample(numbers, replace=F)
shelf2<-sample(numbers, replace=F)
shelf3<-sample(numbers, replace=F)

shelves<-cbind(shelf1,shelf2,shelf3)
write.csv(shelves, "~/Desktop/genotyperandomizationrep3.csv")


###This is to generate positions for the flats on the shelves for each rotation.  Each line is one day of randomizing.
##This generates 3 shelves of random positions for the sun side, mirrored on the shade.  
##I did  6 days worth (2 weeks, since it is MWF, MWF) at a time and pasted into my googlesheet online,
##Google sheet had flat numbers, and just pasted these positions next to it.
randomizer<-c(sample (1:24), sample (1:24), sample (1:24))
randomizer2<-c(sample (1:24), sample (1:24), sample (1:24))
randomizer3<-c(sample (1:24), sample (1:24), sample (1:24))
randomizer4<-c(sample (1:24), sample (1:24), sample (1:24))
randomizer5<-c(sample (1:24), sample (1:24), sample (1:24))
randomizer6<-c(sample (1:24), sample (1:24), sample (1:24))

twoweek<-cbind(randomizer,randomizer2, randomizer3, randomizer4, randomizer5, randomizer6)
write.csv(twoweek, file = "~/Desktop/twoweek.csv")


# rep1<-rep("1H40", 19)
# numbers<-paste(0,1:9, sep="")
# rep1flat<-c(numbers, 10:19)
# rep1final<-paste(rep1, rep1flat, sep="")
# 
# 
# rep2<-rep("1H50", 24)
# numbers<-paste(0,1:9, sep="")
# rep1flat<-c(numbers, 10:24)
# rep2final<-paste(rep2, rep1flat, sep="")
# rep1final
# 
# rep3<-rep("1H60", 24)
# numbers<-paste(0,1:9, sep="")
# rep1flat<-c(numbers, 10:24)
# rep3final<-paste(rep3, rep1flat, sep="")
# rep1final

##This part creates the randomized layout of the genotypes in the flats
##assign a number 1-n (where n is the number of genotypes) and sample
##paste to complete list of flats/positions from above
plantpos<-c(rep1final, rep2final, rep3final)
plantposs<-rep(plantpos, 16)
plantpossort<-sort(plantposs)
inflat<-rep(1:16, 72)
plantposortogether<-cbind(plantpossort,inflat)
genotype<-sample(1:1152)
plantpositions<-cbind(plantposortogether,genotype)
plantpositions
summary(plantpositions)
write.csv(plantpositions, file = "plantpositions.csv")



##just make a randomized set for the special round, with 384 genotypes
numbers<-rep(1:384)
shelf1<-sample(numbers, replace=F)
shelf2<-sample(numbers, replace=F)
shelf3<-sample(numbers, replace=F)

shelves<-cbind(shelf1,shelf2,shelf3)
write.csv(shelves, "~/Desktop/genotyperandomizationrep3.csv")