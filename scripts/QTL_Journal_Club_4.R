#load the qtl mapping library
library(qtl)
setwd("~/Documents/Classes/QTL Mapping/")

hyp.qtl <- read.cross(format="csv",file="Cvl_genotypes_reformatted.csv",
                        genotypes=c("A","B"))
head(hyp.qtl)
class(hyp.qtl)[1] <- "riself"
summary(hyp.qtl) #ap, phenotypes, number of reps/line
plot(hyp.qtl) 

hyp.qtl <- est.rf(hyp.qtl)
plot.rf(hyp.qtl) 

newmap <- est.map(hyp.qtl,verbose=T,error.prob=.01)
plot.map(hyp.qtl,newmap) 
hyp.qtl <- replace.map(hyp.qtl,newmap) #use new map
plot.map(hyp.qtl)

hyp.qtl <- calc.errorlod(hyp.qtl, error.prob=0.01)
top.errorlod(hyp.qtl)
plot.errorlod(hyp.qtl)

#next we calculate genotype probabilities at "virtual" 
#loci between our real markers
?sim.geno()
hyp.qtl <- sim.geno(hyp.qtl,step=1,n.draws=32) 
#uses imputation
  #creates n.draw number of populations where the missing 
  #genotypes are filled in based on recombination frequencies

?calc.genoprob()
#just calculates the probabilities at different locations
hyp.qtl <- calc.genoprob(hyp.qtl, step = 1) 

#####MAPPING FINALLY!!!!########

###marker regression---ANOVA at each marker
out.mr <- scanone(hyp.qtl, method="mr")
out.mr[out.mr$chr==1, ] 
summary(out.mr, threshold=3)
max(out.mr)  ####Find the largest LOD score
plot(out.mr, chr=c(1, 2, 3, 4, 5 ), ylab="LOD Score")  ###the jagged appearnce is from missing marker data

### Interval Mapping
hyp.qtl <- calc.genoprob(hyp.qtl, step=1, error.prob=0.001)
hyp.qtl <- jittermap(hyp.qtl)  ####incase markers are at the same position

out.em <- scanone(hyp.qtl, method="em")
plot(out.em, ylab="LOD score")

plot(out.em, out.mr, chr=c(1, 2, 3, 4, 5), col=c("blue", "red"), ylab="LOD score")
plot(out.em, chr=c(1, 2, 3, 4, 5), col="blue", ylab="LOD score")
plot(out.mr, chr=c(1, 2, 3, 4, 5), col="red", add=TRUE)
plot(out.em, out.mr, chr=c(1, 2, 3, 4, 5), col="black", lty=1:2, ylab="LOD score")


### Haley-KNot Regression
out.hk <- scanone(hyp.qtl, method="hk")
plot(out.em, out.hk, chr=c(1, 2, 3, 4, 5), col=c("blue","red"),ylab="LOD score")


plot(out.hk - out.em, chr=c(1, 2, 3, 4, 5), ylim=c(-0.5, 1.0), ylab=expression(LOD[HK] - LOD[EM]))
abline(h=0, lty=3)

plot(rnorm(100), rnorm(100), xlab=expression(hat(mu)[0]), ylab=expression(alpha^beta), main=expression(paste("Plot of ", alpha^beta, " versus ", hat(mu)[0])))

##### Extended Haley_Knot Regression- this also considers the variances
out.ehk <- scanone(hyp.qtl, method="ehk")


#####To more clearly see the differences among results plot the differences in the LOD scores
plot(out.hk - out.em, out.ehk - out.em, chr=c(1, 2, 3, 4, 5),col=c("blue", "red"), ylim=c(-0.5, 1),ylab=expression(LOD[HK] - LOD[EM]))
abline(h=0, lty=3)

####Multiple Imputation Models ---Better for 2 QTL scans
#### Can change the step size and n.draws to make the computation time less
hyp.qtl <- sim.geno(hyp.qtl, step=1, n.draws=64, error.prob=0.001)
out.imp <- scanone(hyp.qtl, method="imp")
plot(out.em, out.imp, chr=c(1, 2, 3, 4, 5), col=c("blue", "red"),ylab="LOD score")


plot(out.imp - out.em, chr=c(1, 2, 3, 4, 5), ylim=c(-0.5, 0.5), ylab=expression(LOD[IMP] - LOD[EM]))


















#start by doing a one dimensional scan on imputed genotypes

# RQTL 4.2.4
#first permute the data to set a signficance threshold
#randomizes all phenotype data with genotype data to set backgroud LOD scores
scanone.perm.imp <- scanone(hyp.qtl, method = "imp", n.perm = 1000) 
#takes ~1 minute

#look at a summary of our permutations
#shows the 90th and 95th percentile of LOD scores in randomized data
summary(scanone.perm.imp) 
perm95 <- summary(scanone.perm.imp)[1] #keep the 95th percentile for future use.
                                       #This corresponds to p <0.05


?as.numeric
?pull.pheno
?scanone
pull.pheno(hyp.qtl, "Cytoplasm")

#Turn cytoplasm into a numeric vector to use a covariate
cyto <- as.numeric(pull.pheno(hyp.qtl, "Cytoplasm"))
cyto

#note method "em" = interval mapping; "mr" = marker regression
#here we use method "imp" for multiple imputation.
  #advantage of "imp" is better for missing genotype data
  #disadvantage is that it is slower.  Not a big deal on current computers.
  #See table 4.2 in RQTL book page 102
scanone.imp <- scanone(hyp.qtl, pheno.col = 4, method = "imp")
scanone.imp.cov <- scanone(hyp.qtl, addcovar = cyto, 
                            pheno.col = 4, method = "imp")

plot(scanone.imp, bandcol = "gray90")
plot(scanone.imp.cov, bandcol = "gray90")
abline(h=perm95,lty=2) #add permuation threshold

#plot data and data with cytoplasm covariate
plot(scanone.imp, scanone.imp.cov, col=c("blue", "red"), 
       lty=1:2, ylab="LOD score", alternate.chrid=TRUE)
abline(h = perm95, lty = 2)
# the peaks on the chromosomes look almost identical


#Exercise 5: How many chromosomes have significant QTL?
#obtain a summary
summary(scanone.imp, perm = scanone.perm.imp, alpha = .05)
#unfortunately the summary only prints out the top peak per chromosome, 
#even if there is more than one (i.e. chromosome one).


##do cim mapping ----Section 7.4 page 205 RQTL book:
#The advantages of the simultaneous consideration of multiple
#QTL are to (a) reduce residual variation and so better detect loci of more
#modest effect, (b) separate linked QTL, and (c) identify interactions among
#QTL.
#a more sophisticated model is to use the most significant markers as covariates
 #during a 1D scan.  This is known as composite interval mapping (CIM) 
 #and is implmented using the function cim

cim.qtl <- cim(hyp.qtl, n.marcovar = 5, method = "em")
  #here we use the interval mapping method "em" as this is how 
    #cim was originaly implemented.
  #the n.marcovar = argument defines the maximum number of marker 
    #covariates to use.

#run cim permutations.  You should do 1000, but I only do 100 here to save time 
#1,000 will take 15 minutes
cim.perm <- cim(hyp.qtl, n.marcovar = 5, method = "em", n.perm = 100) 
#takes a few minutes

#uncomment line below to run 1000 permutations instead of 100
#cim.perm <- cim(hyp.qtl,n.marcovar=5,method="em",n.perm=1000) 
#takes 15 minutes

summary(cim.perm)
cim.perm95 <- summary(cim.perm)[1] #keep 95% 

plot(cim.qtl,bandcol="gray90")
abline(h=cim.perm95,lty=2)

#Exercise 6
#How do the CIM results compare to the scaneone results?
#A: Interestingly there are now significant QTL on each chromosome when
# accounting for markers as covariates

summary(cim.qtl)
#for comparison, reprint the imp results
summary(scanone.imp)

##a 2D scan
#The above methods only evaluate the evidence for a QTL one location a time.  
#It is possible to instead scan two locations at once.  
#This can be particularly advantageous in identifying QTL if there is epistasis.

scantwo.imp <- scantwo(hyp.qtl,method="imp",verbose=T) 
#takes a couple of minutes

#permutations take a long time, about 2 hours for 1000
#I have pre-run them.
#uncomment this line to run the permutations yourself
#note that if you have a multi-core computer you can spread this out
  #across mutliple processors. To do this you need library(snow) and
  #add the argument n.clusters=8 (or replace 8 with the number of cores)
#scantwo.imp.perm <- scantwo(hyp.qtl,method="imp",verbose=T,n.perm=1000)) 

#load the pre-run permutations
load("scantwoperm.Rdata")
summary(scantwo.imp.perm)

#there are many different models being considered in the 2D models.  
#The LOD scores reflect the following comparisons
#columns:
  #full 2 QTL w/ interaction vs. null model
  #fv1 2 QTL w/ interaction vs 1 QTL
  #int 2 QTL w/interaction vs 2QTL additive.
  #add 2 QTL additive vs null model
  #av1 2 QTL w/interaction vs 1 QTL
  #one 1 QTL model vs null
######Diff between fv1 and av1????


summary(scantwo.imp,perm=scantwo.imp.perm,pvalues=T)

plot(scantwo.imp,zscale=T) 
      #upper compares model for 2 QTL with interaction to 2 QTL additive
          #(evidence for interaction); corresponds to "int" above
      #lower compares 2 QTL interaction with no QTL 
          #(evidence for QTL); corresponds to "full" above

plot(scantwo.imp,zscale=T,lower="fv1")
      #now lower compares 2QTL interaction vs 1QTL model
        #(evidence for a second QTL); corresponds to "fv1" above

plot(scantwo.imp,zscale=T,lower="fv1",upper="av1") #now upper is 
      #evidence for 2 QTL additive model vs 1 QTL;
        #corresopnds to "av1" above.
      #differences in upper and lower refelct possible interactions

#Excercise 7
#Do you see evidence for interaction between QTLs?  
#You should consult the scantwo summary and the plots.

#The next step is to fit a multiple QTL model, where all QTL are accounted 
#for at the same time.  We use the previous results as a starting point.  
#We can then ask if any QTLs should be dropped or added.

#review previous results
summary(scanone.imp)
summary(cim.qtl)
summary(scantwo.imp,perm=scantwo.imp.perm,alpha=.05,pval=T)


#create an object to hold QTL of interest
qtlm <- makeqtl(hyp.qtl,chr=c(1,1,2,4,5),pos=c(25,115.6,34,46,14))

qtlm

#make a model for QTL action.  Since there is no evidence of interaction, 
#start with an additive model
qtl.fit <- fitqtl(cross=hyp.qtl,qtl=qtlm,
                   formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5, get.ests = T)

#examine our model
summary(qtl.fit)

#Exercise 8
# Are all of the QTL well supported?  
# If not, make simpler model using the dropfromqtl() function 

qtlm <- dropfromqtl(qtlm, qtl.name = "Q5")
#alternate form
#qtlm <- dropfromqtl(qtlm,chr=5,pos=14)

#fit a new model 
qtl.fit <- fitqtl(cross=hyp.qtl,qtl=qtlm,formula = y ~ Q1 + Q2 + Q3 + Q4)

summary(qtl.fit)

#we can now ask if all qtls are at their most likely positions
#this also provides a nice map
qtl.refine <- refineqtl(hyp.qtl, qtl = qtlm,
                        formula= y ~ Q1 + Q2 + Q3 + Q4)
#you can see from the LOD scores that the new model 
#is not much better than the old model.  
#Nevertheless we can keep this as our new model

#positions can be plotted on the map
plot(qtl.refine)
#or showing lod scores
plotLodProfile(qtl.refine)

#we can test for interactions specifically among these loci
addint(hyp.qtl, qtl = qtl.refine) # no evidence for interactions

#we can scan the genome for additional QTL to add to this model
qtl.add <- addqtl(hyp.qtl,qtl=qtl.refine)
summary(qtl.add)
plot(qtl.add)
#Exercise 9.
#Do you think any of these should be added?  If so why? 
#(Hint: review the scanone permuation threshold).

#Excercise 10
#If you think new QTL should be added, do so, using addtoqtl()
qtlm2 <- addtoqtl(hyp.qtl,qtl.refine,chr=c("1","5"),
                         pos=c(5,100))
summary(qtlm2)

#refine the positions of our new model
qtl.refine2 <- refineqtl(hyp.qtl,qtl=qtlm2,
                         formula = y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6)

#now look at a summary of the new model
summary(qtl.refine2)
plot(qtl.refine2)
plotLodProfile(qtl.refine2)

#at this point we will do a final fit and summary of that fit
#(although you could repeat looking for interactions)
qtl.fit2 <- fitqtl(hyp.qtl,qtl=qtl.refine2,get.ests=T)
summary(qtl.fit2)

#Excercise 11: are all of the QTL supported?
#is more phenotypic variance explained than in the first model?

#Exercise 12: What is the QTL with the largest effect?  
# Do all QTL act in the same direction?






