library(rstan)
library(rethinking)
setwd("~/git.repos/NAM/output/")
colsha  <- read.table("col_sha_ril_phenotypes.csv", 
                      header=TRUE, sep = ",", na.strings = c("-","NA"),
                      stringsAsFactors = FALSE)

head(colsha, 20)


colsha_m <- colsha[, c("biomass_wet_leaf", "trt", "genotype")]
head(colsha_m)
colnames(colsha_m) <- paste(c("bio", "trt", "geno"))
hist(colsha_m$bio)
plot(colsha_m$bio)
str(colsha_m)
dim(colsha_m)
?ifelse
# have to make STAN happy with variable coding
colsha_m$trt <- ifelse(colsha_m$trt == "H", 0, 1)
colsha_m$geno <- as.numeric(sub("(RIL_)(\\d+)","\\2", colsha_m$geno))
plot(colsha_m$bio)
colsha_m <- na.omit(colsha_m)
colsha_m <- colsha_m[with(colsha_m, order(geno, trt)),]

head(colsha_m, 20)
tail(colsha_m, 20)


#simple model
colsha_stan1 <- map2stan(
	alist(
		bio ~ dnorm(mu, sigma),
		mu <- a + bg*geno + bt*trt, 
		a ~ dnorm(0, 100), 
		bg ~ dnorm(0, 100),
		bt ~ dnorm(0, 10),
		sigma ~ dcauchy(0, 10)
	),
	data = colsha_m, iter = 1e4, chains = 2
)


#build in genotype
colsha_stan1 <- map2stan(
	alist(
		bio ~ dnorm(mu, sigma),
		mu <- a + bg*geno + bt*trt, 
		a ~ dnorm(0, 100), 
		bg ~ dnorm(0, 100),
		bt ~ dnorm(0, 10),
		sigma ~ dcauchy(0, 10)
	),
	data = colsha_m, iter = 1e4, chains = 2
)

#build in genotype
colsha_stan1 <- map2stan(
	alist(
		bio ~ dnorm(mu , sigma),
		mu <- a_geno[geno] + bg_geno[geno]*trt, 
		c(a_geno, bg_geno)[geno] ~ dmvnorm2( c(a, bg), sigma_geno, Rho),
		a ~ dnorm(0, 100), 
		bg ~ dnorm(0, 100),
		sigma_geno ~ dcauchy(0, 2),
		sigma ~ dcauchy(0 , 2),
		Rho ~ dlkjcorr(2)
	),
	data = colsha_m, warmup = 1000, iter = 5000, chains = 2
)

precis(colsha_stan1)
plot(colsha_stan1)
plot(precis(colsha_stan1))
coef(colsha_stan1)



head(colsha)
plot(colsha$boltdays)

colsha_bolt <- colsha[, c("boltdays", "trt", "genotype")]
head(colsha_bolt)
colnames(colsha_bolt) <- paste(c("bolt", "trt", "geno"))
hist(colsha_bolt$bolt)
str(colsha_bolt)
dim(colsha_bolt)
?ifelse
# have to make STAN happy with variable coding
colsha_bolt$trt <- ifelse(colsha_bolt$trt == "H", 0, 1)
colsha_bolt$geno <- as.numeric(sub("(RIL_)(\\d+)","\\2", colsha_bolt$geno))
colsha_bolt <- na.omit(colsha_bolt)
colsha_bolt <- colsha_bolt[with(colsha_bolt, order(geno, trt)),]
head(colsha_bolt)


colsha_stan3 <- map2stan(
	alist(
		bolt ~ dnorm(mu , sigma),
		mu <- a_geno[geno] + bg_geno[geno]*trt, 
		c(a_geno, bg_geno)[geno] ~ dmvnorm2( c(a, bg), sigma_geno, Rho),
		a ~ dnorm(0, 100), 
		bg ~ dnorm(0, 100),
		sigma_geno ~ dcauchy(0, 2),
		sigma ~ dcauchy(0 , 2),
		Rho ~ dlkjcorr(2)
	),
	data = colsha_bolt, warmup = 1000, iter = 5000, chains = 2
)

# this model looks great!
precis
precis(colsha_stan, depth = 2)
plot(colsha_stan3)
precis(colsha_stan3, depth = 2)
coef(colsha_stan3)
str(colsha_stan3)
bolting_out <- precis(colsha_stan3, depth = 2)
bolting_out <- bolting_out@output
head(bolting_out)
str(bolting_out)
bolting_out$name <- rownames(bolting_out)

# to do is figure out ggplot
colnames(bolting_out)[3:4] <- paste(c("low95", "high95"))
head(bolting_out) 
str(bolting_out)
bolting_out$name <- sub("(*)+(\\])", "\\1", bolting_out$name) 
bolting_out$name <- sub("(*)+(,)(*)+", "\\1\\3", bolting_out$name) 
head(bol)
tail(bolting_out, 20)
library(ggplot2)
?geom_pointrange
bolting_out$name

setwd("~/git.repos/NAM_CAM_CHAMBER/output/")
bolt_plot <- ggplot(data = bolting_out)
bolt_plot <- bolt_plot + geom_pointrange(mapping=aes(x = name, y = Mean, ymin=low95, ymax=high95)) + coord_flip()
bolt_plot
?ggsave
ggsave("bolt_hdi.pdf", bolt_plot, width = 10, height = 30)


stancode(colsha_stan3)


head(colsha)
colsha_h3 <- colsha[, c("height3", "trt", "genotype")]
head(colsha_h3)
colnames(colsha_h3) <- paste(c("height", "trt", "geno"))
hist(colsha_h3$height)
str(colsha_h3)
dim(colsha_h3)
?ifelse

# have to make STAN happy with variable coding
colsha_h3$trt <- ifelse(colsha_h3$trt == "H", 0, 1)
colsha_h3$geno <- as.numeric(sub("(RIL_)(\\d+)","\\2", colsha_h3$geno))
colsha_h3 <- na.omit(colsha_h3)
colsha_h3 <- colsha_h3[with(colsha_h3, order(geno, trt)),]



colsha_stan2 <- map2stan(
	alist(
		height ~ dnorm(mu, sigma),
		mu <- a_geno[geno] + bg_geno[geno]*trt, 
		c(a_geno, bg_geno)[geno] ~ dmvnorm2( c(a, bg), sigma_geno, Rho),
		a ~ dnorm(0, 100), 
		bg ~ dnorm(0, 100),
		sigma_geno ~ dcauchy(0, 2),
		sigma ~ dcauchy(0 , 2),
		Rho ~ dlkjcorr(2)
	),
	data = colsha_h3, warmup = 1000, iter = 5000, chains = 2
)

precis(colsha_stan2, depth = 2)
plot(colsha_stan2)
plot(precis(colsha_stan2, depth = 2))


height_out <- precis(colsha_stan2, depth = 2)
height_out <- height_out@output
head(height_out)
str(height_out)
height_out$name <- rownames(height_out)
# to do is figure out ggplot
colnames(height_out)[3:4] <- paste(c("low95", "high95"))
head(height_out) 
str(height_out)
height_out$name <- sub("(*)+(\\])", "\\1", height_out$name) 
height_out$name <- sub("(*)+(\\[)", "\\1", height_out$name) 
height_out$name <- sub("(*)+(,)(*)+", "\\1\\3", height_out$name) 
head(height_out)
tail(height_out, 20)


setwd("~/git.repos/NAM_CAM_CHAMBER/output/")
height_plot <- ggplot(data = height_out)
height_plot <- height_plot + geom_pointrange(mapping=aes(x = name, y = Mean, ymin=low95, ymax=high95)) + coord_flip()
height_plot
ggsave("height_hdi.pdf", height_plot, width = 10, height = 30)


post <- extract.samples(colsha_stan2)
plot(post)
head(post)

pairs(post)