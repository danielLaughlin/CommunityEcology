##################################################################################################
###  R Code and data for                                                                       ###
###  "Species Pools, Traits, Filters, and Assemply Rules: A framework For Community Ecology"   ###
###  by Paul Keddy and Daniel Laughlin                                                         ###
###                                                                                            ###
###  Code written by Daniel Laughlin, University of Wyoming, Botany Dept                       ###
##################################################################################################

# Set your working directory

## load required packages
library(vegan)
library(FD)
library(labdsv)
library(lme4)
library(lmerTest)
library(ade4)
library(piecewiseSEM)
library(MuMIn)
library(dplyr)
library(fields)
library(visreg)
library(mvtnorm)
library(mclust)
library(sads)

### Load datasets ###
load("communityEcologyData.RData")

### Notes on dataframes
# comm is the community dataset (plot by species cover)
# env is the environmental data associated with the comm dataset
# env_p is a vector of environmental data for use in the Traitspace model
# KH_data is the kettle hole wetland dataset for use in the Traitspace model
# traits is a 44 species by 11 traits dataset (average values for each species)

### Notes on important riable names in datasets
# DAYSUB = flooding duration, "days submerged per year"
# poros = aerynchyma root trait, "porosity"
# height = vegatative height



################ Begin Chapter 5: Trait-environment interactions  ######################################

### Community-weighted mean trait linear regression models ###

## compute community weighted mean trait values
ra <- as.matrix(vegan::decostand(comm, method = "total", MARGIN = 1))  # compute relative abundances for the comm dataset
cwm <- FD::functcomp(traits, ra) # compute community weighted mean traits using functcomp
cwm$QuadID <- rownames(cwm) # add rownames to cwm dataframe
cwm <- merge(cwm, env, by = "QuadID") ## merge cwm and env data

### Fit models: CWM traits ~ f(flooding) and summarize results
m1 <- lm(poros ~ DAYSUB, data = cwm)
summary(m1)

m2 <- lm(height ~ DAYSUB, data = cwm)
summary(m2)

### End of Community-weighted mean trait linear regression models ###



### Fourth Corner Analysis: a more robust test of CWM regressions ###

### create data files for fourthcorner function
tabR <- data.frame(DAYSUB = env[, "DAYSUB"]) # create environmental variable dataframe for fourthcorner function
tabL <- comm                                 # rename community dataframe for fourthcorner function
tabQ <- data.frame(height = traits[, c("height")], poros = traits[, c("poros")]) # create trait dataframe for fourthcorner function

### fit fourthcorner model and summarize results
res <- ade4::fourthcorner(tabR, tabL, tabQ)
summary(res)

### End of Fourth Corner Analysis: a more robust test of CWM regressions ###



### Generalized Linear Mixed Effects Model of Trait x Environment Interaction ###

### convert community data to long format for GLM and merge datasets
comm.long <- as.matrix(labdsv::dematrify(comm, thresh = -0.1)) # use dematrify function
colnames(comm.long) <- c("QuadID", "species", "abundance") # rename columns
traits$species <- rownames(traits) # create species column based on rownames
glmm.data <- merge(comm.long, traits, by = c("species")) # merge long form data with trait data
glmm.data <- merge(glmm.data, env, by = c("QuadID")) # merge long form data with environmental data
glmm.data$abun <- as.numeric(as.character(glmm.data$abundance)) # make abundance vector numeric
glmm.data$pa <- glmm.data$abun  # create presence/absence vector based on abundance vector, then...
glmm.data$pa[glmm.data$pa > 0] <- 1   # create presence/absence variable where 0=absent, 1=present

### scale data to unit variance to ensure GLMM convergence
glmm.data$poros <- scale(glmm.data$poros)
glmm.data$DAYSUB <- scale(glmm.data$DAYSUB)
glmm.data$height <- scale(glmm.data$height)

### Fit GLMM Environment Only Model: Binomial Presence/Absence ~ Environment + (Env|Species) + (1|Quadrat)
mod_env <- glmer(pa ~ DAYSUB + (DAYSUB|species) + (1|QuadID), data=glmm.data,
            family = binomial, control = glmerControl(optimizer = "bobyqa"))

### Fit GLMM TraitXEnvironment Model: Binomial Presence/Absence ~ Trait X Environment + (Env|Species) + (1|Quadrat)
m3 <- glmer(pa ~ poros*DAYSUB + (DAYSUB|species) + (1|QuadID), data=glmm.data,
            family = binomial, control = glmerControl(optimizer = "bobyqa"))

summary(m3) # this shows that the TraitXEnvironment interaction term is significant for aerenchyma
anova(mod_env,m3) # this shows that adding the traitXenv interaction significantly improves the model
performance::r2(m3) # this shows that the fixed effects explain 32% of the 'variation'

### Fit GLMM TraitXEnvironment Model: Binomial Presence/Absence ~ Trait X Environment + (Env|Species) + (1|Quadrat)
m4 <- glmer(pa ~ height*DAYSUB + (DAYSUB|species) + (1|QuadID), data=glmm.data,
            family = binomial, control = glmerControl(optimizer = "bobyqa"))

summary(m4) # # this shows that the TraitXEnvironment interaction term is not significant for height
anova(mod_env, m4) # this shows that adding the traitXenv interaction did not improve the model
performance::r2(m4) # this shows that the fixed effects only explain 9% of the 'variation'

### Plot random species effects ###
### If interested, you can plot the random effects for each species to see how
### the 'species effect' is accounted for in the model
flood <- seq(min(glmm.data$DAYSUB), max(glmm.data$DAYSUB), length=100) # create flood variable as a sequence
# create a matrix of predicted probabilities for each species for each value of the flood variable
curves <- matrix(0,100,44,byrow=TRUE)
for (i in 1:44){
  curves[,i] <- 1/(1+exp(-(ranef(m3)$species[i,1] + ranef(m3)$species[i,2]*flood )))
}
colnames(curves) <- rownames(traits)
curves <- as.data.frame(curves)
curves$flood <- flood
# plot probability curves for each species
plot(curves$flood, curves$AGRcap, type="l", col="white", xlab="Flooding duration", ylab="Modeled probability")
for(i in 1:44){
lines(curves$flood, curves[,i], col=1, lwd=2)
}
### End of plotting random effects ###

### End of Generalized Linear Mixed Effects Model of Trait x Environment Interaction ###



### Figure 5.4 illustrating comparison of CWM regression, fourth corner, and GLMM analysis ###

par(mfrow=c(2,2), mar=c(5,5,4,4)) # create 2x2 panel

plot(cwm$poros ~ cwm$DAYSUB, pch=19, col="grey50", ylim=c(5,60), cex.lab=1.4,
     xlab="Flooding duration (days)", ylab="Root aerenchyma (%)",
     main="a) Community-weighted mean regression
     (supported by 4th corner analysis)")
xvec<-seq(0,250,1)
legend(-5,63,legend=expression(paste("CWM r=0.77, ",italic(p), " < 0.001")),bty="n")
legend(-5,58,legend=expression(paste("4th corner r=0.67, ",italic(p), " = 0.022")),bty="n")
#mtext("A", side=3, adj=0, line=0, cex=1.1, font=1)
predframe <- data.frame(predict(m1,data.frame(DAYSUB=xvec),type="response",re.form=NA,se.fit=TRUE))
polygon(c(xvec,rev(xvec)),c(predframe$fit-2*predframe$se.fit,rev(predframe$fit+2*predframe$se.fit)),
        border = FALSE,col="grey")
lines(xvec,predict(m1,data.frame(DAYSUB=xvec)),lwd=2)
points(cwm$DAYSUB,cwm$poros,pch=19, col="grey50")

plot(cwm$height ~ cwm$DAYSUB, pch=19, col="grey60",ylim=c(0,35), cex.lab=1.4,
     xlab="Flooding duration (days)", ylab="Height (cm)",
     main="b) Community-weighted mean regression
     (unsupported by 4th corner analysis)")
legend(-5,37,legend=expression(paste("CWM r=0.33, ",italic(p), " < 0.001")),bty="n")
legend(-5,34,legend=expression(paste("4th corner r=0.26, ",italic(p), " = 0.421")),bty="n")
#mtext("A", side=3, adj=0, line=0, cex=1.1, font=1)
predframe <- data.frame(predict(m2,data.frame(DAYSUB=xvec),type="response",re.form=NA,se.fit=TRUE))
polygon(c(xvec,rev(xvec)),c(predframe$fit-2*predframe$se.fit,rev(predframe$fit+2*predframe$se.fit)),
        border = FALSE,col="grey")
lines(xvec,predict(m2,data.frame(DAYSUB=xvec)),lwd=2)
points(cwm$DAYSUB,cwm$height,pch=19, col="grey50")

par(mar=c(2,2,2,2))
grid.l<- list( seq(-3,3,length=100), seq(-3,3,length=100)) 
xy <- make.surface.grid(grid.l)
xy.prob.pred <- as.surface(xy, 1/(1 + exp(-( -5.6 + 1.49*xy[,1] - 1.51*xy[,2] + 2.12*xy[,1]*xy[,2]))))
persp(xy.prob.pred$x, xy.prob.pred$y, xy.prob.pred$z, theta = 30, phi = 20,
      xlab="\nRoot aerenchyma",ylab="\nFlooding duration",zlab="\n\nProbability
      of occurrence", ticktype="simple",
      main="c) Strong Trait-by-Environment interaction
      TxE P < 0.001, Marginal R-square = 0.33", cex.lab=1.4, lwd=0.5, border="grey40", zlim=c(0,1))

xy.prob.pred <- as.surface(xy, 1/(1 + exp(-( -5.7 + 0.51*xy[,1] - 1.53*xy[,2] -0.13*xy[,1]*xy[,2]))))
persp(xy.prob.pred$x, xy.prob.pred$y, xy.prob.pred$z, theta = 30, phi = 20,
      xlab="\nHeight",ylab="\nFlooding duration",zlab="\n\nProbability
      of occurrrence", ticktype="simple",
      main="d) No Trait-by-Environment interaction
      TxE P = 0.83, Marginal R-square = 0.10", cex.lab=1.4, lwd=0.5, border="grey40", zlim=c(0,1))

### End of Figure 5.4 illustrating comparison of CWM regression, fourth corner, and GLMM analysis ###

################ End Chapter 5: Trait-environment interactions  ######################################





################ Begin Chapter 7: Models  ################################################

### Simple CATS example in Figure 7.1 ###

avgLeafN <- 2  # set average trait value of the community
sppLeafN <- c(0.5, 1.5, 2.5) # provide average trait values of each species
output <- FD::maxent(constr=avgLeafN, states=sppLeafN) # run the maxent model
output$prob # read the output probabilities for maximum entropy solution
output$entropy # read the entropy value for the relative abundance distribution

#compute entropy for five possible solutions (convert zeroes to 0.0001)
solutionLine <- c(0, 0.25, 0.5, 0.75, 1)
entropy <- - c( (0.0001*log(0.0001))+(0.5*log(0.5))+(0.5*log(0.5)),
                (0.0581525*log(0.0581525))+(0.3837963*log(0.3837963))+(0.5581012*log(0.5581012)),
                (0.1162050*log(0.1162050))+(0.2675926*log(0.2675926))+(0.6162024*log(0.6162024)),
                (0.1831025*log(0.1831025))+(0.1338463*log(0.1338463))+(0.6831012*log(0.6831012)),
                (0.25*log(0.25))+(0.0001*log(0.0001))+(0.75*log(0.75)) )

### Plot Figure 7.1b
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(solutionLine, entropy, ylab="Entropy (Evenness)", ylim=c(0.5,1.15), xlim=c(0,1), pch=19,
     xlab="Solution set: the line in 3-D species space", lwd=5, xaxt="n", col="darkgrey")
lines(solutionLine, entropy, lty=3, col="darkgrey", lwd=3)
abline(h=1.098612, lty=2, lwd=2)
text(0.5,1.13, "maximum possible entropy for 3 species = log(3) = 1.098612", cex=0.8)
par(xpd=TRUE)
arrows(0,0.4,1,0.4, code=3, lwd=5, col=1)
par(xpd=FALSE)

### End of Simple CATS example in Figure 7.1 ###



### Demonstrating use of minEntropy function and the problem of local minimums ###
source("minEntropy.R") # this function adapts the selectSpecies.R function from Laughlin et al. (2018) Methods in E&E

### Use minEntropy function to find minimum entropy solution in Figure 7.1
# Note that we need to provide a starting value so that the descent algorithm finds the
# global rather than local solution because entropy is a convex downward function
avgLeafN <- 2  # set average trait value of the community
sppLeafN <- c(0.5, 1.5, 2.5) # provide average trait values of each species

### minimum entropy with uninformative starting solution, p={1/3,1/3,1/3}
minEntropy(t2c=as.matrix(sppLeafN), constraints=avgLeafN, t2d=as.matrix(sppLeafN),
           obj='H', pars=c(1/3,1/3,1/3))$prob

### minimum entropy with a more informative starting solution, p={0.4,0.2,0.4}
### so that global minimum is found
minEntropy(t2c=as.matrix(sppLeafN), constraints=avgLeafN, t2d=as.matrix(sppLeafN),
           obj='H', pars=c(0.4,0.2,0.4))$prob

### minimum entropy with a starting solution near max entropy, p={c0.12,0.26,0.62}
### so that global minimum is found, note this is the same global minimum
minEntropy(t2c=as.matrix(sppLeafN), constraints=avgLeafN, t2d=as.matrix(sppLeafN),
           obj='H', pars=c(0.12,0.26,0.62))$prob

### End of Demonstrating use of minEntropy function and the problem of local minimums ###



### Simulating the traitspace model ###

# load source functions
source("traitspaceFunctions.R")

### Data simulation to implement the Traitspace model
### Generate env gradient values between -3 and +3
env <- seq(-3,3,0.1)
Ntr <- length(env)
env_pred <- data.frame(env)
### Randomly generate traits to satisfy a positive trait-env relationships
t1.pos <- rnorm(Ntr,env*1.5,2)
### Randomly generate 5 species with traits within the range of the trait~env relationship
# traits range from and average of -5 to +5 across the species
# set trait variance as 1 for all species
N <- 1000
var <- 1
sp1 <- rnorm(N, -5, var)
sp2 <- rnorm(N, -2.5, var)
sp3 <- rnorm(N, 0, var)
sp4 <- rnorm(N, 2.5, var)
sp5 <- rnorm(N, 5, var)
sim.trait.data <- data.frame(
  trait = c(sp1,sp2,sp3,sp4,sp5),
  species = c(rep("sp1", N), rep("sp2", N), rep("sp3", N), rep("sp4", N), rep("sp5", N))
)
### END of Data Simulation ###

### Model calibration
# Fit a linear model on trait to quantify the environmental filter
lm1 <- lm(t1.pos ~ env)

# Fit probability density functions of traits for each of the species.
# Note: warnings in the pdf fitting are common and generally not fatal but should be examined
pdf_species_trait <- lapply(levels(sim.trait.data$species), function(x){
  Mclust(na.omit(sim.trait.data[sim.trait.data$species == x, c("trait")]),warn = TRUE)
})
names(pdf_species_trait) <- levels(sim.trait.data$species)

### Model inference stage: apply the traitspace function using the lm1 and pdf_species_trait using 500 simulated trait values
simresult <- traitspace(data = sim.trait.data, multi.mod = lm1, env_p = env_pred, PT.Sk =  pdf_species_trait, N = 500, avg.out = T, parallel = F)

### Plot Figure 7.4 Traitspace simulation ###
par(mfrow=c(1,3))
cols=c("grey10","grey30","grey50","grey70","grey90")

# T-E
plot(env,t1.pos, pch=16, xlab="Environmental gradient", ylab="Trait value", font=2, cex=1, font.lab=2, cex.lab=1.3, ylim=c(-6,6),
     main="a) Trait-environment relationship", cex.main=1.5)
pred1=predict.lm(lm1,env_pred,se.fit=TRUE,interval="prediction",level=0.95)
lines(env,pred1$fit[,1],lwd=3,col=1, lty=1)
lines(env,pred1$fit[,2],lwd=3,col=1, lty=2)
lines(env,pred1$fit[,3],lwd=3,col=1, lty=2)

# species traits
tseq <- seq(-6, 6, 0.1)
plot(dnorm(tseq,-5, 1), tseq, type="l", ylim=c(-6,6), xlab="Probability density", ylab="Trait value", font.lab=2, cex.lab=1.3, col=cols[1],
     lwd=3, main="b) Species trait distributions", cex.main=1.5, xlim=c(0.02,0.4))
lines(dnorm(tseq,-2.5, 1), tseq, col=cols[2], lwd=3)
lines(dnorm(tseq,0, 1), tseq, col=cols[3], lwd=3)
lines(dnorm(tseq,2.5, 1), tseq, col=cols[4], lwd=3)
lines(dnorm(tseq,5, 1), tseq, col=cols[5], lwd=3)

# preds
plot(env, xlim=range(env), ylim=c(0,1), xlab="Environmental gradient", main="c) Traitspace predictions",
     ylab="Predicted relative abundance", cex.lab=1.4, col="white", font.lab=2, font=2, cex.main=1.5)
points(env,simresult[,1], col=cols[1], pch=16)
points(env,simresult[,2], col=cols[2], pch=16)
points(env,simresult[,3], col=cols[3], pch=16)
points(env,simresult[,4], col=cols[4], pch=16)
points(env,simresult[,5], col=cols[5], pch=16)
legend("top",c("Species A","Species B","Species C", "Species D", "Species E"),
       col=cols, lty=c(1,1,1,1,1), lwd=c(4,4,4,4,4), cex=1.3, bty="n")

### End of Figure 7.4 ###
### End of Simulating the traitspace model ###





### CATS model of the kettlehole wetland ###

## compute community weighted mean trait values (Note: already computed in Analysis 1 above)
ra <- as.matrix(vegan::decostand(comm, method = "total", MARGIN = 1)) # compute relative abundances
cwm <- FD::functcomp(traits[,1:10], ra) # compute CWMs
cwm$QuadID <- rownames(cwm) # add rownames

### set up data set to test aereychyma
constraints <- data.frame(poros = cwm$poros) # the constraints are CWM traits
rownames(constraints) <- cwm$QuadID # add rownames to constraints dataframe
states <- data.frame(poros = traits$poros) # the states are the species trait values
rownames(states) <- rownames(traits) # add rownames to states dataframe
observed <- decostand(comm, method ="total") # create observed species abundance matrix for model testing

# fit CATS model using aerenchyma and test it
cats.m1 <- maxent(constr = constraints, states = states)
# this shows that aerenchyma is associated with species abundances
test.m1 <- maxent.test(cats.m1, obs = observed)

### set up data set to test height
constraints2 <- data.frame(height = cwm$height)
rownames(constraints2) <- cwm$QuadID
states2 <- data.frame(height = traits$height)
rownames(states2) <- rownames(traits)

# fit CATS model using height and test it
cats.m2 <- maxent(constr = constraints2, states = states2)
# this shows that height is not associated with species abundances
test.m2 <- maxent.test(cats.m2, obs = observed)

### End of CATS model of the kettlehole wetland ###





### Traitspace model of the kettlehole wetland ###

# load source functions
source("traitspaceFunctions.R")

### Test the importance of height using Traitspace

# Fit environmental filter: Individual-level traits ~ flooding
poros.lm <- lm(cbind(log(poros)) ~ sub, data = KH_data, weights = avgbio)
height.lm <- lm(cbind(log(height)) ~ sub, data = KH_data, weights = avgbio)

## Fit species trait distributions
# Note that warnings are common - refer to mclust library guide
# pdfs of aerenchyma
pdf_species_poros <- lapply(levels(KH_data$species), function(x){
  Mclust(log(na.omit(KH_data[KH_data$species==x, c("poros")])),warn=TRUE)
})
names(pdf_species_poros) <- levels(KH_data$species)

# pdfs of height
pdf_species_height <- lapply(levels(KH_data$species), function(x){
  Mclust(log(na.omit(KH_data[KH_data$species==x, c("height")])),warn=TRUE)
  })
names(pdf_species_height) <- levels(KH_data$species)


### Traitspace predictions

# fit Traitspace model for aerenchyma, Note: a warning message will appear because we are weighting the regression model by avgbio
traitspace.poros <- traitspace(data = KH_data, multi.mod = poros.lm, env_p = env_p, PT.Sk =  pdf_species_poros, N = 100, avg.out = T, parallel = F)

# fit Traitspace model for height, Note: a warning message will appear because we are weighting the regression model by avgbio
traitspace.height <- traitspace(data = KH_data, multi.mod = height.lm, env_p = env_p, PT.Sk =  pdf_species_height, N = 100, avg.out = T, parallel = F)

# load observed abundances and create dataframe for model testing
obscomm <- comm[row.names(traitspace.height), colnames(traitspace.height)] # select correct rows and columns from community dataset
obscomm <- sweep(as.matrix(obscomm), 1, rowSums(obscomm), "/") # convert to relative abundances

# test aerenchyma model significance using a permutation test
# A lower value of SES indicates a better fit.
# this test shows that aerenchyma is significantly associated with species abundances 
nullModel(obs = as.matrix(obscomm), result = traitspace.poros, permutations = 100) 

# test height model significance using a permutation test
# A lower value of SES indicates a better fit.
# this test shows that height is not significantly associated with species abundances 
nullModel(obs = as.matrix(obscomm), result = traitspace.height, permutations = 100) 

### End of Traitspace model of the kettlehole wetland ###

################ End Chapter 7: Models  ################################################








################ Chapter 8  ################################################

### CATS on Functional Groups ###

# Clustering species into functional groups by their continuous traits
# In this simple example, we create 6 functional groups based on variation in aerenchyma

# compute euclidean distance matrix with aerenchyma ("poros")
d <- vegdist(scale(traits[,c("poros")]), method="euclidean")
#d <- vegdist(scale(traits[,5]), method="euclidean")
caver <- hclust(d, method = "aver")
par(mfrow=c(1,1))
plot(caver, hang=-1)
rect.hclust(caver, 6)
grp <- cutree(caver, 6)
boxplot(poros ~ grp, ylab="Aerenchyma", xlab="Functional group", data = traits)
# add these functional group classifications to the traits dataframe
traits$fg <- grp

# make average trait dataframe for functional groups
fg_traits <- traits %>% group_by(fg) %>% summarize(sla = mean(sla), height = mean(height), poros = mean(poros))

# make new observed community matrix by functional group by summing cover within each functional group category
comm_fg <- comm
colnames(comm_fg) <- traits$fg
comm_fg_long <- dematrify(comm_fg)
comm_fg_long_summed <- comm_fg_long %>% group_by(sample,species) %>% summarize(abundance = sum(abundance))
comm_fg_mat <- matrify(data.frame(comm_fg_long_summed))

## compute community weighted mean trait values (Note: already computed in Analysis 1 above)
ra_fg <- as.matrix(vegan::decostand(comm_fg_mat, method = "total", MARGIN = 1)) # relative abundances
cwm_fg <- FD::functcomp(fg_traits[,4], ra_fg) # CWMs
cwm_fg$QuadID <- rownames(cwm_fg) # add rownames

### set up data set to test aereychyma
constraints3 <- data.frame(poros = cwm_fg$poros)
rownames(constraints3) <- cwm_fg$QuadID
states3 <- data.frame(poros = fg_traits$poros)
rownames(states3) <- rownames(fg_traits)

# fit CATS model using functional groups/aerenchyma and test it
cats.m3 <- maxent(constr = constraints3, states = states3)
test.m3 <- maxent.test(cats.m3, obs = ra_fg)
# this model for Functional Groups explains more variation than the CATS model for species

# Compare to the CATS model for species (this code was implemented above for Chapter 7)
cats.m1 <- maxent(constr = constraints, states = states)
test.m1 <- maxent.test(cats.m1, obs = observed)



### Plot of Figure 8.2 Species VS Functional group CATS predictions ###

# convert matrices of abundances into longform to plot scatterplots
obs_long <- stack(data.frame(observed))
maxent_long <- stack(data.frame(cats.m2[[1]]))
obs_fg_long <- stack(data.frame(ra_fg))
maxent_fg_long <- stack(data.frame(cats.m3[[1]]))

par(mfrow=c(1,2),mar=c(4,5,2,2))
plot(obs_long[,1], maxent_long[,1], xlab="Observed relative abundances", ylab="Maximum entropy solutions",
     main="(a) Species predictions", pch=19, col="grey50")
abline(0,1)
text(0.8,0.87,"1:1")
legend("topleft", expression(paste(R^2 ,"= 0.29")), bty="n")
plot(obs_fg_long[,1], maxent_fg_long[,1], xlab="Observed relative abundances", ylab="Maximum entropy solutions",
     main="(b) Functional group predictions", pch=19, col="grey50")
abline(0,1)
text(0.8,0.87,"1:1")
legend("topleft", expression(paste(R^2 ,"= 0.45")), bty="n")

### End of Plot of Figure 8.2 Species VS Functional group CATS predictions ###

### End of CATS on Functional Groups ###




### Minimizing entropy ###

source("minEntropy.R")
# this function adapted the selectSpecies.R function from Laughlin et al. (2018) Methods in E&E
# the only difference is that entropy is minimized rather than maximized by removing a negative sign from the function

## compute community weighted mean trait values (Note: already computed in Chapter 5 above)
ra <- as.matrix(vegan::decostand(comm, method = "total", MARGIN = 1)) # compute relative abundances
cwm <- FD::functcomp(traits[,1:10], ra) # compute CWMs
cwm$QuadID <- rownames(cwm) # add rownames

### set up data set to test aereychyma using minimum entropy
constraints <- as.matrix(cwm$poros)
rownames(constraints) <- cwm$QuadID
states <- as.matrix(traits$poros)
rownames(states) <- rownames(traits)
observed <- decostand(comm, method ="total")

### fit Minimum Entropy model
minEnt <- data.frame(matrix(0,nrow=187,ncol=44))
colnames(minEnt) <- rownames(states)
noSolution <- c() # convergence check: 0=converged
# Note the following code will take a few minutes, the model is being fit to all 187 rows of data
for(i in 1:length(constraints)){
  minEnt[i,] <- minEntropy(t2c=states, constraints=constraints[i], t2d=states, obj='H')$prob
  noSolution[i] <- minEntropy(t2c=states, constraints=constraints[i], t2d=states, obj='H')$convergence # convergence check: 0=converged
}


### Find maximum entropy solution using CATS (same code as used above in Chapter 7)
cats.maxEnt <- maxent(constr = constraints, states = states)

### use the rad fucntion from the sads package to create ranked abundances for observed communities, MaxEnt results, and MinEnt results
comm_sum <- apply(observed, 2, sum)
rad_comm_sum <- rad(as.numeric(comm_sum))
maxent_sum <- apply(cats.maxEnt[[1]], 2, sum)
rad_maxent_sum <- rad(as.numeric(maxent_sum))
minent_sum <- apply(minEnt, 2, sum)
rad_minent_sum <- rad(as.numeric(minent_sum))

# convert matrices of abundances into longform to plot scatterplots
obs_long <- stack(data.frame(observed))
maxent_long <- stack(data.frame(cats.maxEnt[[1]]))
minent_long <- stack(data.frame(minEnt))
cor.test(obs_long[,1],maxent_long[,1])

### Plot of Figure 8.4 scatterplots and RADs ###
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE),
       widths=c(1,1), heights=c(2,2.5))

plot(obs_long[,1], maxent_long[,1], xlab="Observed relative abundances", ylab="Maximum entropy solutions",
     main="(a) Maximum entropy solutions", pch=19, col="grey50")
abline(0,1)
legend("topleft", expression(paste(R^2 ,"= 0.29")), bty="n")
plot(obs_long[,1], minent_long[,1], xlab="Observed relative abundances", ylab="Minimum entropy solutions",
     main="(b) Minimum entropy solutions", pch=19, col="grey50")
abline(0,1)
legend(0.6,0.2, expression(paste(R^2 ,"= 0.08")), bty="n")

plot(rad_comm_sum, pch=19, ylim=c(0.001, 100), main="(c) Species abundance distributions (summed across all plots)")
points(rad_maxent_sum, pch=19, col="grey70")
points(rad_minent_sum, col=1)
legend("bottomleft", legend=c("Observed relative abundances","Maximum entropy solutions","Minimum entropy solutions"),
       pch=c(19,19,1), col=c(1,"grey70",1), bty="n")
### End of Plot of Figure 8.4 scatterplots and RADs ###

### End of Minimizing entropy ###

################ End of Chapter 8  ################################################

################ End of code  ################################################
