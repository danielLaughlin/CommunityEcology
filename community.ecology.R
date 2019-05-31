
#############################################
###  Code and data for Community Ecology  ###
###  Daniel Laughlin                      ###
#############################################

setwd("~/OneDrive - University of Wyoming/Data/Community Ecology/Chapter6/Code for Repo")

## load required packages
library(vegan)
library(FD)
library(mgcv)
library(labdsv)
library(lme4)
library(lmerTest)
library(piecewiseSEM)
library(ade4)
library(fields)
library(sjPlot)
library(FD)
library(visreg)
library(MuMIn)
library(mvtnorm)
library(mclust)

### Load datasets ###
load("communityEcologyData.RData")

### Notes on variable names
# DAYSUB = flooding variable, "days submerged per year"
# poros = aerynchyma variable, "porosity"
# height = vegatative height


### ANALYSIS 1: Community-weighted mean trait regression models ###

## compute community weighted mean trait values
ra <- as.matrix(vegan::decostand(comm, method = "total", MARGIN = 1)) # relative abundances
cwm <- FD::functcomp(traits, ra) # CWMs
cwm$QuadID <- rownames(cwm) # add rownames
cwm <- merge(cwm, env, by = "QuadID") ## merge cwm and env data

### Fit models: CWM traits ~ f(Env) and summarize results
m1 <- lm(poros ~ DAYSUB, data = cwm)
summary(m1)
piecewiseSEM::rsquared(m1)

m2 <- lm(height ~ DAYSUB, data = cwm)
summary(m2)
piecewiseSEM::rsquared(m2)



### ANALYSIS 2: Fourth corner analysis as more robust test of CWM regressions ###

### create data files for fourthcorner function
tabR <- data.frame(DAYSUB = env[, "DAYSUB"])
tabL <- comm
tabQ <- data.frame(height = traits[, c("height")], poros = traits[, c("poros")])

### fit fourthcorner model and summarize results
res <- ade4::fourthcorner(tabR, tabL, tabQ)
summary(res)




### ANALYSIS 3: GLMM TxE (including zeros from sparse matrix)

### convert community data to long format and merge datasets
comm.long <- as.matrix(labdsv::dematrify(comm, thresh = -0.1))
colnames(comm.long) <- c("QuadID", "species", "abundance")
traits$species <- rownames(traits)
glmm.data <- merge(comm.long, traits, by = c("species"))
glmm.data <- merge(glmm.data, env, by = c("QuadID"))
glmm.data$abun <- as.numeric(as.character(glmm.data$abundance)) # make abundance is numeric
glmm.data$pa <- glmm.data$abun  # create presence/absence variable
glmm.data$pa[glmm.data$pa > 0] <- 1   # create presence/absence variable

### scale data to ensure GLMM convergence
glmm.data$poros <- scale(glmm.data$poros)
glmm.data$DAYSUB <- scale(glmm.data$DAYSUB)
glmm.data$height <- scale(glmm.data$height)
glmm.data$DAYSUB.squared <- (glmm.data$DAYSUB - mean(glmm.data$DAYSUB))^2

### Fit GLMMs: Binomial Presence/Absence ~ Trait X Environment + (Env|Species)
m3 <- glmer(pa ~ poros*DAYSUB + (DAYSUB + DAYSUB.squared|species) + (1|QuadID), data=glmm.data,
            family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m3)
MuMIn::r.squaredGLMM(m3)

m4 <- glmer(pa ~ height*DAYSUB + (DAYSUB + DAYSUB.squared|species) + (1|QuadID), data=glmm.data,
            family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m4)
MuMIn::r.squaredGLMM(m4)

### Plot random effects for species
sjPlot::plot_model(m3, type = "re", sort.est = TRUE)
ranef(m3)$species
#...



### Figure illustrating comparison of CWM regression, fourth corner, and GLMM analysis ###

par(mfrow=c(2,2), mar=c(5,5,4,4))

plot(cwm$poros ~ cwm$DAYSUB, pch=19, col="grey50", ylim=c(5,60), cex.lab=1.4,
     xlab="Flooding duration (days)", ylab="Root aerenchyma (%)",
     main="A) Community-weighted mean regression
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
     main="B) Community-weighted mean regression
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
      xlab="\nRoot aerynchyma",ylab="\nFlooding duration",zlab="\n\nProbability
      of occurrence", ticktype="simple",
      main="C) Strong Trait-by-Environment interaction
      TxE P < 0.001, Marginal R-square = 0.33", cex.lab=1.4, lwd=0.5, border="grey40", zlim=c(0,1))

xy.prob.pred <- as.surface(xy, 1/(1 + exp(-( -5.7 + 0.51*xy[,1] - 1.53*xy[,2] -0.13*xy[,1]*xy[,2]))))
persp(xy.prob.pred$x, xy.prob.pred$y, xy.prob.pred$z, theta = 30, phi = 20,
      xlab="\nHeight",ylab="\nFlooding duration",zlab="\n\nProbability
      of occurrrence", ticktype="simple",
      main="D) No Trait-by-Environment interaction
      TxE P = 0.83, Marginal R-square = 0.10", cex.lab=1.4, lwd=0.5, border="grey40", zlim=c(0,1))

### End Figure #######################################################################################




### Analysis 4: CATS model ###

### set up data set to test aereychyma
constraints <- data.frame(poros = cwm$poros)
rownames(constraints) <- cwm$QuadID
states <- data.frame(poros = traits$poros)
rownames(states) <- rownames(traits)
observed <- decostand(comm, method ="total")

# fit CATS model and test it
cats.m1 <- maxent(constr = constraints, states = states)
test.m1 <- maxent.test(cats.m1, obs = observed)


### set up data set to test height
constraints2 <- data.frame(height = cwm$height)
rownames(constraints2) <- cwm$QuadID

states2 <- data.frame(height = traits$height)
rownames(states2) <- rownames(traits)

# fit CATS model and test it
cats.m2 <- maxent(constr = constraints2, states = states2)
test.m2 <- maxent.test(cats.m2, obs = observed)





### Analysis 5: Traitspace model ###

# load source functions
source("traitspace.functions.R")

### Test the importance of height using Traitspace

# Fit environmental filter: Individual-level height ~ flooding
poros.lm <- lm(cbind(log(poros)) ~ sub, data = KH.data, weights = avgbio)
height.lm <- lm(cbind(log(height)) ~ sub, data = KH.data, weights = avgbio)

## Fit species trait distributions
# Note that warnings are common - refer to mclust library guide

# pdfs of aerenchyma
pdf_species_poros <- lapply(levels(KH.data$species), function(x){
  Mclust(log(na.omit(KH.data[KH.data$species==x, c("poros")])),warn=TRUE)
})
names(pdf_species_poros) <- levels(KH.data$species)

# pdfs of height
pdf_species_height <- lapply(levels(KH.data$species), function(x){
  Mclust(log(na.omit(KH.data[KH.data$species==x, c("height")])),warn=TRUE)
  })
names(pdf_species_height) <- levels(KH.data$species)


### Traitspace predictions

# fit model for aerenchyma
traitspace.poros <- traitspace(data = KH.data, multi.mod = poros.lm, env_p = env_p, PT.Sk =  pdf_species_poros, N = 100, avg.out = T, parallel = F)

# fit model for height
traitspace.height <- traitspace(data = KH.data, multi.mod = height.lm, env_p = env_p, PT.Sk =  pdf_species_height, N = 100, avg.out = T, parallel = F)

# load observed data
obscomm <- comm[row.names(traitspace.height), colnames(traitspace.height)] # select correct rows and columns from community dataset
obscomm <- sweep(as.matrix(obscomm), 1, rowSums(obscomm), "/") # convert to relative abundances

# test model significance
nullModel(obs = as.matrix(obscomm), result = traitspace.poros, permutations = 100) # A lower value of SES indicates a better fit.
nullModel(obs = as.matrix(obscomm), result = traitspace.height, permutations = 100) 



