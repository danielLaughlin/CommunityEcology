
#############################################
###  Code and data for Community Ecology  ###
###  Daniel Laughlin                      ###
#############################################

# Set working directory
setwd("~/OneDrive - University of Wyoming/Data/Community Ecology/Chapter6/Code for Repo")

## load required packages
library(vegan)
library(FD)
library(labdsv)
library(lme4)
library(lmerTest)
library(ade4)
library(piecewiseSEM)
library(MuMIn)
library(sjPlot)
library(fields)
library(visreg)
library(mvtnorm)
library(mclust)

### Load datasets ###
load("communityEcologyData.RData")

### Notes on variable names
# DAYSUB = flooding duration, "days submerged per year"
# poros = aerynchyma root trait, "porosity"
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

m2 <- lm(height ~ DAYSUB, data = cwm)
summary(m2)



### ANALYSIS 2: Fourth corner analysis is a more robust test of CWM regressions ###

### create data files for fourthcorner function
tabR <- data.frame(DAYSUB = env[, "DAYSUB"])
tabL <- comm
tabQ <- data.frame(height = traits[, c("height")], poros = traits[, c("poros")])

### fit fourthcorner model and summarize results
res <- ade4::fourthcorner(tabR, tabL, tabQ)
summary(res)




### ANALYSIS 3: GLMM Trait x Environment (including zeros from sparse matrix) ###

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
m3 <- glmer(pa ~ poros*DAYSUB + (DAYSUB|species) + (1|QuadID), data=glmm.data,
            family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m3)
MuMIn::r.squaredGLMM(m3)

m4 <- glmer(pa ~ height*DAYSUB + (DAYSUB|species) + (1|QuadID), data=glmm.data,
            family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m4)
MuMIn::r.squaredGLMM(m4)

### If interested, you can plot the random effects for each species to see how
### the 'specie effect' is accounted for in the model
sjPlot::plot_model(m3, type = "re", sort.est = TRUE)
flood <- seq(min(glmm.data$DAYSUB), max(glmm.data$DAYSUB), length=100)

curves <- matrix(0,100,44,byrow=TRUE)
for (i in 1:44){
  curves[,i] <- 1/(1+exp(-(ranef(m3)$species[i,1] + ranef(m3)$species[i,2]*flood )))
}
colnames(curves) <- rownames(traits)
curves <- as.data.frame(curves)
curves$flood <- flood

cols <- viridis::viridis_pal(option="D")(ncol(curves))
plot(curves$flood, curves$AGRcap, type="l", col="white", xlab="Flooding duration", ylab="Modeled probability")
for(i in 1:44){
lines(curves$flood, curves[,i], col=cols[i], lwd=2)
}

### End of plotting random effects ###



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

################ End Chapter 6: Trait-environment interactions  ######################################




################ Begin Chapter 7: Trait-based models  ################################################


### Analysis 4: Simulating the traitspace model ###

# load source functions
source("traitspace.functions.R")

###Data simulation to implement Traitspace
###Generate env gradient values
env <- seq(-3,3,0.1)
Ntr <- length(env)
env_pred = data.frame(env)
###Generate traits to satisfy trait-env relationships
t1.pos <- rnorm(Ntr,env*1.5,2)
###Generate species with traits within the range
N=1000
var=1
sp1 <- rnorm(N,-5,var)
sp2 <- rnorm(N,-2.5,var)
sp3 <- rnorm(N,0,var)
sp4 <- rnorm(N,2.5,var)
sp5 <- rnorm(N,5,var)
sim.trait.data <- data.frame(
  trait = c(sp1,sp2,sp3,sp4,sp5),
  species = c(rep("sp1",N),rep("sp2",N),rep("sp3",N),rep("sp4",N),rep("sp5",N))
)
###END of Data Simulation ###

### Model calibration
# Fit a linear model on trait
lm1 = lm(t1.pos~env)

# Fit pdfs of traits
pdf_species_trait <- lapply(levels(sim.trait.data$species), function(x){
  Mclust(na.omit(sim.trait.data[sim.trait.data$species==x, c("trait")]),warn=TRUE)
})
names(pdf_species_trait) <- levels(sim.trait.data$species)

### Model inference stage: apply the traitspace function
simresult <- traitspace(data = sim.trait.data, multi.mod = lm1, env_p = env_pred, PT.Sk =  pdf_species_trait, N = 100, avg.out = T, parallel = F)

### Plot figure of simulation results ###
par(mfrow=c(1,3))
cols=c("red","orange","green","blue","purple")

# T-E
plot(env,t1.pos, pch=16, xlab="Environmental gradient", ylab="Trait value", font=2, cex=1, font.lab=2, cex.lab=1.3, ylim=c(-6,6),
     main="A) Trait-environment relationship", cex.main=1.5)
pred1=predict.lm(lm1,env_pred,se.fit=TRUE,interval="prediction",level=0.95)
lines(env,pred1$fit[,1],lwd=3,col=1, lty=1)
lines(env,pred1$fit[,2],lwd=3,col=1, lty=2)
lines(env,pred1$fit[,3],lwd=3,col=1, lty=2)

# species traits
tseq <- seq(-6,6,0.1)
plot(dnorm(tseq,-5, 1), tseq, type="l", ylim=c(-6,6), xlab="Probability density", ylab="Trait value", font.lab=2, cex.lab=1.3, col=cols[1],
     lwd=3, main="B) Species trait distributions", cex.main=1.5)
lines(dnorm(tseq,-2.5, 1), tseq, col=cols[2], lwd=3)
lines(dnorm(tseq,0, 1), tseq, col=cols[3], lwd=3)
lines(dnorm(tseq,2.5, 1), tseq, col=cols[4], lwd=3)
lines(dnorm(tseq,5, 1), tseq, col=cols[5], lwd=3)

# preds
plot(env, xlim=range(env), ylim=c(0,1), xlab="Environmental gradient", main="C) Traitspace predictions",
     ylab="Predicted relative abundance", cex.lab=1.4, col="white", font.lab=2, font=2, cex.main=1.5)
points(env,simresult[,1], col=cols[1], pch=16)
points(env,simresult[,2], col=cols[2], pch=16)
points(env,simresult[,3], col=cols[3], pch=16)
points(env,simresult[,4], col=cols[4], pch=16)
points(env,simresult[,5], col=cols[5], pch=16)
legend("top",c("Species A","Species B","Species C", "Species D", "Species E"),
       col=cols, lty=c(1,1,1,1,1), lwd=c(4,4,4,4,4), cex=1.3, bty="n")

### End of Figure ###







### Analysis 5: CATS model ###

## compute community weighted mean trait values (Note: already computed in Analysis 1 above)
ra <- as.matrix(vegan::decostand(comm, method = "total", MARGIN = 1)) # relative abundances
cwm <- FD::functcomp(traits, ra) # CWMs
cwm$QuadID <- rownames(cwm) # add rownames

### set up data set to test aereychyma
constraints <- data.frame(poros = cwm$poros)
rownames(constraints) <- cwm$QuadID
states <- data.frame(poros = traits$poros)
rownames(states) <- rownames(traits)
observed <- decostand(comm, method ="total")

# fit CATS model using aerenchyma and test it
cats.m1 <- maxent(constr = constraints, states = states)
test.m1 <- maxent.test(cats.m1, obs = observed)

### set up data set to test height
constraints2 <- data.frame(height = cwm$height)
rownames(constraints2) <- cwm$QuadID

states2 <- data.frame(height = traits$height)
rownames(states2) <- rownames(traits)

# fit CATS model using height and test it
cats.m2 <- maxent(constr = constraints2, states = states2)
test.m2 <- maxent.test(cats.m2, obs = observed)






### Analysis 6: Traitspace model ###

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

# test aerenchyma model significance using permutation
nullModel(obs = as.matrix(obscomm), result = traitspace.poros, permutations = 100) # A lower value of SES indicates a better fit.

# test height model significance using permutation
nullModel(obs = as.matrix(obscomm), result = traitspace.height, permutations = 100) 



