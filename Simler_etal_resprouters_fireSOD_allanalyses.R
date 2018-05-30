## Analyses for: "Novel disturbance interactions between fire and an emerging disease impact survival and growth of resprouting trees"

###### CONTENTS #####################################

## PART 1: Examining effect of disease-related fuels on belowground mortality in fire
## PART 2: Examining effect of local measurements of fire severity on belowground mortality
## PART 3: Assessing ties between coarse woody debris and fire severity metrics
## PART 4: Examining effect of disturbance severity on post-fire resprouting
## PART 5: Examining effect of pre-fire SOD-related mortality on post-fire relative dominance.
## PART 6: Comparison models from unburned plots: Resprouter mortality and resprouting vigor

library(rstan)
library(ggplot2)
library(loo)
library(shinystan)
library(rstanarm)
library(dplyr)
library(bayesplot)

## PART 1: Examining effect of disease-related fuels on belowground mortality in fire ##########

###### 1.1) Data organization: #######
## For access to datasets, contact absimler@ucdavis.edu ####
setwd("~/Dropbox/Data Analysis/resproutingfireSOD")
Plot = read.csv("datasets/Postfire_Plots_incl2009.csv", header=T, na.strings=c("", " "))
Mortality = read.csv("datasets/TreeMortality_incl2009.csv", header=T, na.strings=c("", " "))


## For the belowground mortality analysis you will need the "Mortality" and the "Plot" datasets created in the data setup script.
#Individual variables are:
# Plot = PlotID
# TotKilled = 0,1 variable for whether tree experienced belowground mortality
# Species = Species ID
#LiveBAc = pre-fire basal area (log)
#Canker06 = prefire canker presence

# PLot-level variables are:
# HostsA = plot-level basal area of standing dead hosts retaining fine fuels
# HostsC = plot-level basal area of standing dead hosts not retaining fine fuels
# HostCWD = plot-level coarse woody debris from hosts
# PFtemp = Post-fire mean annual temperature
# BigSurPlot = Plot ID

##Assemble individual data for tanoaks from Mortality file
deadL <- subset(Mortality, Species=="LIDE" & Burned==TRUE & ForestType=="Redwood" &  Live.BA.0607 <0.30 ) ## reduce data set to trees below certain size to compare redwoods and tanoaks, and to burned redwood plots only
deadL <- deadL[c("Plot", "TotKilled", "Species", "Live.BA.0607", "Treelist", "Canker06")]
colnames(deadL)[4] <- "liveBA"
deadL$LiveBAlog <- log(deadL$liveBA + 0.01)
deadL$LiveBAc <- (deadL$LiveBAlog - mean(deadL$LiveBAlog) ) / sd( deadL$LiveBAlog)
deadL$plotindex <- group_indices(deadL, Plot)
deadL$TotKilled2[deadL$TotKilled==1] <- 0 ## Change mortality to survival, per reviewer recommendations.
deadL$TotKilled2[deadL$TotKilled==0] <- 1
deadL$TotKilled <- as.numeric(deadL$TotKilled2)
deadL$Canker06 <- as.numeric(as.character(deadL$Canker06))
deadL$Canker06 <- (deadL$Canker06 - mean(deadL$Canker06) ) / sd( deadL$Canker06)
deadL[is.na(deadL)]<-0

## Plot-level data for tanoaks:
#Subset to only plots present in tanoak data, scale and center data:
MortPlot <- subset(Plot, BigSurPlot %in% unique(deadL$Plot))
MortPlot <- MortPlot[c("BigSurPlot", "hostCWD", "hostsA", "hostsC", "PFtemp")]
MortPlot$hostCWDc <- (MortPlot$hostCWD - mean(MortPlot$hostCWD) ) / sd( MortPlot$hostCWD)
MortPlot$hostsAc<- (MortPlot$hostsA - mean(MortPlot$hostsA) ) / sd( MortPlot$hostsA)
MortPlot$hostsCc<- (MortPlot$hostsC - mean(MortPlot$hostsC) ) / sd( MortPlot$hostsC)
MortPlot$tempc <- (MortPlot$PFtemp - mean(MortPlot$PFtemp) ) / sd( MortPlot$PFtemp)


##Assemble individual data for redwoods from Mortality file
deadRW <- subset(Mortality, Species=="SESE" & Burned==TRUE & ForestType=="Redwood")
deadRW <- subset(deadRW, Live.BA.0607 <0.30 ) ## reduce data set to trees below certain size to compare redwoods and tanoaks, and to burned redwood plots only
deadRW <- deadRW[c("Plot", "TotKilled", "Species", "Live.BA.0607")]
colnames(deadRW)[4] <- "liveBA"
deadRW$LiveBAlog <- log(deadRW$liveBA + 0.01)
deadRW$LiveBAc <- (deadRW$LiveBAlog - mean(deadRW$LiveBAlog) ) / sd( deadRW$LiveBAlog)
deadRW$plotindex <- group_indices(deadRW, Plot)
deadRW[is.na(deadRW)]<-0
deadRW$TotKilled2[deadRW$TotKilled==1] <- 0
deadRW$TotKilled2[deadRW$TotKilled==0] <- 1
deadRW$TotKilled <- as.numeric(deadRW$TotKilled2)

## Assemble plot level data for redwoods
MortPlot2 <- subset(Plot, BigSurPlot %in% unique(deadRW$Plot))
MortPlot2 <- MortPlot2[c("BigSurPlot", "hostCWD", "hostsA", "hostsC", "PFtemp")]
MortPlot2$hostCWDc <- (MortPlot2$hostCWD - mean(MortPlot2$hostCWD) ) / sd( MortPlot2$hostCWD)
MortPlot2$hostsAc<- (MortPlot2$hostsA - mean(MortPlot2$hostsA) ) / sd( MortPlot2$hostsA)
MortPlot2$hostsCc<- (MortPlot2$hostsC - mean(MortPlot2$hostsC) ) / sd( MortPlot2$hostsC)
MortPlot2$tempc <- (MortPlot2$PFtemp - mean(MortPlot2$PFtemp) ) / sd( MortPlot2$PFtemp)


###### 1.2) Null model, includes only a random intercept for plot #######
##### Run null model for LIDE: ##
data <- list( TotKilled = deadL$TotKilled, Plot = deadL$plotindex, N = nrow(deadL), P = length(unique(deadL$plotindex)))

null.mortality.lide <- stan(file = "nullmortality.stan",
                            data = data, 
                            chains=4 , iter=2000 , warmup=1000, thin=1,
                            control=list(adapt_delta = 0.99, stepsize = 0.01))

#traceplot(null.mortality.lide)

#extract log liklihood for loo scores
loglik <- extract_log_lik(null.mortality.lide, parameter_name = "log_lik")
null.lide.loo <- loo(loglik, k_threshold=0.7)


##### Run null model for SESE: #
data <- list( TotKilled = deadRW$TotKilled, Plot = deadRW$plotindex, N = nrow(deadRW), P = length(MortPlot2$BigSurPlot))

null.mortality.rw <- stan(file = "nullmortality.stan",
                          data = data, 
                          chains=4 , iter=2000 , warmup=1000, thin=1,
                          control=list(adapt_delta = 0.99, 
                                       stepsize = 0.01))

#traceplot(null.mortality.rw)

#extract log liklihood for loo scores
loglik <- extract_log_lik(null.mortality.rw, parameter_name = "log_lik")
null.sese.loo <- loo(loglik, k_threshold=0.7)


###### 1.3) Run Full model, includes standing dead fuels (A and C stages) + plot level CWD + pre-fire tree size. #######

###### Run full model for LIDE 
## Simulation data for plotting model predictions in gen. quantities -- right now, all set to mean except for dbh effect, but change to plot other effects.
cwd.seq <-rep(0, 300)
dbh.seq <- rep(seq(from = min(deadL$LiveBAc), to = max(deadL$LiveBAc), length.out = 100), 3)
Aseq <- rep(0, 300)
Bseq <- rep(0, 300)

# Assemble data frame for tanoak model fit:
data <- list( TotKilled = deadL$TotKilled, Plot = deadL$plotindex, LiveBAc = deadL$LiveBAc, hostCWD = MortPlot$hostCWDc, N = nrow(deadL), P = length(unique(deadL$plotindex)), climate= MortPlot$tempc, standingA = MortPlot$hostsAc, standingB = MortPlot$hostsCc, Canker=deadL$Canker06, Aseq=Aseq, Bseq=Bseq, cwdseq = cwd.seq, dbhseq = dbh.seq, Nsims= length(cwd.seq))

#fit model for tanoak:
full.mortality.lide <- stan(file = "fullmortalityhost.stan",
                            data = data, 
                            chains=4 , iter=2000 , warmup=1000, thin=1,
                            control=list(adapt_delta = 0.99, stepsize = 0.01))


#traceplot(full.mortality.lide)

plot(full.mortality.lide, pars=c("a", "bdbh", "ba", "bb", "bcwd", "bclim", "bsymp", "sigma_plot"), outer_level=0.90)
full.mort.LIDE.summary <- round(summary(full.mortality.lide, pars=c("a", "bdbh", "ba", "bb", "bcwd", "bclim", "bsymp", "sigma_plot"), probs = c(0.05, 0.95))$summary, digits=2)


#extract logliklihood from generated quantities to calculate LOO
loglik <- extract_log_lik(full.mortality.lide, parameter_name = "log_lik")
full.lide.loo <- loo(loglik, k_threshold=0.7)


######  Run full model for SESE
cwd.seq2 <- rep(0, 300)
dbh.seq2 <- rep(seq(from = min(deadRW$LiveBAc), to = max(deadRW$LiveBAc), length.out = 100), 3)
Aseq2 <- rep(0, 300)
Bseq2 <- rep(0, 300)
Cseq2 <- rep(0, 300)
data <- list( TotKilled = deadRW$TotKilled, Plot = deadRW$plotindex, LiveBAc = deadRW$LiveBAc, hostCWD = MortPlot2$hostCWDc, N = nrow(deadRW), P = length(unique(deadRW$plotindex)), standingA = MortPlot2$hostsAc, standingB = MortPlot2$hostsCc, disstage= MortPlot2$Pram, Aseq=Aseq2, Bseq=Bseq2, cwdseq = cwd.seq2, dbhseq = dbh.seq2, Nsims= length(cwd.seq2), climate= MortPlot2$tempc)

full.mortality.sese <- stan(file = "fullmortality.stan",
                             data = data, 
                             chains=4 , iter=2000 , warmup=1000, thin=1,
                             control=list(adapt_delta = 0.99, 
                                          stepsize = 0.01))

plot(full.mortality.sese, pars=c("a", "bdbh", "ba", "bb", "bcwd", "bclim"), outer_level=0.9)
#traceplots(full.mortality.sese)
full.mort.SESE.summary <- summary(full.mortality.sese, pars=c("a", "bdbh", "ba", "bb", "bcwd", "bclim"), probs = c(0.05, 0.95))$summary

#extract log liklihood for loo scores
loglik <- extract_log_lik(full.mortality.sese, parameter_name = "log_lik")
full.sese.loo <- loo(loglik, k_threshold=0.7)

###### 1.4) Model comparison #######

#compare to the full model:

compare.lide.mort <- compare(null.lide.loo, full.lide.loo)
compare.sese.mort <- compare(null.sese.loo, full.sese.loo)

compare.lide.mort
compare.sese.mort

## PART 2: Examining effect of local measurements of fire severity on belowground mortality #######

#### 2.1) Data organization #####
setwd("~/Dropbox/Data Analysis/resproutingfireSOD")
BurnTrees = read.csv("datasets/Fireseverity_Trees.csv", header=T, na.strings=c("", " "))

## The belowground mortality analysis requires the "BurnTrees" dataset from the setup script.##
##litterc = litter burn severity, centered
##soilc = soil burn severity, centered
##duffc = intact duff level (in cm), centered
##LiveBA = pre-fire basal area of individual
##PLot = Plot ID

#subset to just LIDEs beneath our size threshold, within 5m of sampling point
BurnSevLIDE <- subset(BurnTrees, Species=="LIDE" & sampdist<5 & liveBA<0.3)
#Data setup -- scale and center
BurnSevLIDE$litterc <- scale(BurnSevLIDE$litterCBI)
BurnSevLIDE$duffc <- scale(BurnSevLIDE$intactduff)
BurnSevLIDE$soilc <- scale(BurnSevLIDE$burnfactor)
BurnSevLIDE$LiveBAc <- scale(BurnSevLIDE$LiveBAlog)
#### 2.2) Run Models ########

## Null model
null <- stan_glmer(TotKilled2 ~  (1|Plot),
                   data = BurnSevLIDE, family = binomial,
                   prior_intercept = normal(0, 5, autoscale = FALSE),
                   chains = 3, cores = 1, iter = 2000, warmup=1000 )

loonull <- loo(null, k_threshold=0.7)

## Effect of soil severity
soil <- stan_glmer(TotKilled2 ~ LiveBAc + soilc + (1|Plot),
                   data = BurnSevLIDE, family = binomial,
                   prior = normal(0, 5, autoscale = FALSE),
                   prior_intercept = normal(0, 5, autoscale = FALSE),
                   chains = 3, cores = 1, iter = 2000, warmup=1000 )

loosoil <- loo(soil, k_threshold=0.7)
summary(soil, digits=2)
plot(soil)

### Effect of duff severity
duff <- stan_glmer(TotKilled2 ~ LiveBAc  + duffc + (1|Plot),
                   data = BurnSevLIDE, family = binomial,
                   prior = normal(0, 5, autoscale = FALSE),
                   prior_intercept = normal(0, 5, autoscale = FALSE),
                   chains = 3, cores = 1, iter = 2000, warmup=1000 )
looduff <- loo(duff, k_threshold = 0.7)


### Effect of litter severity
litter <- stan_glmer(TotKilled2 ~ LiveBAc + litterc +  (1|Plot),
                     data = BurnSevLIDE, family = binomial,
                     prior = normal(0, 5, autoscale = FALSE),
                     prior_intercept = normal(0, 5, autoscale = FALSE),
                     chains = 3, cores = 1, iter = 2000, warmup=1000 )

plot(litter)
loolitter <- loo(litter, k_threshold=0.7)

duff2 <- stan_glmer(TotKilled2 ~ LiveBAc + duffc +  (1|Plot),
                     data = BurnSevLIDE, family = binomial,
                     prior = normal(0, 5, autoscale = FALSE),
                     prior_intercept = normal(0, 5, autoscale = FALSE),
                     chains = 3, cores = 1, iter = 2000, warmup=1000 )

### 2.3) Model comparison ####

compare(loosoil, looduff, loolitter, loonull) 


## PART 3: Relationship between coarse woody debris and fire severity #######

## These models requires the "SevPlot" dataset from the setup script.#
setwd("~/Dropbox/Data Analysis/resproutingfireSOD")
SevPlot = read.csv("datasets/Fireseverity_Plots.csv", header=T, na.strings=c("", " "))
### 3.1) Run models #####

cwd.duff <- stan_glmer(meanduff ~ hostCWD +(1|plot),
                       data = SevPlot, family = gaussian,
                       prior = normal(0, 5, autoscale = FALSE),
                       prior_intercept = normal(0, 1, autoscale = FALSE),
                       chains = 3, cores = 1, iter = 2000, warmup=1000, adapt_delta=0.999999)
summary(cwd.duff, digits =2)
plot((cwd.duff), par="hostCWD")

cwd.litter <- stan_glmer(litter.mean ~ hostCWD +(1|plot),
                         data = SevPlot, family = gaussian,
                         prior = normal(0, 5, autoscale = FALSE),
                         prior_intercept = normal(0, 1, autoscale = FALSE),
                         chains = 3, cores = 1, iter = 2000, warmup=1000, adapt_delta=0.999999 )
summary(cwd.litter)
plot((cwd.litter), par="hostCWD")

cwd.soil <- stan_glmer(soil.mean ~ hostCWD +(1|plot),
                       data = SevPlot, family = gaussian,
                       prior = normal(0, 5, autoscale = FALSE),
                       prior_intercept = normal(0, 1, autoscale = FALSE),
                       chains = 3, cores = 1, iter = 2000, warmup=1000, adapt_delta=0.9999 )
summary(cwd.soil)
plot((cwd.soil), par="hostCWD")


## PART 4: Examining effect of disturbance severity on post-fire resprouting ######

### 4.1) Data organization ####

## For the resprouting vigor analysis you will need the "RSTrees" and the "GrowthPlot" datasets created in the data setup script.
#Individual variables are:
# Plot = PlotID
# TotNewSTems = total number of new stems (<1cm dbh) that have recruited post-fire
# Species = Species ID
#LiveBAc = pre-fire basal area (log)
#DeadBAc = post-fire dead basal area due to fire (log)
#Pr2013 = for tanoak, post-fire infection by P.ramorum
#Canker06 = prefire canker presence

# PLot-level variables are:
# FireMortc = change in living basal area pre- and post-fire
# DeadBA06c = Pre-fire dead basal area of hosts
# tempc = Post-fire mean annual temperature

##Load the datasets that were formatted by setup file:
setwd("~/Dropbox/Data Analysis/resproutingfireSOD")
Plot = read.csv("datasets/Postfire_Plots_longterm.csv", header=T, na.strings=c("", " "))
Uniquetree = read.csv("datasets/BStrees_longterm.csv", header=T, na.strings=c("", " "))

### Data organization for tanoak:

#TANOAK Data:
#Subset individual data to tanoaks that were alive pre-fire and survived to resprout:
spLIDE = subset(Uniquetree, Species=="LIDE" & Status2006=='L' & Status2010=='L' & Sprouting=="Y" & Burned==TRUE)
spLIDE <- spLIDE[c("Plot", "TotNewStems", "Species", "Burned", "Pr2013", "Canker06", "LBA0607", "DBA10")]
spLIDE$LiveBAc <- (spLIDE$LBA0607 - mean(spLIDE$LBA0607) ) / sd( spLIDE$LBA0607)
spLIDE$DeadBAc <- (spLIDE$DBA10 - mean(spLIDE$DBA10) ) / sd( spLIDE$DBA10)
spLIDE <- transform(spLIDE,plotindex=as.numeric(factor(spLIDE$Plot)))
spLIDE$Canker06 <- as.numeric(as.character(spLIDE$Canker06))
spLIDE$Pr2013 <- as.numeric(as.character(spLIDE$Pr2013))
spLIDE[is.na(spLIDE)]<-0

# Plot-level variables:
# Subset to just plots present in tanoak data & variables needed:
GrowthPlot <- subset(Plot, BigSurPlot %in% unique(spLIDE$Plot))
GrowthPlot <- GrowthPlot[c("BigSurPlot", "FireMort", "DeadBA0607", "DeadBA10", "PFprecip5", "PFtemp5", "Burned_2008")]

#Scale, center, make variables numeric, etc:
GrowthPlot <- transform(GrowthPlot,plotindex=as.numeric(factor(GrowthPlot$BigSurPlot))) #Make plot-level index
GrowthPlot$DeadBA0607 <- as.numeric(as.character(GrowthPlot$DeadBA0607))
GrowthPlot$DeadBA06c <- (GrowthPlot$DeadBA0607 - mean(GrowthPlot$DeadBA0607) ) / sd( GrowthPlot$DeadBA0607)
GrowthPlot$tempc <- (GrowthPlot$PFtemp5- mean(GrowthPlot$PFtemp5) ) / sd( GrowthPlot$PFtemp5)
GrowthPlot$FireMortc<-  (GrowthPlot$FireMort - mean(GrowthPlot$FireMort ) ) / sd( GrowthPlot$FireMort)

# REDWOOD data:
#Subset individual data to redwoods that were alive pre-fire and survived to resprout:
spSESE = subset(Uniquetree, Species=="SESE" & Status2006=='L' & Status2010=='L' & Sprouting=="Y" & Live.BA.0607>0)

spSESE <- spSESE[c("Plot", "TotNewStems", "Species", "LBA0607", "DBA10", "Burned", "Treelist")]
spSESE <- subset(spSESE, Burned==TRUE)
spSESE$LiveBAc <- (spSESE$LBA0607 - mean(spSESE$LBA0607) ) / sd( spSESE$LBA0607)
spSESE$DeadBAc <- (spSESE$DBA10 - mean(spSESE$DBA10) ) / sd( spSESE$DBA10)
spSESE$plotindex <- group_indices(spSESE, Plot)
spSESE[is.na(spSESE)]<-0

## Subset plot-level data to only plots included in redwood individual dataset:
GrowthPlot2 <- subset(Plot, BigSurPlot %in% unique(spSESE$Plot))
GrowthPlot2 <- GrowthPlot2[c("BigSurPlot", "FireMort", "DeadBA0607", "DeadBA10", "PFprecip5", "PFtemp5", "Burned_2008")]

#Scale, center, make variables numeric, etc
GrowthPlot2$FireMortc<-  (GrowthPlot2$FireMort - mean(GrowthPlot2$FireMort ) ) / sd( GrowthPlot2$FireMort)
GrowthPlot2$DeadBA06c <- (GrowthPlot2$DeadBA0607- mean(GrowthPlot2$DeadBA0607) ) / sd( GrowthPlot2$DeadBA0607)
GrowthPlot2$tempc <- (GrowthPlot2$PFtemp- mean(GrowthPlot2$PFtemp) ) / sd( GrowthPlot2$PFtemp)

##### 1) Null model (only has plot level random effect) #####

data <- list( TotNewStems = spLIDE$TotNewStems, Plot = spLIDE$plotindex, N = nrow(spLIDE), P = length(unique(spLIDE$plotindex)))

null.RS.LIDE <- stan(file = "nullresprouting.stan",
                     data = data, 
                     chains=4 , iter=2000 , warmup=1000, thin=1,
                     control=list(adapt_delta = 0.99, stepsize = 0.01))


#traceplot(null.RS.LIDE)
#summary(null.RS.LIDE, depth=1)

loglik <- extract_log_lik(null.RS.LIDE, parameter_name = "log_lik")
null.rslide.loo <- loo(loglik, k_threshold=0.7)


### Run model for SESE
data <- list( TotNewStems = spSESE$TotNewStems, Plot = spSESE$plotindex, N = nrow(spSESE), P = length(unique(GrowthPlot2$BigSurPlot)))

null.RS.SESE <- stan(file = "nullresprouting.stan",
                     data = data, 
                     chains=4 , iter=2000 , warmup=1000, thin=1,
                     control=list(adapt_delta = 0.99, stepsize = 0.01))

loglik <- extract_log_lik(null.RS.SESE, parameter_name = "log_lik")
null.rssese.loo <- loo(loglik)

#####2) Full model (plot level random effect, individual size predictor, and plot level effects for fire mortality, disease mortality) ######

## Run model for LIDE:
dbh.seq = rep(seq(from = min(spLIDE$LiveBAc), to = max(spLIDE$LiveBAc), length.out = 100), 3)
top.seq = rep(0, 300)
sod.seq = rep(c(-0.25, 0.5, 1.1), each = 100 )
fire.seq= rep(0, 300)
spLIDE$Pr <- as.numeric(as.character(spLIDE$Pr2013))
Pr.seq = rep(mean(spLIDE$Pr), 300)
Symp.seq = rep(mean(spLIDE$Canker06), 300)

data <- list(TotNewStems = spLIDE$TotNewStems, Plot = spLIDE$plotindex, N = nrow(spLIDE), P = length(unique(GrowthPlot$BigSurPlot)), LiveBAc = spLIDE$LiveBAc, SODdead=GrowthPlot$DeadBA06c, Firedead = GrowthPlot$FireMortc, Pr=spLIDE$Pr, symptom06=spLIDE$Canker06, interact = (GrowthPlot$FireMortc)*(GrowthPlot$DeadBA06c), TopkillBA = spLIDE$DeadBAc, climate=GrowthPlot$tempc, dbhsim = dbh.seq, topkillsim = top.seq , DeadSODsim = sod.seq, DeadFiresim = fire.seq, Nsims = length(dbh.seq), Prsim = Pr.seq, Sympsim = Symp.seq)

full.RS.LIDE <- stan(file = "fullresproutinghost.stan",
                     data = data, 
                     chains=4 , iter=2000 , warmup=1000, thin=1,
                     control=list(adapt_delta = 0.99, stepsize = 0.01))


#traceplot(full.RS.LIDE)
plot(full.RS.LIDE, pars=c("a", "bdbh", "btop", "bSOD", "bfire", "bclim", "bsymp", "bPr" ),  outer_level=0.90)
full.RS.LIDE.summary <- round(summary(full.RS.LIDE, pars=c("a", "bdbh", "btop", "bSOD", "bfire", "bclim", "bsymp", "bPr", "sigma_plot" ), probs = c(0.05, 0.95))$summary, digits=2)

## Log likelihood for loo score
loglik <- extract_log_lik(full.RS.LIDE, parameter_name = "log_lik")
full.rslide.loo <- loo(loglik, k_threshold=0.7)


###### SESE model run:

# Simulated data for plotting model predictions -- right now, all set to mean except tree size. Change these to plot effects.
dbh.seq = rep(seq(from = min(spSESE$LiveBAc), to = max(spSESE$LiveBAc), length.out = 100), 3)
top.seq = rep(0, 300)
sod.seq= rep(0, 300)
fire.seq= rep(0, 300 )

data <- list( TotNewStems = spSESE$TotNewStems, Plot = spSESE$plotindex, N = nrow(spSESE), P = length(unique(spSESE$plotindex)), LiveBAc = spSESE$LiveBAc, SODdead=GrowthPlot2$DeadBA06c, Firedead = GrowthPlot2$FireMortc, climate=GrowthPlot2$tempc, TopkillBA = spSESE$DeadBAc,  dbhsim = dbh.seq, topkillsim = top.seq , DeadSODsim = sod.seq, DeadFiresim = fire.seq, Nsims = length(dbh.seq) )

full.RS.SESE <- stan(file = "fullresprouting.stan",
                     data = data, 
                     chains=4 , iter=2000 , warmup=1000, thin=1,
                     control=list(adapt_delta = 0.99, stepsize = 0.01))

#traceplot(full.RS.SESE)
plot(full.RS.SESE, pars=c("a", "bdbh", "btop", "bSOD", "bfire", "bclim"),  outer_level=0.90)
full.RS.SESE.summary <- round(summary(full.RS.SESE, pars=c("a", "bdbh", "btop", "bSOD", "bfire", "bclim", "sigma_plot" ), probs = c(0.05, 0.95))$summary, digits=2)

# Likelihood for calculating loo
loglik <- extract_log_lik(full.RS.SESE, parameter_name = "log_lik")
full.rssese.loo <- loo(loglik)


### 3) Compare models #######
compare.sese.rs <- compare(null.rssese.loo, full.rssese.loo)
compare.lide.rs <- compare(null.rslide.loo, full.rslide.loo)


## PART 5: Examining effect of pre-fire SOD-related mortality on post-fire relative dominance. #######

setwd("~/Dropbox/Data Analysis/resproutingfireSOD")
##Load the datasets that were formatted by Chapter 1 data setup for SPROUTING file.
Plot = read.csv("datasets/Postfire_Plots_longterm.csv", header=T, na.strings=c("", " "))

ImpPlot = subset(Plot, Burned_2008==TRUE)
ImpPlot$LIDE.BA.LIVE.0607 <- as.numeric(as.character(ImpPlot$LIDE.BA.LIVE.0607))
ImpPlot$UMCA.BA.LIVE.0607 <- as.numeric(as.character(ImpPlot$UMCA.BA.LIVE.0607))
ImpPlot$SESE.BA.LIVE.0607 <- as.numeric(as.character(ImpPlot$SESE.BA.LIVE.0607))
ImpPlot$LIDE.BA.LIVE.13 <- as.numeric(as.character(ImpPlot$LIDE.BA.LIVE.13))
ImpPlot$UMCA.BA.LIVE.13 <- as.numeric(as.character(ImpPlot$UMCA.BA.LIVE.13))
ImpPlot$SESE.BA.LIVE.13 <- as.numeric(as.character(ImpPlot$SESE.BA.LIVE.13))
ImpPlot[is.na(ImpPlot)] <- 0

#Total basal area of dominant tree species
ImpPlot$TotWoody06 <- ImpPlot$LIDE.BA.LIVE.0607 + ImpPlot$SESE.BA.LIVE.0607 + ImpPlot$UMCA.BA.LIVE.0607
ImpPlot$TotWoody13 <- ImpPlot$LIDE.BA.LIVE.13 + ImpPlot$SESE.BA.LIVE.13 + ImpPlot$UMCA.BA.LIVE.13
ImpPlot$ImpLIDE06 <- ImpPlot$LIDE.BA.LIVE.0607/ImpPlot$TotWoody06 #Relative importance
ImpPlot$ImpLIDE13 <- ImpPlot$LIDE.BA.LIVE.13/ImpPlot$TotWoody13 #Relative importance

ImpPlot$ImpSESE06 <- ImpPlot$SESE.BA.LIVE.0607/ImpPlot$TotWoody06 #Relative importance
ImpPlot$ImpSESE13 <- ImpPlot$SESE.BA.LIVE.13/ImpPlot$TotWoody13 #Relative importance

SESEimp <- stan_glm(ImpSESE13 ~ ImpSESE06 + DeadBA06,
                    data = ImpPlot, family = gaussian,
                    prior = normal(0, 5, autoscale = FALSE),
                    prior_intercept = normal(0, 1, autoscale = FALSE),
                    chains = 3, cores = 1, iter = 2000, warmup=1000, adapt_delta=0.9999)
plot(SESEimp)
SESEimpsummary<- summary(SESEimp, digits=2, probs = c(0.05, 0.95))

LIDEimp <- stan_glm(ImpLIDE13 ~ ImpLIDE06 + DeadBA06,
                    data = ImpPlot, family = gaussian,
                    prior = normal(0, 5, autoscale = FALSE),
                    prior_intercept = normal(0, 1, autoscale = FALSE),
                    chains = 3, cores = 1, iter = 2000, warmup=1000, adapt_delta=0.99)
LIDEimpsummary <- summary(LIDEimp, digits=2, probs=c(0.05, 0.95))
plot((LIDEimp))


## PART 6: Comparison models from unburned plots: Resprouter mortality and resprouting vigor #######

### 6.1) Belowground mortality unburned comparison model: #######

setwd("~/Dropbox/Data Analysis/resproutingfireSOD")
Plot = read.csv("datasets/Postfire_Plots_incl2009.csv", header=T, na.strings=c("", " "))
Mortality = read.csv("datasets/TreeMortality_incl2009.csv", header=T, na.strings=c("", " "))
m=match(Mortality$Plot, Plot$BigSurPlot) #match back to tree IDs
Mortality$longterm <- Plot$longterm[m]

##Assemble individual data for tanoaks from Mortality file
deadL <- subset(Mortality, Species=="LIDE" & Burned==FALSE & Live.BA.0607 <0.3 & longterm==TRUE)
deadL <- deadL[c("Plot", "TotKilled", "Species", "Live.BA.0607", "Treelist", "Canker06")]
colnames(deadL)[4] <- "liveBA"
deadL$LiveBAlog <- log(deadL$liveBA + 0.01)
deadL$LiveBAc <- (deadL$LiveBAlog - mean(deadL$LiveBAlog) ) / sd( deadL$LiveBAlog)
deadL$plotindex <- group_indices(deadL, Plot)
deadL$TotKilled2[deadL$TotKilled==1] <- 0 ## Change mortality to survival, per reviewer recommendations.
deadL$TotKilled2[deadL$TotKilled==0] <- 1
deadL$TotKilled <- as.numeric(deadL$TotKilled2)
deadL$Canker06 <- as.numeric((deadL$Canker06))
deadL$Canker06 <- (deadL$Canker06 - mean(deadL$Canker06) ) / sd( deadL$Canker06)
deadL[is.na(deadL)]<-0

## Plot-level data for tanoaks:
#Subset to only plots present in tanoak data, scale and center data:
MortPlot <- subset(Plot, BigSurPlot %in% unique(deadL$Plot))
MortPlot <- MortPlot[c("BigSurPlot", "hostCWD", "hostsA", "hostsC", "PFtemp")]
MortPlot$hostCWDc <- (MortPlot$hostCWD - mean(MortPlot$hostCWD) ) / sd( MortPlot$hostCWD)
MortPlot$hostsAc<- (MortPlot$hostsA - mean(MortPlot$hostsA) ) / sd( MortPlot$hostsA)
MortPlot$hostsCc<- (MortPlot$hostsC - mean(MortPlot$hostsC) ) / sd( MortPlot$hostsC)
MortPlot$tempc <- (MortPlot$PFtemp - mean(MortPlot$PFtemp) ) / sd( MortPlot$PFtemp)


##Assemble individual data for redwoods from Mortality file
deadRW <- subset(Mortality, Species=="SESE" & Burned==FALSE & longterm==TRUE & Live.BA.0607 <0.3)
deadRW <- deadRW[c("Plot", "TotKilled", "Species", "Live.BA.0607")]
colnames(deadRW)[4] <- "liveBA"
deadRW$LiveBAlog <- log(deadRW$liveBA + 0.001)
deadRW$LiveBAc <- (deadRW$LiveBAlog - mean(deadRW$LiveBAlog) ) / sd( deadRW$LiveBAlog)
deadRW$plotindex <- group_indices(deadRW, Plot)
deadRW[is.na(deadRW)]<-0
deadRW$TotKilled2[deadRW$TotKilled==1] <- 0
deadRW$TotKilled2[deadRW$TotKilled==0] <- 1
deadRW$TotKilled <- as.numeric(deadRW$TotKilled2)

## Assemble plot level data for redwoods
MortPlot2 <- subset(Plot, BigSurPlot %in% unique(deadRW$Plot))
MortPlot2 <- MortPlot2[c("BigSurPlot", "hostCWD", "hostsA", "hostsC", "PFtemp")]
MortPlot2$hostCWDc <- (MortPlot2$hostCWD - mean(MortPlot2$hostCWD) ) / sd( MortPlot2$hostCWD)
MortPlot2$hostsAc<- (MortPlot2$hostsA - mean(MortPlot2$hostsA) ) / sd( MortPlot2$hostsA)
MortPlot2$hostsCc<- (MortPlot2$hostsC - mean(MortPlot2$hostsC) ) / sd( MortPlot2$hostsC)
MortPlot2$tempc <- (MortPlot2$PFtemp - mean(MortPlot2$PFtemp) ) / sd( MortPlot2$PFtemp)

###### Run full model for LIDE 
## Simulation data for plotting model predictions in gen. quantities -- right now, all set to mean except for dbh effect, but change to plot other effects.
cwd.seq <-rep(0, 300)
dbh.seq <- rep(seq(from = min(deadL$LiveBAc), to = max(deadL$LiveBAc), length.out = 100), 3)
Aseq <- rep(0, 300)
Bseq <- rep(0, 300)

# Assemble data frame for tanoak model fit:
data <- list( TotKilled = deadL$TotKilled, Plot = deadL$plotindex, LiveBAc = deadL$LiveBAc, hostCWD = MortPlot$hostCWDc, N = nrow(deadL), P = length(unique(deadL$plotindex)), climate= MortPlot$tempc, standingA = MortPlot$hostsAc, standingB = MortPlot$hostsCc, Canker=deadL$Canker06, Aseq=Aseq, Bseq=Bseq, cwdseq = cwd.seq, dbhseq = dbh.seq, Nsims= length(cwd.seq))

#fit model for tanoak:
full.mortality.lide.ub <- stan(file = "fullmortalityhost.stan",
                            data = data, 
                            chains=4 , iter=2000 , warmup=1000, thin=1,
                            control=list(adapt_delta = 0.99, stepsize = 0.01))


#traceplot(full.mortality.lide.ub)

plot(full.mortality.lide.ub, pars=c("a", "bdbh", "ba", "bb", "bcwd", "bclim", "bsymp"), outer_level=0.90)
ub.full.mort.LIDE.summary <- round(summary(full.mortality.lide.ub, pars=c("a", "bdbh", "ba", "bb", "bcwd", "bclim", "bsymp", "sigma_plot"), probs = c(0.05, 0.95))$summary, digits=2)

#extract logliklihood from generated quantities to calculate LOO
loglik <- extract_log_lik(full.mortality.lide.ub, parameter_name = "log_lik")
full.lide.loo.ub <- loo(loglik, k_threshold=0.7)


######  Run full model for SESE
cwd.seq2 <- rep(0, 300)
dbh.seq2 <- rep(seq(from = min(deadRW$LiveBAc), to = max(deadRW$LiveBAc), length.out = 100), 3)
Aseq2 <- rep(0, 300)
Bseq2 <- rep(0, 300)
Cseq2 <- rep(0, 300)
data <- list( TotKilled = deadRW$TotKilled, Plot = deadRW$plotindex, LiveBAc = deadRW$LiveBAc, hostCWD = MortPlot2$hostCWDc, N = nrow(deadRW), P = length(unique(deadRW$plotindex)), standingA = MortPlot2$hostsAc, standingB = MortPlot2$hostsCc, disstage= MortPlot2$Pram, Aseq=Aseq2, Bseq=Bseq2, cwdseq = cwd.seq2, dbhseq = dbh.seq2, Nsims= length(cwd.seq2), climate= MortPlot2$tempc)

full.mortality.sese.ub <- stan(file = "fullmortality.stan",
                            data = data, 
                            chains=4 , iter=2000 , warmup=1000, thin=1,
                            control=list(adapt_delta = 0.99, 
                                         stepsize = 0.01))

plot(full.mortality.sese.ub, pars=c("a", "bdbh", "ba", "bb", "bcwd", "bclim"), outer_level=0.90)

#traceplots(full.mortality.sese.ub)
ub.full.mort.SESE.summary <- round(summary(full.mortality.sese.ub, pars=c("a", "bdbh", "ba", "bb", "bcwd", "bclim", "sigma_plot"), probs = c(0.05, 0.95))$summary, digits=2)

#extract log liklihood for loo scores
loglik <- extract_log_lik(full.mortality.sese.ub, parameter_name = "log_lik")
ub.full.sese.loo <- loo(loglik, k_threshold=0.7)


#### 6.2) Resprouting vigor unburned comparison model: #######

##Load the datasets that were formatted by setup file:
setwd("~/Dropbox/Data Analysis/resproutingfireSOD")
Plot = read.csv("datasets/Postfire_Plots_longterm.csv", header=T, na.strings=c("", " "))
Uniquetree = read.csv("datasets/BStrees_longterm.csv", header=T, na.strings=c("", " "))

### Data organization for tanoak:

#TANOAK Data:
#Subset individual data to tanoaks that were alive pre-fire and survived to resprout:
spLIDE = subset(Uniquetree, Species=="LIDE" & Status2006=='L' & Status2010=='L' & Sprouting=="Y" & Burned==FALSE)
spLIDE <- spLIDE[c("Plot", "TotNewStems", "Species", "Burned", "Pr2013", "Canker06", "LBA0607", "DBA10")]
spLIDE$LiveBAc <- (spLIDE$LBA0607 - mean(spLIDE$LBA0607) ) / sd( spLIDE$LBA0607)
spLIDE$DeadBAc <- (spLIDE$DBA10 - mean(spLIDE$DBA10) ) / sd( spLIDE$DBA10)
spLIDE <- transform(spLIDE,plotindex=as.numeric(factor(spLIDE$Plot)))
spLIDE$Canker06 <- as.numeric(as.character(spLIDE$Canker06))
spLIDE$Pr2013 <- as.numeric(as.character(spLIDE$Pr2013))
spLIDE[is.na(spLIDE)]<-0

# Plot-level variables:
# Subset to just plots present in tanoak data & variables needed:
GrowthPlot <- subset(Plot, BigSurPlot %in% unique(spLIDE$Plot))
GrowthPlot <- GrowthPlot[c("BigSurPlot", "FireMort", "DeadBA0607", "DeadBA10", "PFprecip5", "PFtemp5", "Burned_2008")]

#Scale, center, make variables numeric, etc:
GrowthPlot <- transform(GrowthPlot,plotindex=as.numeric(factor(GrowthPlot$BigSurPlot))) #Make plot-level index
GrowthPlot$DeadBA0607 <- as.numeric(as.character(GrowthPlot$DeadBA0607))
GrowthPlot$DeadBA06c <- (GrowthPlot$DeadBA0607 - mean(GrowthPlot$DeadBA0607) ) / sd( GrowthPlot$DeadBA0607)
GrowthPlot$tempc <- (GrowthPlot$PFtemp5- mean(GrowthPlot$PFtemp5) ) / sd( GrowthPlot$PFtemp5)
GrowthPlot$FireMortc<-  (GrowthPlot$FireMort - mean(GrowthPlot$FireMort ) ) / sd( GrowthPlot$FireMort)

# REDWOOD data:
#Subset individual data to redwoods that were alive pre-fire and survived to resprout:
spSESE = subset(Uniquetree, Species=="SESE" & Status2006=='L' & Status2010=='L' & Sprouting=="Y" & Live.BA.0607>0)

spSESE <- spSESE[c("Plot", "TotNewStems", "Species", "LBA0607", "DBA10", "Burned", "Treelist")]
spSESE <- subset(spSESE, Burned==FALSE)
spSESE$LiveBAc <- (spSESE$LBA0607 - mean(spSESE$LBA0607) ) / sd( spSESE$LBA0607)
spSESE$DeadBAc <- (spSESE$DBA10 - mean(spSESE$DBA10) ) / sd( spSESE$DBA10)
spSESE$plotindex <- group_indices(spSESE, Plot)
spSESE[is.na(spSESE)]<-0

## Subset plot-level data to only plots included in redwood individual dataset:
GrowthPlot2 <- subset(Plot, BigSurPlot %in% unique(spSESE$Plot))
GrowthPlot2 <- GrowthPlot2[c("BigSurPlot", "FireMort", "DeadBA0607", "DeadBA10", "PFprecip5", "PFtemp5", "Burned_2008")]

#Scale, center, make variables numeric, etc
GrowthPlot2$FireMortc<-  (GrowthPlot2$FireMort - mean(GrowthPlot2$FireMort ) ) / sd( GrowthPlot2$FireMort)
GrowthPlot2$DeadBA06c <- (GrowthPlot2$DeadBA0607- mean(GrowthPlot2$DeadBA0607) ) / sd( GrowthPlot2$DeadBA0607)
GrowthPlot2$tempc <- (GrowthPlot2$PFtemp- mean(GrowthPlot2$PFtemp) ) / sd( GrowthPlot2$PFtemp)

##### Run Full model (identical to burned resprouting vigor model, except no variable for fire mortality)

## Run model for LIDE:
dbh.seq = rep(seq(from = min(spLIDE$LiveBAc), to = max(spLIDE$LiveBAc), length.out = 100), 3)
top.seq = rep(0, 300)
sod.seq = rep(c(-0.25, 0.5, 1.1), each = 100 )
fire.seq= rep(0, 300)
spLIDE$Pr <- as.numeric(as.character(spLIDE$Pr2013))
Pr.seq = rep(mean(spLIDE$Pr), 300)
Symp.seq = rep(mean(spLIDE$Canker06), 300)

data <- list(TotNewStems = spLIDE$TotNewStems, Plot = spLIDE$plotindex, N = nrow(spLIDE), P = length(unique(GrowthPlot$BigSurPlot)), LiveBAc = spLIDE$LiveBAc, SODdead=GrowthPlot$DeadBA06c, Firedead = GrowthPlot$FireMortc, Pr=spLIDE$Pr, symptom06=spLIDE$Canker06, interact = (GrowthPlot$FireMortc)*(GrowthPlot$DeadBA06c), TopkillBA = spLIDE$DeadBAc, climate=GrowthPlot$tempc, dbhsim = dbh.seq, topkillsim = top.seq , DeadSODsim = sod.seq, DeadFiresim = fire.seq, Nsims = length(dbh.seq), Prsim = Pr.seq, Sympsim = Symp.seq)

UB.full.RS.LIDE <- stan(file = "unburnedresproutinghost.stan",
                     data = data, 
                     chains=4 , iter=2000 , warmup=1000, thin=1,
                     control=list(adapt_delta = 0.99, stepsize = 0.01))


#traceplot(full.RS.LIDE)
plot(UB.full.RS.LIDE, pars=c("a", "bdbh", "btop", "bSOD", "bclim", "bsymp", "bPr" ),  outer_level=0.90)
UB.full.RS.LIDE.summary <- round(summary(UB.full.RS.LIDE, pars=c("a", "bdbh", "btop", "bSOD", "bclim", "bsymp", "bPr", "sigma_plot" ), probs = c(0.05, 0.95))$summary, digits=2)


###### SESE model run:

# Simulated data for plotting model predictions
dbh.seq = rep(seq(from = min(spSESE$LiveBAc), to = max(spSESE$LiveBAc), length.out = 100), 3)
top.seq = rep(0, 300)
sod.seq= rep(0, 300)
fire.seq= rep(0, 300 )

data <- list( TotNewStems = spSESE$TotNewStems, Plot = spSESE$plotindex, N = nrow(spSESE), P = length(unique(spSESE$plotindex)), LiveBAc = spSESE$LiveBAc, SODdead=GrowthPlot2$DeadBA06c, Firedead = GrowthPlot2$FireMortc, climate=GrowthPlot2$tempc, TopkillBA = spSESE$DeadBAc,  dbhsim = dbh.seq, topkillsim = top.seq , DeadSODsim = sod.seq, DeadFiresim = fire.seq, Nsims = length(dbh.seq) )

UB.full.RS.SESE <- stan(file = "unburnedresprouting.stan",
                     data = data, 
                     chains=4 , iter=2000 , warmup=1000, thin=1,
                     control=list(adapt_delta = 0.99, stepsize = 0.01))

#traceplot(full.RS.SESE)
plot(UB.full.RS.SESE, pars=c("a", "bdbh", "btop", "bSOD", "bclim"),  outer_level=0.90)
UB.full.RS.SESE.summary <- round(summary(UB.full.RS.SESE, pars=c("a", "bdbh", "btop", "bSOD", "bclim", "sigma_plot" ), probs = c(0.05, 0.95))$summary, digits=2)


