###################################
##Script to examine Dunwiddies macrofossil data 
##Use gjam fit to 2015 data to recreate MAT from macrofossils
##Created December 5, 2019; last modified ------
##Janneke HilleRisLambers (jhrl@uw.edu)
###################################

##############################################################
###LOAD PACKAGES, READ IN DATA, DEFINE EXPLANATORY VARIBLES### 
##############################################################

##THINGS TO DO
#Reduce tree file to species in macrofossil data
#Rerun models for trees, then predict macrofossil data

##First, read in 2015 tree data
SpeciesRaw2015 <- read.csv("MtRainierTrees2015.csv", header=TRUE)
plot.id1 <- strsplit(as.character(SpeciesRaw2015[,3]),"-")
plot.id2 <- matrix(unlist(plot.id1),ncol=2,byrow=TRUE)
Plotno <- as.numeric(plot.id2[,2])
SpeciesRaw2015 <- cbind(Plotno,SpeciesRaw2015)

#Now, make into matrix - number of trees
SpeciesDat2015 <- tapply(SpeciesRaw2015[,6],list(Plotno,SpeciesRaw2015[,5]), length)
SpeciesDat2015[is.na(SpeciesDat2015)==TRUE] <- 0

#Now, change column headers to 4 letter codes
sp1978_2015 <- read.csv("MtRainierSpeciesTranslation.csv", header=TRUE)
newnames <- dimnames(SpeciesDat2015)[[2]]
for(i in 1:length(newnames)){
  indd <- which(sp1978_2015[,5]==dimnames(SpeciesDat2015)[[2]][i])
  newnames[i] <- as.character(sp1978_2015[indd,3])
}
SpeciesDat2015 <- cbind(as.numeric(dimnames(SpeciesDat2015)[[1]]), SpeciesDat2015)
dimnames(SpeciesDat2015)[[2]] <- c("Plotno",newnames)

#POTR - two species - add together both columns and make new column; order everything alphabetically
POTRind <- which(dimnames(SpeciesDat2015)[[2]]=="POTR")
POTR <- rowSums(SpeciesDat2015[,POTRind])
SpeciesDat2015 <- cbind(SpeciesDat2015[,-POTRind], POTR)
neworder <- order(dimnames(SpeciesDat2015[,-1])[[2]]) + 1
SpeciesDat2015 <- SpeciesDat2015[,c(1,neworder)]

##Next 1978 historic tree data
SpeciesRaw1978 <- read.csv("MtRainierCommunities1978.csv", header=TRUE)

#Combine tree codes: consider when saplings and above (i.e. ABAM2, ABAM3, ABAM4; but not ABAM1)
treesp <- dimnames(SpeciesDat2015)[[2]][-1]
spdelete <- c(); colsadd <- c()

for(i in 1:length(treesp)){
  focsp <- treesp[i]
  foccols <- which(substr(dimnames(SpeciesRaw1978)[[2]],1,4)==focsp)
  if(length(foccols)==0){spdelete <- c(spdelete, i)}
  if(length(foccols)==0){next}
  focpres <- SpeciesRaw1978[,foccols]
  if(length(foccols)==1){
    if(substr(dimnames(SpeciesRaw1978)[[2]][foccols],5,5)!="1"){
      colsadd <- cbind(colsadd, focpres)}
  } 
  
  if(length(foccols)>1){
    focpres <- focpres[,substr(dimnames(focpres)[[2]],5,5)!="1"]
    focpres[focpres<1] <- 0
    colsadd <- cbind(colsadd,rowSums(focpres))}
}

#Add dimnames, turn into presence / absence if desired, merge with plot number
dimnames(colsadd)[[2]] <- treesp[-spdelete]
colsadd[colsadd<=0.2] <-0 #turn any nearby and trace into zero (can't worry about # individuals)
SpeciesDat1978 <- cbind(SpeciesRaw1978[,2],colsadd)
dimnames(SpeciesDat1978)[[2]][1] <- "Plotno"

#Turn both into presence absence
presabs2015 <- SpeciesDat2015[,-1]; presabs2015[presabs2015>0] <-1
presabs2015 <- cbind(SpeciesDat2015[,1], presabs2015)
presabs1978 <- SpeciesDat1978[,-1]; presabs1978[presabs1978>0] <-1
presabs1978 <- cbind(SpeciesDat1978[,1], presabs1978)

# This turns data into presence / absence
SpeciesDat1978 <- presabs1978; dimnames(SpeciesDat1978)[[2]][1] <- "Plotno"
SpeciesDat2015 <- presabs2015; dimnames(SpeciesDat2015)[[2]][1] <- "Plotno"


##Remove columns where species were NOT observed in macrofossil record
#note - in original script, species were only removed based on being present in > 20 plots
focspp <- c("ABAM","ABLA","ABGR","ABPR", "CHNO", "PICO", "PIEN", "PIMO", "PSME","THPL","TSHE","TSME")
keep1978 <- c()
keep2015 <- c()
for(i in 1:length(focspp)){
  tmp1978 <- which(focspp[i]==dimnames(SpeciesDat1978)[[2]])
  if(length(tmp1978)>0){keep1978 <- c(keep1978,tmp1978)}
  tmp2015 <- which(focspp[i]==dimnames(SpeciesDat2015)[[2]])
  if(length(tmp2015)>0){keep2015 <- c(keep2015,tmp2015)}
}

#Prune 1978 and 2015 datasets
SpeciesDat1978 <- SpeciesDat1978[,c(1,keep1978)]
SpeciesDat2015 <- SpeciesDat2015[,c(1,keep2015)]


############################
##Now extract environmental data for 2015 data
###########################
library(raster)  
library(sp)
library(rgdal)

##Now extract elevation and PRISM precipitation data (get MAT, Snow from microclimate models)
EnvDat <- read.csv("MtRainierSiteInfo.csv", header=TRUE)
EnvDat <- EnvDat[,-c(6,9)]

##Now extract LIDAR data
LidarDat <- read.csv("MtRainierLIDAR.csv", header=TRUE)
LidarDat <- LidarDat[LidarDat[,1]=="Franklin",]
LidarDat <- LidarDat[,c(2,4,9,13,18,20)]

##Merge PRISM and LidarDat: Env Dat
EnvDat <- merge(EnvDat, LidarDat,by="Plotno")

##Now extract MAT, GDD and SCD for data
MATfile <- raster("MORAclim_all_summaries_v09/air_temp_summaries/ann_tavg_sum.tif")
GDDfile <- raster("MORAclim_all_summaries_v09/temp_sum_summaries/gdd_avg_2010_2015.tif")
SCDfile <- raster("MORAclim_all_summaries_v09/snow_cover_summaries/scd_avg_wy_2010_2015.tif")

XY <- EnvDat[,4:5] #UTME; UTMN

##Extract spatial data
MAT <- extract(MATfile,XY,method='simple')
GDD <- extract(GDDfile,XY,method='simple')
SCD <- extract(SCDfile,XY,method='simple')

##Now merge above into EnvDat
EnvDat <- cbind(EnvDat[,1:5],MAT,GDD,SCD,EnvDat[,8:15])
rm(MAT); rm(GDD); rm(SCD)

############
###Merge EnvDat and 1978 pres / abs; define X and Y
TotDat1978 <- merge(EnvDat,SpeciesDat1978, by="Plotno")
TotDat1978 <- TotDat1978[is.na(TotDat1978$MAT)==FALSE,]
X <- TotDat1978[,1:10]; Y <- TotDat1978[,17:dim(TotDat1978)[2]]
pltarea <- TotDat1978[,11]

##Load appropriate packages for gjam, MCMC
library(gjam); library(coda)

##Models to fit during model selection; preliminary model fitting suggests MAT best fit
MAT2 <- as.formula(~MAT+I(MAT^2)) #positive coefficient = concave up; negative = concave down
MAT2WINT <- as.formula(~MAT+I(MAT^2)+Wintppt) ##best fitting
MAT2SUM <- as.formula(~MAT+I(MAT^2)+Sumppt)
MAT2WINTSUM <- as.formula(~MAT+I(MAT^2)+Wintppt+Sumppt)
forms <- list(MAT2,MAT2WINT, MAT2SUM, MAT2WINTSUM)

DICs <- c(); minDICs <- c()
fldr <- c("MAT2","MAT2WINT","MAT2SUM","MAT2WINTSUM")

##specify chains, thinning, burnin - used 20K, 10K and 10
samples <- 20000; brnin <- 10000; thn <- 10

## For loop for model selection
for(i in 1:length(forms)){
  bchainsout <- list(); schainsout  <- list(); DICc <- c()
  minDIC <- c()
  
  #specify priors
  #specifing priors; for models with squared terms
  S <- dim(Y)[2] 
  Evars <- attr(terms(forms[[i]]),"term.labels")
  Q <- length(Evars) + 1
  loBeta <- matrix(-Inf,Q,S)         # initialize priors
  hiBeta <- matrix(Inf,Q,S)
  rownames(loBeta) <- rownames(hiBeta) <- c('intercept',Evars)
  
  #set squared term to negative
  sq <- Evars[substr(Evars,1,1)=="I"]
  if(length(sq)>0){hiBeta[sq,] <- 0} # maximum zero for squared term
  
  #set wintppt to negative
  wp <- Evars[Evars[]=="Wintppt"]
  if(length(wp)>0){hiBeta[wp,] <- 0} # maximum zero for wintppt
  
  #set sumppt to negative
  sp <- Evars[Evars[]=="Sumppt"]
  if(length(sq)>0){loBeta[sp,] <- 0} # maximum zero for wintppt
  
  bp <- list(lo = loBeta, hi = hiBeta)
  
   for(j in 1:3){ #one for each chain
    gjamout <- gjam(forms[[i]], X, Y, 
                   modelList=list(ng = samples, burnin = brnin, betaPrior=bp,
                   PREDICTX = F, typeNames='PA', effort = list(pltarea)))
    bchainsout[[j]] <- mcmc(gjamout$chains$bgibbs, start=brnin, end=samples, thin=thn)
    schainsout[[j]] <- mcmc(gjamout$chains$sgibbs, start=brnin, end=samples, thin=thn)
    DICc <- c(DICc, gjamout$fit$DIC)
    minDIC <- c(minDIC, gjamout$fit$DIC)
    if(min(minDIC)==gjamout$fit$DIC){
      pl <- list(sigONLY = F, PLOTY = T, SMALLPLOTS = F, PLOTX = F, 
                 GRIDPLOTS=T, SAVEPLOTS = T, outfolder= fldr[i])
      gjamPlot(gjamout,plotPars = pl)}
    } #, effort = list(pltarea)))} #presence absence
  bchainsmcmc <- mcmc.list(bchainsout[[1]],bchainsout[[2]],bchainsout[[3]])#,bchainsout[[4]])
  schainsmcmc <- mcmc.list(schainsout[[1]],schainsout[[2]],schainsout[[3]])#,schainsout[[4]])
  print(fldr[i])
  print(gelman.diag(bchainsmcmc, multivariate = FALSE))
  print(gelman.diag(schainsmcmc, multivariate = FALSE))
  print(DICc)
  DICs <- c(DICs,mean(DICc)); minDICs <- c(minDICs, min(DICc))
}

##Check DICs; both average and minimum
names(DICs) <- fldr; names(minDICs) <- fldr 
#Note - best model for PA trees = Mat2 or Mat3 (Mat 2 lower mean; Mat3 lowest min)

#############Explore model fit; parameter estimation
##for best fitting model, demonstrate holdouts predict data well
finmodel <- which(DICs==min(DICs)) #This is best fitting model




#########FROM TREE MODEL

####HOLDOUT######
##Now fit best fitting model, holdout 40 samples to check model fit
#gjamout <- gjam(forms[[finmodel]], X, Y, 
#                modelList=list(ng = samples, burnin = brnin, 
#                               PREDICTX = T, holdoutN = 40,  
#                               typeNames='PA', effort = list(pltarea))) #presence absence

#xMu  <- gjamout$modelSummary$xpredMu
#xSd  <- gjamout$modelSummary$xpredSd
#yMu  <- gjamout$modelSummary$ypredMu
#hold <- gjamout$holdoutIndex

#determine predicted and observed for climate in holdouts
#obsenv <- gjamout$x[hold,-1]
#predenv <- xMu[hold,-1]

#Create plot demonstrating how well model fits
#X11(width=8,height=4)
#par(mfrow=c(1,dim(obsenv)[2]-1), pty="s", omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.5), mgp=c(1.25,0.5,0), tck=-0.02) 
#for(j in 1:dim(obsenv)[2]){
#  pred <- predenv[,j]; obs <- obsenv[,j]
#  if(substr(dimnames(predenv)[[2]][j],1,1)=="I"){next}
#  plot(pred, obs, pch=21, bg="lightblue", xlab="predicted", ylab="observed")
#  title(dimnames(predenv)[[2]][j]); abline(0,1)
#  test <- cor.test(obs,pred)
#  mtext(paste("r=", signif(test$estimate,3)), side=3, adj=0.05, line=-1, cex=0.75)
#}

#prespred <- yMu[hold,]; obspred <- gjamout$y[hold,]
#X11(width=8,height=4); par(mfrow=c(2,5), pty="s", omi=c(0,0,0,0), mai=c(0.3,0.3,0.3,0.2), mgp=c(1.25,0.4,0), tck=-0.02)
#for(j in 1:dim(prespred)[2]){
#  obs <- obspred[,j]; pred <- prespred[,j]
#  if(prsabs==1){
#    obs[obs==0] <- "abs"; obs[obs=="1"] <- "pres"
#    plot(as.numeric(pred)~as.factor(obs), col="yellowgreen", xlab="observed", ylab="predicted")
#    title(dimnames(prespred)[[2]][j])}
#  if(prsabs==0){
#    plot(pred,obs, pch=21, bg="yellowgreen", xlab="predicted", ylab="observed")
#    title(dimnames(prespred)[[2]][j])
#    abline(0,1)}
#}


####################
#Read in Macrofossil data
macro <- read.csv("Macrofossils_dunwiddie/DunwiddieMacrofossils.csv", header=TRUE)
#remove PICO1 and PIMO1; change PICO2 and PIMO2 to PICO and PIMO
macro <- macro[macro$Species_code!="PIMO1",]
macro <- macro[macro$Species_code!="PICO1",]
macro[,4] <- substr(macro[,4],1,4)

#Reformat into pres / abs matrix Y2
locdepth <- paste(macro$Location,macro$Depth) #create a vector of location / depth
locdepth2 <- unique(locdepth)
Y2 <- matrix(0, nrow = length(unique(locdepth)), ncol = dim(Y)[2])
dimnames(Y2)[[1]] <- locdepth2
dimnames(Y2)[[2]] <- dimnames(Y)[[2]]

for(i in 1:dim(Y2)[1]){
  tmp <- macro[locdepth[]==dimnames(Y2)[[1]][i],]
  for(j in 1:dim(Y2)[2]){
    tmp2 <- tmp$Percentage[tmp$Species_code==dimnames(Y2)[[2]][j]]
    if(length(tmp2)==0){next}
    if(tmp2>0){Y2[i,j] <- 1}
   }
}

X2 <- matrix(NA, nrow=dim(Y2)[1], ncol=dim(X)[2])
dimnames(X2)[[2]] <- dimnames(X)[[2]]
dimnames(X2)[[1]] <- locdepth2

###Refit model, with climate of past as 'unknown', to check for thermophilization
X3 <- rbind(X2, X)
Y3 <- rbind(Y2, Y)
pltarea3 <- c(rep(max(pltarea), times=dim(Y2)[[1]]), pltarea)
n_macro <- dim(Y2)[[1]]

#Now fit model with unknowns for explanatory variables
gjamout <- gjam(forms[[finmodel]], X3, Y3, 
                  modelList=list(ng = samples, burnin = brnin,
                  typeNames='PA', effort = list(pltarea3))) #presence absence

###########NOW PLOT MACROFOSSIL DATA
#Read in ages, create age depth relationship
agedepth <- read.csv("Macrofossils_dunwiddie/DunwiddieAgeDepth.csv", header=TRUE)
lakes <- unique(agedepth$Location)
lakemod <- list()

##modify depth
Yash <- matrix(NA, nrow=3, ncol=2)
dimnames(Yash)[[1]] <- lakes
dimnames(Yash)[[2]] <- c("depth1","depth2")
Yash[1,] <- c(87,105)
Yash[2,] <- c(40,98)
Yash[3,] <- c(45,64)

X11(width=5,height=9)
par(mfrow=c(3,1))

#Now create an age model for each lake
for(i in 1:length(lakes)){
  lakedepth <- agedepth$Depth[agedepth$Location==lakes[i]]
  gap <- Yash[dimnames(Yash)[[1]]==lakes[i],]
  lakedepth2 <- lakedepth
  lakedepth2[lakedepth2[]>=gap[2]] <-  lakedepth2[lakedepth2[]>=gap[2]] - (gap[2]-gap[1])
  lakeage <- agedepth$Age[agedepth$Location==lakes[i]]
  DAtest <- smooth.spline(lakedepth2, lakeage)
  tmpX <- seq(0, max(lakedepth2))
  tmpY <- predict(DAtest,tmpX)$y
  plot(lakedepth2,lakeage, pch=21, bg="grey", cex=1.5)
  lines(tmpX,tmpY)
  title(lakes[i])
  lakemod[[i]] <- DAtest
}

#Now plot predicted MAT, wintppt and sumppt going back in time, at 3 different ponds
predX <- gjamout$prediction$xpredMu[1:n_macro,]

#Now plot by lake, age; first determine depth and ages
Yash2 <- Yash
Yash2[1,] <- c(85,104)
Yash2[2,] <- c(38,100)
Yash2[3,] <- c(49,64)

lake <- c()
depth <- c(); depth2 <- c()
for(i in 1:n_macro){
  tmp <- strsplit(dimnames(Y2)[[1]][i], split=" ")
  lake <- c(lake, tmp[[1]][1])
  tmpdepth <- as.numeric(tmp[[1]][2])
  depth <- c(depth, tmpdepth)
  tmpgap <- Yash2[dimnames(Yash2)[[1]]==tmp[[1]][1],]
  if(tmpdepth<=tmpgap[1]){depth2 <- c(depth2, tmpdepth)}
  if(tmpdepth>tmpgap[1]){depth2 <- c(depth2, tmpdepth - (tmpgap[2]-tmpgap[1]))}
}


#Now plot time and three climatic variables
X11(width=6, height=4)
par(mfrow=c(1,1),omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.5), tck=-0.02, mgp=c(1.2,0.5,0))
lakeages <- c()

for(i in 1:length(lakes)){
  lakedat <- predX[lake==lakes[i],]
  lakedat <- as.data.frame(lakedat)
  lakedepth <- depth2[lake==lakes[i]]
  lakeage <- predict(lakemod[[i]],lakedepth)$y
  lakeages <- c(lakeages, lakeage)
  print(lakes[i]); print(cor.test(lakeage, lakedat$MAT))
  if(i==1){plot(lakeage, lakedat$MAT, xlab="Age (Years BP)", ylab="predicted MAT (C)", 
                xlim=c(0,6500), ylim=c(2,8), type="l", col="firebrick1", lty=1, lwd=2)}
  if(i>1){lines(lakeage, lakedat$MAT,col="firebrick1", lty=i, lwd=2)}
}

X11(width=6, height=4)
par(mfrow=c(1,1),omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.5), tck=-0.02, mgp=c(1.2,0.5,0))

for(i in 1:length(lakes)){
  lakedat <- predX[lake==lakes[i],]
  lakedat <- as.data.frame(lakedat)
  lakedepth <- depth2[lake==lakes[i]]
  lakeage <- predict(lakemod[[i]],lakedepth)$y
  print(lakes[i]); print(cor.test(lakeage, lakedat$Sumppt))
  if(i==1){plot(lakeage, lakedat$Sumppt, xlab="Age (Years BP)", ylab="Predicted Summer precip (cm)", 
                xlim=c(0,6500), ylim=c(15,25), type="l", col="green3", lty=1, lwd=2)}
  if(i>1){lines(lakeage, lakedat$Sumppt,col="green3", lty=i, lwd=2)}
}



###plot macrofossil community data in CCA (old analysis from community turnover)
## Read in data from 1976-1980 census
Legacy_Trees <- read.csv("Macrofossils_dunwiddie/LegacyPlots_1978.csv", header=TRUE) #Franklin_Communityl is Legacy data (463 plots censused between 1976-1980)

## Change data into a matrix with relative abundance, species (columns) x site (rows)
plots <- unique(Legacy_Trees[,1]); plots <- plots[order(plots)]
spp <- unique(Legacy_Trees[,2]); spp <- spp[order(spp)]
Franklin_Communityl <- matrix(0, nrow=length(plots), ncol=length(spp))
dimnames(Franklin_Communityl) <- list(plots, spp)

## For loop to fill matrix
for(i in 1:length(plots)){
  plottrees <- Legacy_Trees[Legacy_Trees[,1]==plots[i],]
  
  for(j in 1:length(spp)){
    basp <- 0 #reset basal area to zero
    plottreessp <- plottrees[plottrees[,2]==spp[j],]
    if(dim(plottreessp)[1]==0){next} #if not entries, go to next species
    
    ## trees counted in dbh classes, so add basal area by midpoint
    for(k in 3:14){
      numtrees <- plottreessp[1,k]
      if(is.na(numtrees)==TRUE){next}
      basp <- basp + numtrees*(pi*(((-25 + k*10)/2)^2))
    }
    
    ## DBH of large trees measured exactly, add basal area of these species
    nbig <- as.numeric(plottreessp[15])
    if(is.na(nbig)==FALSE){
      for(k in 1:nbig){
        basp <- basp + pi*(plottreessp[15+k]/2)^2
        basp <- basp[1,1]
      }
    }
    
    ## now put final basal area, divided by plot size into matrix
    Franklin_Communityl[i,j] <- basp / plottreessp[1,23]
    #print(j); print(dim(Franklin_Communityl))
  }
}

## turn matrix into relative abundance
for(i in 1:dim(Franklin_Communityl)[1]){
  Franklin_Communityl[i,] <- Franklin_Communityl[i,] / sum(Franklin_Communityl[i,])
}

#Reduced matrix for species in macrofossil data
Franklin_Communityl2 <- Franklin_Communityl[,c(1,3,4,6,7,8,10,13,15,16,17)]

## Read in environmental data
Legacy_clim <- read.csv("Macrofossils_dunwiddie/Legacyplots_clim.csv", header=TRUE)
Legacy_clim <- Legacy_clim[order(Legacy_clim[,1]),] #order by plot


## First load packages
library(vegan)

## Now identify explanatory variables (Response variable is Franklin_Communityl)
MAT <- Legacy_clim[,10]
Wintppt <- Legacy_clim[,11]
Sumppt <-Legacy_clim[,12]

## Now analysis, and various goodness of fit
CCA_Franklinl <- cca(Franklin_Communityl ~ MAT + Wintppt + Sumppt, scale=FALSE)
goodness(CCA_Franklinl, display=c("species","sites"), statistic="explained")
vif.cca(CCA_Franklinl)
anova(CCA_Franklinl, by="mar", perm=1000) #assesses significance of explanatory vars, type III (sig wit all others present)
anova(CCA_Franklinl, by="axis", perm=1000) #assesses signficance of each constrained axis
spenvcor(CCA_Franklinl) #species environmental correlation - not a great measure of fit though
goodcca_site<-goodness(CCA_Franklinl, display="sites", model="CCA", statistic="explained") #goodness of fit for sites
print(quantile(goodcca_site[,2], probs=c(0.025,0.5,0.975))) # #print median and range of site goodness of fit (CCA1 and CCA2 only)
goodcca_spp<-goodness(CCA_Franklinl, display="species", model="CCA", statistic="explained") #goodness of fit for sites
print(goodcca_spp) # #print median and range of site goodness of fit (CCA1 and CCA2 only)

## Reduced community, rescale and remove 
#MAT2 <- MAT[rowSums(Franklin_Communityl2)>0]
#Wintppt2 <- Wintppt[rowSums(Franklin_Communityl2)>0]
#Sumppt2 <- Sumppt[rowSums(Franklin_Communityl2)>0]
#Franklin_Communityl2 <- Franklin_Communityl2[rowSums(Franklin_Communityl2)>0,]
#for(i in 1:dim(Franklin_Communityl2)[1]){
#  Franklin_Communityl2[i,] <- Franklin_Communityl2[i,] / sum(Franklin_Communityl2[i,])
#}
#CCA_Franklinl <- cca(Franklin_Communityl2 ~ MAT2 + Wintppt2 + Sumppt2, scale=FALSE)


###################################################################################
##### Make figures of plots, species, climate (from legacy data set): Fig. 5A #####
###################################################################################

## Extract necessary values from CCA output for plotting
CCAres <- plot(CCA_Franklinl, type="n", scaling=3) #just to extract values
PlotScores <- CCAres$sites
SpeciesScores <- CCAres$species
EnvScores <- CCAres$biplot

## Make a two panel plot, first shows the data
X11(width=5, height=5)
par(mfrow=c(1,1), tck=-0.02, omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.4), mgp=c(1.25,0.5,0))
xrng <- c(min(floor(PlotScores[,1])), max(ceiling(PlotScores[,1])))
yrng <- c(min(floor(PlotScores[,2])), max(ceiling(PlotScores[,2])))
plot(PlotScores[,1], PlotScores[,2], xlab="CCA1", ylab="CCA2", pch=16, col="grey", 
     cex=0.9, xlim=xrng, ylim=c(-7,7))

## Add horizontal, vertical line
abline(h=0, lty=3)
abline(v=0, lty=3)

## Plot species scores, but only of select species
SpeciesScores2<-c()
speciestoinclude<-c("ABAM","ABLA","ABPR", "CHNO", "PICO", "PIMO", "PSME","TSHE","TSME")

for(j in 1:length(speciestoinclude)){
  SpeciesScores2<-rbind(SpeciesScores2,SpeciesScores[dimnames(SpeciesScores)[[1]]==speciestoinclude[j],])
}

dimnames(SpeciesScores2)[[1]]<-speciestoinclude

## now add to plot
text(SpeciesScores2[,1], SpeciesScores2[,2], labels=dimnames(SpeciesScores2)[[1]], cex=0.75)

## next plot environmental scores
arrowcol <- c("red","blue","yellowgreen")
for(i in 1:dim(EnvScores)[1]){
  arrows(0,0,EnvScores[i,1], EnvScores[i,2], length=0.1, angle=45, code=2, col=arrowcol[i], lwd=2.5)
}


## Change macrofossil data to same format at Franklin_Communityl
macro_communityl <- matrix(0, nrow=dim(Y2)[1], ncol=dim(Franklin_Communityl)[2])
dimnames(macro_communityl)[[1]] <- dimnames(Y2)[[1]]
dimnames(macro_communityl)[[2]] <- dimnames(Franklin_Communityl)[[2]]

## Now fill macrofossil matrix
for(i in 1:dim(macro_communityl)[1]){
  tmp <- macro[locdepth[]==dimnames(macro_communityl)[[1]][i],]
  for(j in 1:dim(macro_communityl)[2]){
    tmp2 <- tmp$Percentage[tmp$Species_code==dimnames(macro_communityl)[[2]][j]]
    if(length(tmp2)==0){next}
    if(tmp2>0){macro_communityl[i,j] <- tmp2}
  }
}


## turn into data frame
macro_communityl <- data.frame(macro_communityl)

# second data frame - averaged by 1000 years (use lake ages); except for 1st 500 years
macro_communityl2 <- c()
ages <- seq(1000,5000, by=1000)
ages <- c(0,250,500,ages,6500)
lake2 <- c()
lbl <- c()

for(i in 1:length(lakes)){
  lake_communityl <- macro_communityl[lake==lakes[i],]
  tmpages <- lakeages[lake==lakes[i]]
  for(j in 1:(length(ages)-1)){
    lake_communitylage <- lake_communityl[tmpages[]>=ages[j]&tmpages[]<ages[j+1],]
    lake_communityavg <- colMeans(lake_communitylage)
    macro_communityl2 <- rbind(macro_communityl2,lake_communityavg)
    tmplbl <- paste(lakes[i],ages[j])
    lbl <- c(lbl, tmplbl)
    lake2 <- c(lake2, as.character(lakes[i]))
    print(j); print(tmplbl)

  }
}

dimnames(macro_communityl2)[[1]] <- lbl

## Now use CCA results (from spatial analysis) to predict climate space of each year of composition
testchange <- predict(CCA_Franklinl, newdata=macro_communityl, type="wa")

#Now plot macrofossil trajectories by lakes (backdrop of Franklin data)
## Make a two panel plot, first shows the data
X11(width=5, height=5)
par(mfrow=c(1,1), tck=-0.02, omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.4), mgp=c(1.25,0.5,0))
xrng <- c(min(floor(PlotScores[,1])), max(ceiling(PlotScores[,1])))
yrng <- c(min(floor(PlotScores[,2])), max(ceiling(PlotScores[,2])))
plot(PlotScores[,1], PlotScores[,2], xlab="CCA1", ylab="CCA2", 
     pch=16, col="grey", cex=0.9, xlim=xrng, ylim=c(-7,7))

## Add horizontal, vertical line
abline(h=0, lty=3)
abline(v=0, lty=3)

## Add macrossil data all data
for(i in 1:length(lakes)){
  lakecca <- testchange[lake==lakes[i],]
  lines(lakecca[,1], lakecca[,2])
  points(lakecca[1,1], lakecca[1,2], pch=21, bg="plum1", cex=2)
  points(lakecca[dim(lakecca)[1],1], lakecca[dim(lakecca)[1],2], pch=24, bg="mediumpurple4", cex=2)
}



## Remake, a bit easier to follow
X11(width=5, height=5)
par(mfrow=c(1,1), tck=-0.02, omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.4), mgp=c(1.25,0.5,0))
xrng <- c(min(floor(PlotScores[,1])), max(ceiling(PlotScores[,1])))
yrng <- c(min(floor(PlotScores[,2])), max(ceiling(PlotScores[,2])))
plot(PlotScores[,1], PlotScores[,2], xlab="CCA1", ylab="CCA2", 
     pch=16, col="grey", cex=0.9, xlim=xrng, ylim=c(-7,7))

## Add horizontal, vertical line
abline(h=0, lty=3)
abline(v=0, lty=3)

##Color code by mid holocene (3000-6000), late holocene (0-3000)
#ageint <- c(-1,200,3000, 6500)
ageint <- c(-1,3000, 6500)
agecols <- c("plum1", "mediumpurple4")

## Add macrossil data all data
for(i in 1:2){
  lakecca <- testchange[lakeages>=ageint[i]&lakeages<ageint[i+1],]
  points(lakecca[,1], lakecca[,2], pch=21, bg=agecols[i], cex=2)
}

## next plot environmental scores
arrowcol <- c("red","blue","yellowgreen")
for(i in 1:dim(EnvScores)[1]){
  arrows(0,0,EnvScores[i,1], EnvScores[i,2], length=0.1, angle=45, code=2, col=arrowcol[i], lwd=2.5)
}

  
##Dissimilarity
alldat <- rbind(macro_communityl, Franklin_Communityl2)
diss <- vegdist(alldat, binary=TRUE, method="bray") #add binary = TRUE for pres / abs
diss <- as.matrix(diss)

#get mean dissimilarities by age
dissmacro <- diss[1:dim(macro_communityl)[1], (dim(macro_communityl)[1]+1):dim(diss)[1]]
X11(width=8, height=6)
par(mfrow=c(4,6), omi=c(0,0,0,0), mai=c(0.4,0.4,0.4), mgp=c(1.25,0.5,0))
for(i in 1:dim(dissmacro)[1]){
  #hist(dissmacro[i,])
  #title(dimnames(dissmacro)[[1]][i])
  print(dimnames(dissmacro)[[1]][i])
  print(mean(dissmacro[i,])); print(min(dissmacro[i,])); print(max(dissmacro[i,]))
}

for(i in 1:dim(diss)[1]){
  print(min(diss[i,])); print(max(diss[i,]))
}