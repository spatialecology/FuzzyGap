# R Script to support reproducing the analysis in the manuscript:

# "Refined Gap Analysis for Biodiversity Conservation Under Climate Change"

# By: Elham Ebrahimi, Faraham Ahmadzadeh, Asghar Abdoli, Miguel B Araújo, Babak Naimi

##############################################################

# Species Distribution Modelling (SDMs) 
# developed using the sdm R package
#----------------------------------------

library(sdm)

spNames <- c('') # a vector of species names used in the study!

.path <- "~/Documents/Data/GapAnalysis" # path to the working directory!

# Future GCMs (climate models) considered:
gcms <- c( "CanESM5","CNRM-CM6-1" ,"CNRM-ESM2-1","GISS-E2-1-H","INM-CM4-8" ,"MIROC-ES2L", "MIROC6","MPI-ESM1-2-HR")

# future climate change scenarios:
ssps <- c("126","245","370","585")

pr <- rast(paste0(.path,'/bioclim.tif')) # loading predictor variables!

#parallel::detectCores()

for (n in spNames) { # loop over species list
  
  cat('\n Modelling for species No.',i,'which is',n,'is started:')
  sp <- read.csv(paste0(.path,'/species/',n,'.csv')) # read the csv file of species coordinates
  
  sp <- vect(sp,geom=c('lon','lat')) # convert species points to SpatVector!
  
  # preparing sdmdata object (with 1000 background records):
  d <- sdmData(species~.,sp,predictors = pr,bg=list(n=1000,method='gRandom'))
  
  # fitting and evaluating SDMs using 10 modelling methods each with 10 replications:
  m <- sdm(species~., d, methods=c('glmp','gam','brt','rf','mars','cart','fda','maxent','maxlike','bioclim.dismo'),
           replication='boot',n=10,parallelSetting=list(method='parallel',ncore=8))
  
  write.sdm(m,paste0(.path,'/models/MODEL_',n,'.sdm')) # write the model object!
  
  gc()
  
  cat('\n ########### Model for species',n,'is done...!')
  #--------------------------------------- # present-time:
  
  # predictions:
  p <- predict(m, pr, filename=paste0(.path,'/predictions/predicts_present_',n,'.tif'))
  
  # ensemble:
  en <- ensemble(m, p, filename=paste0(.path,'/ensembles/ensemble_present_',n,'.tif'),
                 setting = list(method='weighted',stat='AUC',power=2))
  #---------------
  # Future projections across different years and for different GCMs and SSPs:
  for (y in c(2040,2060,2080,2100)) {
    for (g in gcms) {
      for (s in ssps) {
        prf <- rast(paste0(.path,'/bioclim_',g,'_ssp',s,'_Year-',y,'.tif')) # loading future climate layers
        pf <- predict(m, prf, filename=paste0(.path,'/predictions/predicts_FUTURE_',g,'_ssp',s,'_Year-',y,'__',n,'.tif'))
        enf <- ensemble(m, pf, filename=paste0(.path,'/ensembles/predicts_FUTURE_',g,'_ssp',s,'_Year-',y,'__',n,'.tif'),
                       setting = list(method='weighted',stat='AUC',power=2))
      }
    }
  }
}
#################################################



#################################################
###---- Projection incorporating Dispersal Capacity
#################################################

# function to assign NA to pixels byound the dispersal range (used in .dispersal function):
.getDispersal <- function(r1,r2,distance) {
  w1 <- which(r1[] == 1)
  w <- which(d[] > distance & r2[] == 1)
  r2[w] <- NA
  r2
  
  
}


# assess the dispersal distance per year to identify pixels should be filter out:
.dispersal <- function(pp,y,dis=2000,barrier=NULL) {
  if (nlyr(pp) != length(y)) stop('number of layers should be the same as the length of y (years)')
  
  if (!is.null(barrier)) {
    wb <- which(barrier[] == 1)
  } else wb <- NULL
  
  
  r2 <- pp[[2]]
  r1 <- pp[[1]]
  
  if (!is.null(wb)) r1[wb] <- NA
  
  d <- gridDist(r1,target=1)
  
  .d <- (y[2] - y[1]) * dis
  
  w <- which(r2[] == 1 & d[] > .d)
  
  if (length(w) > 0) r2[w] <- 0
  
  #----
  if (length(y) > 2)
    for (i in 2:(length(y)-1)) {
      .d <- (y[i+1] - y[i]) * dis
      
      r1 <- r2
      r2 <- pp[[i+1]]
      
      if (!is.null(wb)) r1[wb] <- NA
      
      d <- gridDist(r1,target=1)
      
      w <- which(r2[] == 1 & d[] > .d)
      
      if (length(w) > 0) r2[w] <- 0
    }
  
  r2
}
#---------
######
# Example shows how to use .dispersal:

.dis=3000 # Bufo (dispersal capacity for certain species in meter per year)

r <- rast('SPECIES_PROBABILITY_CURRENT_TIME.tif') # output of ensemble for a species for the current time 
rf1 <- rast('SPECIES_PROBABILITY_FUTURE_2040.tif') # output of ensemble for a species for the future time (2040)
rf2 <- rast('SPECIES_PROBABILITY_FUTURE_2060.tif') # 2060
rf3 <- rast('SPECIES_PROBABILITY_FUTURE_2080.tif') # 2080
rf4 <- rast('SPECIES_PROBABILITY_FUTURE_2100.tif') # 2100

# get the best threshold to convert probabilities to presence-absence
th <- getThreshold(m, id='ensemble',opt=2) # use the getThreshold function in the sdm package to get the best threshold for the certain species (m is the output of the sdm function)

rr <- c(r,rf1,rf2,rf3,rf4) # stack all rasters to a single SpatRaster object

# convert to P/A:
pp <- ifel(rr >= th,1,0) 

d1 <- .dispersal(pp[[1:2]],y=c(2000,2040),dis=.dis) # dispersal between the current time and 2040
d2 <- .dispersal(pp[[1:3]],y=c(2000,2040,2060),dis=.dis) # dispersal between the current time and 2060
d3 <- .dispersal(pp[[1:4]],y=c(2000,2040,2060,2080),dis=.dis) # dispersal between the current time and 2080
d4 <- .dispersal(pp[[1:5]],y=c(2000,2040,2060,2080,2100),dis=.dis) # dispersal between the current time and 2100

# repeat this function with the argument dis=0 generates the scenario of No_Dispersal (no dispersal capacity is assumed for the species!)

#########


#################################################
###---- Uncertainty Assessment
#################################################



# variance partitioning using 3-way ANOVA
# factors include: SDMs; GCMs, SSPs



# based on each SDM method and for each GCM and SSP, the richness among all species is calculated.
# here, we get the list of richness rasters for different SDMs, GCMs, and SSPs"

# list of raster files (for Year 2100):
lst <- c()
for (g in gcms) {
  lst <- c(lst,paste0(.path,"/_richness/pr_rich_Disp_",g,"_",ssps,"_2100.tif"))
}

# load all rasters in a single SpatRaster:

r <- rast(lst) # a raster object with all layers (3200= 100 SDM x 8 GCMs x 4 SSPs)


# generate a data.frame with corresponding factors:

.m <- rep(rep(c("glmpoly","gam","brt","rf","mars","fda","svm","maxent","maxlike","bioclim"),each=10),times=32)
.g <- rep(gcms,each=400) # gcms
.s <- rep(rep(ssps,each=100),times=8)

.df <- data.frame(value=NA,gcm=.g,ssp=.s,method=.m)

head(.df)

.df$gcm <- as.factor(.df$gcm)
.df$ssp <- as.factor(.df$ssp)
.df$method <- as.factor(.df$method)



# function for applying 3-way ANOVA at each pixel of the raster object using app function in terra:
.vp <- function(x,.d) {
  if (all(is.na(x))) return(c(NA,NA,NA))
  .d$value <- x
  .aov <- aov(value ~ gcm + ssp + method, data = .d)
  s <- summary(.aov)[[1]]$`Sum Sq`
  s <- s / sum(s)
  s[-length(s)]
}


# ANOVA:
v2100 <- app(r,fun=.vp,.d=.df,cores=10)


names(v2100) <- c('GCMs','SSPs','SDM_Methods')
plot(v2100) # 3 rasters representing spatial distribution of uncertainty linked to GCMs, SSPs, and SDMs
##################


