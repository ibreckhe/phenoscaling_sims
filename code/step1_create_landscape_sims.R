####Script to create landscape phenology simulations.
##Author: Ian Breckheimer
##Created: September 3rd, 2019
##Modified: February 13th, 2020

#library("ggmcmc")
library("readr")
library("tidyr")
library("dplyr")
library("stringr")
library("raster")
library("foreach")
library("doParallel")


##Project Path
proj_path <- "~/code/phenoscaling_sims"
setwd(proj_path)

##Model coefficients for virtual species.
coef_means <- as.matrix(read_csv("./data/virtual_species_coef_means.csv"))

##Covariate scaling and centering coefficients
doy_scale <- 36.224
doy_center <- 198.6824

##Function to predict number of flowers given covariates and fit coefficients.
predict_fun <- function(x,spp,coefs){
  opt <- coefs[spp,1] + x[,2] * coefs[spp,2] + x[,3] * coefs[spp,3] + x[,4] * coefs[spp,4] + x[,6] * coefs[spp,5]
  width <- exp(coefs[spp,6] + x[,2] * coefs[spp,7] + x[,4] * coefs[spp,8] + x[,5] * coefs[spp,9] + x[,6] * coefs[spp,10]) * -1
  height <- coefs[spp,11] + x[,2] * coefs[spp,12] + x[,4] * coefs[spp,13] + x[,6] * coefs[spp,14]
  mu <- exp(width * (x[,1] - opt)^2 + height)
  return(mu)
}

predict_pheno_rast <- function(xrast_path = "output/Fit_RF_snow_gdd_2017_rep_1.tif",
                               xdoy = seq(130,200,by=10),
                               xdoy_scale = doy_scale,
                               xdoy_center = doy_center,
                               xspp = c(2,4,18)){
  out <- as.list(rep(NA,length(xspp)))
  
  xrast <- readAll(brick(xrast_path))
  sdd <- values(xrast[[1]])
  sdd_sq <- sdd^2
  gdd <- values(xrast[[2]])
  sm <- rep(0,length(gdd))
  sm_sq <- rep(0,length(gdd))
  doy_sc <- ((xdoy-xdoy_center)/xdoy_scale)
  
  #par(mfrow=c(length(spp),length(doy)))
  for(i in 1:length(xspp)){
    print(paste("Computing predictions for species",xspp[i]))
    spp <- rep(xspp[i],length(gdd))
    yrasts <- stack()
    for(j in 1:length(xdoy)){
      print(paste("Day",xdoy[j]))
      xmat <- cbind(doy_sc[j],sdd,sdd_sq,gdd,sm,sm_sq)
      yrast <- xrast[[1]]
      values(yrast) <- predict_fun(x=xmat,
                                   spp=spp[i],
                                   coefs=coef_means)
      names(yrast) <- paste("DOY",xdoy[j],sep="")
      yrasts <- stack(yrasts,yrast)
    }
    ybrick <- brick(yrasts)
    out[[i]] <- ybrick
  }
  names(out) <- xspp
  return(out)
}

test <- predict_pheno_rast(xrast_path = "output/Fit_RF_snow_gdd_2017_rep_1.tif",
                   xdoy = seq(135,170,by=5),
                   xdoy_scale = doy_scale,
                   xdoy_center = doy_center,
                   xspp = c(1))
plot(test[[1]])

##Sums flowers in lower-resolution bins.
aggregate_pheno_rast <- function(ybrick,
                                 factors,
                                 agg_fun){
  downsampled <- as.list(rep(NA,length(factors)+1))
  downsampled[[1]] <- ybrick
  for(i in 1:length(factors)){
    print(paste("Aggregating rasters with factor",factors[i]))
    downsampled[[(i+1)]] <- aggregate(ybrick,fact=factors[i],fun=agg_fun,expand=TRUE)
  }
  return(downsampled)
}


####Loops through every species, year, and replicate, creating fully scaled and downscaled landscapes.
spp <- 1:45
spp_names <- paste("species",spp,sep="")
years <- 2007:2017
reps <- 1:10
days <- seq(90,306,by=1)
overwrite <- TRUE

##Parallel backend.
library(foreach)
library(doParallel)
cores <- parallel::detectCores()
cl <- makeCluster(floor(cores/2))
registerDoParallel(cl)

foreach(i=1:length(spp),.packages=c("raster")) %dopar% {
  print(paste("Computing Predictions for ", spp_names[i],"(", i,"of",length(spp)))
  for (j in 1:length(years)){
    print(paste("Year",years[j]))
    for (k in 1:length(reps)){
      print(paste("Repicate",k))
      rast_name <- paste("Fit_RF_snow_gdd_",years[j],"_rep_",reps[k],".tif",sep="")
      rast_path <- paste("output/", rast_name,sep="")
      out_name <- paste("./scratch/landscape_preds/",spp_names[i],"-",rast_name,sep="")
      if (file.exists(out_name) & overwrite == FALSE){
        print(paste("Output", out_name, "exists, skipping..."))
        next
      }else{
        fullres <- predict_pheno_rast(xrast_path = rast_path,
                                      xdoy = days,
                                      xdoy_scale = attributes(scale(traind$doy))$'scaled:scale',
                                      xdoy_center = attributes(scale(traind$doy))$'scaled:center',
                                      xspp = spp[i])[[1]]
        writeRaster(fullres,filename=out_name,datatype="INT2U",overwrite=TRUE)
      }
    }
  }
}
stopCluster(cl)

##Resamples all datasets and extracts a random subset of time-series.

landscape_files <- list.files("./scratch/landscape_preds/",pattern=".tif$", full.names = TRUE)
#landscape_files <- landscape_files[1:32]
landscape_info <- str_split_fixed(landscape_files,pattern="//",n=2)[,2]
landscape_spp <- str_split_fixed(landscape_info,pattern="-",n=2)[,1]
landscape_other <- str_split_fixed(landscape_info,pattern="-",n=2)[,2]
landscape_year <- str_split_fixed(landscape_other,pattern="_",n=7)[,5]
landscape_rep <- gsub(".tif","",str_split_fixed(landscape_other,pattern="_",n=7)[,7])

landscape_attribs <- data.frame(Species=landscape_spp,
                                Year=landscape_year,
                                Replicate=landscape_rep)

#cores <- parallel::detectCores()
cl <- makeCluster(24)
registerDoParallel(cl)


##Takes about 45 minutes with 32 processes.
start_time <- Sys.time()
samples_all <- foreach(i=1:length(landscape_files),.packages=c("dplyr","raster"),.combine="rbind") %dopar% {
  #RevoUtilsMath::setMKLthreads(1)
  land_brick <- brick(landscape_files[i])
  agg_fact <- c(2,4,8,16,32,64,128,
                256,512)
  aggregated <- aggregate_pheno_rast(ybrick=land_brick,
                       factors=agg_fact,
                       agg_fun=sum)
  sample_list <- as.list(rep(NA,length(aggregated)))
  for (j in 1:length(aggregated)){
    samp_size=100
    if(ncell(aggregated[[j]]) > samp_size){
      sampled <- as.data.frame(sampleRandom(aggregated[[j]],size=samp_size))
    }else{
      sampled <- as.data.frame(getValues(aggregated[[j]]))
    }
    
    colnames(sampled) <- paste("CountDay",days,sep="")
    sampled$Species <- landscape_attribs$Species[i]
    sampled$Year <- landscape_attribs$Year[i]
    sampled$Replicate <- landscape_attribs$Replicate[i]
    sampled$Cellsize <- agg_fact[j]
    sample_list[[j]] <- sampled
  }
  sample_frame <- bind_rows(sample_list)
  sample_frame
}
stopCluster(cl)
end_time <- Sys.time()
duration <- end_time - start_time
duration

samples_all$Cellsize[is.na(samples_all$Cellsize)] <- 1024
samples_tbl <- as.tbl(samples_all)

write_csv(samples_tbl,path="./output/scaling_allspp_2007_2017_daily.csv")
saveRDS(samples_tbl,file="./output/scaling_allspp_2007_2017_daily.Rdata")
#write_csv(samples_tbl,path="~/Dropbox/scaling_allspp_2007_2017_daily.csv")
#saveRDS(samples_tbl,file="~/Dropbox/scaling_allspp_2007_2017_daily.Rdata")
