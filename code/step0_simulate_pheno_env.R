####Script to simulate realistic surfaces of snow melt date, gdd, and soil moisture based on field data.
##Author:Ian Breckheimer
##Created: 3 March 2019
##Updated: 13 February 2020

library(readr)
library(dplyr)
library(raster)
library(RandomFields)
library(svglite)
library(gstat)

project_dir <- "~/code/phenoscaling_sims/"
#project_dir <- "~/host_data/pheno_scaling/"
setwd(project_dir)

##Imports RMBL coordinates and environmental covariates.
locations <- read_csv(paste(project_dir,"data/rmbl_phenology_centers.csv",sep=""))
env <- read_csv(paste(project_dir,"data/rmbl_spatial_pheno_clim_zerofilled_2007_2017.csv",sep=""))[,-c(4,5)]
env$BB_date_lastsnow <- as.numeric(format(as.Date(env$BB_date_lastsnow,
                                              format="%d-%b-%y"),format="%j"))

##Summarises covariates by year.
env_sum <- env %>% group_by(year,plot,plotyear) %>% summarise_all(.funs="mean")
env <- as.data.frame(left_join(env_sum,locations,by=c("plot"="Plot")))
coordinates(env) <- ~longitude+latitude
proj4string(env) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
env_utm <- spTransform(env,CRSobj=CRS("+proj=utm +zone=13 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
#plot(env_utm)
env_utm@data$UTM_X <- coordinates(env_utm)[,1]
env_utm@data$UTM_Y <- coordinates(env_utm)[,2]

##Creates an empty raster template
env <- env_utm@data
plot(env$UTM_X,env$UTM_Y)
xmin <- round(min(env$UTM_X - 300)/2)*2
ymin <- round((min(env$UTM_Y)-38)/2)*2

grid_coords <- expand.grid(x=seq(xmin,xmin+1024,by=2),
                           y=seq(ymin,ymin+1024,by=2))
grid_coords$z <- runif(n=nrow(grid_coords),0,1)
grid_raster <- rasterFromXYZ(grid_coords,res=c(2,2),
                             crs=CRS("+proj=utm +zone=13 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
grid_coords$x_scaled <- (grid_coords$x - mean(grid_coords$x)) / sd(grid_coords$x)
grid_coords$y_scaled <- (grid_coords$y - mean(grid_coords$y)) / sd(grid_coords$y)


writeRaster(grid_raster,filename=paste(project_dir,"data/rastertemplate.tif",sep=""),
            overwrite=TRUE)
plot(grid_raster)
points(env_utm)

##Fits univariate spatial models using gstat
env$snow_scaled <- scale(env$HOBO_doy_lastsnow)
env$sfgdd_scaled <- scale(env$HOBO_sfgdd_90d)
x_center <- mean(env$UTM_X)
y_center <- mean(env$UTM_Y)
env$X_UTM_jit <- (env$UTM_X-x_center)/1000 + rnorm(length(env$UTM_X),0,0.01)
env$Y_UTM_jit <- (env$UTM_Y-y_center)/1000 + rnorm(length(env$UTM_Y),0,0.01)


years <- 2007:2017
snow_mean <- rep(NA,length(years))
snow_var <- rep(NA,length(years))
snow_nugs <- rep(NA,length(years))
snow_ranges <- rep(NA,length(years))
snow_sills <- rep(NA,length(years))
gdd_mean <- rep(NA,length(years))
gdd_var <- rep(NA,length(years))
gdd_nugs <- rep(NA,length(years))
gdd_ranges <- rep(NA,length(years))
gdd_sills <- rep(NA,length(years))
covar_gamma1 <- rep(NA,length(years))
covar_gamma2 <- rep(NA,length(years))
covar_gamma3 <- rep(NA,length(years))

for(i in 1:length(years)){
  env_yr <- subset(env,year == years[i])
  snow_mean[i] <- mean(env_yr$snow_scaled)
  snow_var[i] <- sd(env_yr$snow_scaled)^2
  gdd_mean[i] <- mean(env_yr$sfgdd_scaled)
  gdd_var[i] <- sd(env_yr$sfgdd_scaled)^2
  cv <- variogram(snow_scaled+sfgdd_scaled~1,width=0.05, cutoff=0.7,
                  locations=~X_UTM_jit+Y_UTM_jit,data=env_yr,covariogram=TRUE)
  covar_gamma1[i] <- cv$gamma[length(cv$gamma)]
  covar_gamma2[i] <- cv$gamma[1]
  covar_gamma3[i] <- cv$gamma[2]
  plot(cv$dist,cv$gamma,cex=cv$np/max(cv$np))
  
  vg <- variogram(snow_scaled~1,width=0.01, cutoff=0.7,
                  locations=~X_UTM_jit+Y_UTM_jit,data=env_yr,covariogram=FALSE)
  
  plot(vg$dist,vg$gamma)
  vg_fit <- fit.variogram(vg,vgm(psill=0.5, model="Mat", nugget=0, range=0.05,kappa=5),
                          fit.sills=c(TRUE,TRUE))
  plot(vg,vg_fit)
  snow_nugs[i] <- vg_fit[1,2]
  snow_ranges[i] <- vg_fit[2,3]
  snow_sills[i] <- vg_fit[2,2]
  
  
  vg2 <- variogram(sfgdd_scaled~1,width=0.01, cutoff=0.6,
                   locations=~X_UTM_jit+Y_UTM_jit,data=env_yr,covariogram=FALSE)
  plot(vg2)
  vg2_fit <- fit.variogram(vg2,vgm(psill=0.1, model="Mat", nugget=0.1, range=0.1,kappa=1),
                           fit.sills=c(TRUE,TRUE))
  gdd_nugs[i] <- vg2_fit[1,2]
  gdd_ranges[i] <- vg2_fit[2,3]
  gdd_sills[i] <- vg2_fit[2,2]
  plot(vg2,vg2_fit)
}

##Corrects a few missed fits.
snow_ranges[snow_ranges>1] <- 0.7
snow_ranges[snow_ranges<=0] <- 1e-4
snow_sills[snow_sills>1] <- 0.6
snow_sills[snow_sills<0.05] <- 0.1
gdd_ranges[gdd_ranges>1] <- 0.7
gdd_ranges[gdd_ranges<0.01] <- 0.05
gdd_sills[gdd_sills>1] <- 0.6
covar_gamma <- rowMeans(cbind(covar_gamma1,covar_gamma2,covar_gamma3))


##Simulates using the Parsimonious Multivariate Whittle Matern Model 
##Gneiting, T., Kleiber, W., Schlather, M. (2010) Matern covariance functions for multivariate
##random fields JASA

RFoptions(cPrintlevel=2,
          linesimustep=0.001,
          modus_operandi="easygoing")

rho_year <- matrix(nc=2, c(0.3, 0.1, 0.1, 0.3))
nu_year <- c(0.1,0.3)
model <- RMparswmX(nudiag=nu_year, rho=rho_year)
plot(model)
year_reps <- 10
years <- 2007:2017

##Creates snow disapearance day and growing degree-day landscapes.
for (i in 1:length(years)){
  print(paste("Simulating landscape for year",years[i]))
  rho_year <- matrix(nc=2, c(snow_sills[i], covar_gamma[i]/4, covar_gamma[i]/4, gdd_sills[i]))
  nu_year <- c(snow_ranges[i],gdd_ranges[i])
  model <- RMparswmX(nudiag=nu_year, rho=rho_year)
  #plot(model)
  
  x.seq <-  grid_coords_sc[,1]
  y.seq <- grid_coords_sc[,2]
  
  for(j in 1:year_reps){
    print(paste("Replicate",j))
    z <- RFsimulate(model = model, x=x.seq, y=y.seq)
    z1 <- z@data$variable1
    z2 <- z@data$variable2
    zr1 <- grid_raster
    values(zr1) <- z1 + snow_mean[i]
    zr2 <- grid_raster
    values(zr2) <- z2 + gdd_mean[i]
    zr_brick <- brick(zr1,zr2)
    zr_brick <- brick(zr1,zr2,filename=paste("./output/Fit_RF_snow_gdd_",years[i],"_rep_",j,".tif",sep=""),
                      overwrite=TRUE)
    names(zr_brick) <- c("Last Snow DOY","GDD90")
    svglite(paste("./figs/sdd_gdd_rf_",years[i],".svg",sep=""),width=6,height=2.6)
    plot(zr_brick,zlim=c(-3,3))
    dev.off()
    #plot(zr1,zr2)
  }
}