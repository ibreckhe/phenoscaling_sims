####Script to summarise landscape phenology sims.
##Author: Ian Breckheimer
##Created: 15 September 2019
##Updated: 14 February 2020

####Sets up workspace.
library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)

setwd("~/code/phenoscaling_sims")

####Loads data
days <- seq(90,306,by=1)
numeric_cols <- length(days)
colspec <- paste(c(rep("i",numeric_cols),rep("c",3),rep("i",1)),collapse="")
colnames <- c(paste("Day",days,sep=""),"Species","Year","Replicate","Cellsize")
samples_all <- readRDS("./output/scaling_allspp_2007_2017_daily.Rdata")

##Function to calculate beginning and end of flowering using the "area under the curve" method. ##
obs_intervals <- function(x,y,threshold=0.1586553){
  
  stopifnot(length(x)==length(y))
  
  # Temporary variables.
  y_ord <- y[order(x)]
  x_ord <- x[order(x)]
  
  # Gets the bin width from the first pred interval.
  bin_width <- x_ord[2] - x_ord[1]
  
  # Assume all NA predictions are zero
  y_ord[is.na(y_ord)] <- 0
  
  # Total area under the curve.
  total_area <- sum(y_ord*bin_width)
  
  # Computes cumulative proportions in both directions
  cumprop_up <- cumsum(y_ord)*bin_width/total_area
  cumprop_down <- rev(cumsum(rev(y_ord))*bin_width)/total_area
  
  # Finds the indices of the first and last values greater than 0.158
  lwr_index <- min(which(cumprop_up >= threshold))
  upr_index <- max(which(cumprop_down >= threshold))
  
  # Finds the corresponding values of x.
  lwr_bound <- x_ord[lwr_index]
  upr_bound <- x_ord[upr_index]
  max_x <- x_ord[which.max(y_ord)]
  
  bounds <- c(lwr_bound,max_x,upr_bound)
  names(bounds) <- c("lwr_bound","max","upr_bound")
  
  # Output
  return(bounds)
}

##computes start and end of flowering for all samples.
bounds_mat <- matrix(NA,nrow=nrow(samples_all),ncol=3)

for(i in 1:nrow(samples_all)){
  if(i/1000==round(i/1000)){
    print(paste("Sample",i,"of",nrow(samples_all)))
  }
  bounds_mat[i,] <- obs_intervals(x=days,y=as.numeric(samples_all[i,1:217]),threshold=0.001)
}

bounds_df <- as.data.frame(bounds_mat)
bounds_df$Species <- samples_all$Species
bounds_df$Year <- samples_all$Year
bounds_df$Replicate <- samples_all$Replicate
bounds_df$Cellsize <- samples_all$Cellsize
colnames(bounds_df)[1:3] <- c("DOY_start","DOY_max","DOY_end")

##Removes values at the start and end of the season.
bounds_df$DOY_max[bounds_df$DOY_max==90] <- NA
bounds_df$DOY_max[bounds_df$DOY_max==306] <- NA
bounds_df$DOY_start[bounds_df$DOY_start==90] <- NA
bounds_df$DOY_start[bounds_df$DOY_start==306] <- NA
bounds_df$DOY_end[bounds_df$DOY_start==90] <- NA
bounds_df$DOY_end[bounds_df$DOY_end==306] <- NA



bounds_df$length <- bounds_df$DOY_end - bounds_df$DOY_start


##Summarises by grain size for each species.
bounds_spp <- bounds_df %>% group_by(Species, Cellsize, Year, Replicate) %>% 
                            summarise(mean_DOY_start = mean(DOY_start),
                                      sd_DOY_start = sd(DOY_start,na.rm=TRUE),
                                      mean_DOY_peak =  mean(DOY_max,na.rm=TRUE),
                                      sd_DOY_peak = sd(DOY_max,na.rm=TRUE),
                                      mean_DOY_end = mean(DOY_end,na.rm=TRUE),
                                      sd_DOY_end = sd(DOY_end,na.rm=TRUE),
                                      mean_DOY_length = mean(length,na.rm=TRUE),
                                      sd_DOY_length = sd(length,na.rm=TRUE))
bounds_spp$scenario <- paste(bounds_spp$Species,bounds_spp$Year,bounds_spp$Replicate,sep="-")


##Subsets only 2m grain sizes.
bounds_2m <- filter(bounds_df, Cellsize==2)

bounds_2m_sum <- bounds_2m %>% group_by(Species,Year, Replicate) %>% 
  summarise(mean_DOY_start_2m = mean(DOY_start,na.rm=TRUE),
            sd_DOY_start_2m = sd(DOY_start,na.rm=TRUE),
            mean_DOY_peak_2m =  mean(DOY_max,na.rm=TRUE),
            sd_DOY_peak_2m = sd(DOY_max,na.rm=TRUE),
            mean_DOY_end_2m = mean(DOY_end,na.rm=TRUE),
            sd_DOY_end_2m = sd(DOY_end,na.rm=TRUE),
            mean_DOY_length_2m = mean(length,na.rm=TRUE),
            sd_DOY_length_2m = sd(length,na.rm=TRUE))

bounds_2m_sum$scenario <- paste(bounds_2m_sum$Species,bounds_2m_sum$Year,bounds_2m_sum$Replicate,sep="-")

##Joins resampled results with 2m summaries and computes anomalies.
bounds_spp_anom <- left_join(bounds_spp,bounds_2m_sum,by=c("scenario","Species","Year","Replicate"))

bounds_spp_anom$DOY_start_anom <- bounds_spp_anom$mean_DOY_start - bounds_spp_anom$mean_DOY_start_2m
bounds_spp_anom$DOY_peak_anom <- bounds_spp_anom$mean_DOY_peak - bounds_spp_anom$mean_DOY_peak_2m
bounds_spp_anom$DOY_end_anom <- bounds_spp_anom$mean_DOY_end - bounds_spp_anom$mean_DOY_end_2m
bounds_spp_anom$length_anom <- bounds_spp_anom$mean_DOY_length - bounds_spp_anom$mean_DOY_length_2m

anom_sum <- bounds_spp_anom %>% group_by(Cellsize) %>% summarise(start_anom_q10=quantile(DOY_start_anom,probs=0.1,na.rm=TRUE),
                                                                 start_anom_q25=quantile(DOY_start_anom,probs=0.25,na.rm=TRUE),
                                                                 start_anom_q50=quantile(DOY_start_anom,probs=0.50,na.rm=TRUE),
                                                                 start_anom_q75=quantile(DOY_start_anom,probs=0.75,na.rm=TRUE),
                                                                 start_anom_q90=quantile(DOY_start_anom,probs=0.90,na.rm=TRUE),
                                                                 peak_anom_q10=quantile(DOY_peak_anom,probs=0.1,na.rm=TRUE),
                                                                 peak_anom_q25=quantile(DOY_peak_anom,probs=0.25,na.rm=TRUE),
                                                                 peak_anom_q50=quantile(DOY_peak_anom,probs=0.50,na.rm=TRUE),
                                                                 peak_anom_q75=quantile(DOY_peak_anom,probs=0.75,na.rm=TRUE),
                                                                 peak_anom_q90=quantile(DOY_peak_anom,probs=0.90,na.rm=TRUE),
                                                                 end_anom_q10=quantile(DOY_end_anom,probs=0.1,na.rm=TRUE),
                                                                 end_anom_q25=quantile(DOY_end_anom,probs=0.25,na.rm=TRUE),
                                                                 end_anom_q50=quantile(DOY_end_anom,probs=0.50,na.rm=TRUE),
                                                                 end_anom_q75=quantile(DOY_end_anom,probs=0.75,na.rm=TRUE),
                                                                 end_anom_q90=quantile(DOY_end_anom,probs=0.90,na.rm=TRUE),
                                                                 length_anom_q10=quantile(length_anom,probs=0.1,na.rm=TRUE),
                                                                 length_anom_q25=quantile(length_anom,probs=0.25,na.rm=TRUE),
                                                                 length_anom_q50=quantile(length_anom,probs=0.50,na.rm=TRUE),
                                                                 length_anom_q75=quantile(length_anom,probs=0.75,na.rm=TRUE),
                                                                 length_anom_q90=quantile(length_anom,probs=0.90,na.rm=TRUE))
                                                                 
                                                                 

p1 <- ggplot(anom_sum)+
  geom_abline(aes(intercept=0,slope=0),linetype="dotted")+
  geom_linerange(aes(x=Cellsize,ymin=start_anom_q25,ymax=start_anom_q75),color="black",lwd=1.5)+
  geom_linerange(aes(x=Cellsize,ymin=start_anom_q10,ymax=start_anom_q90),color="black")+
  geom_point(aes(x=Cellsize,y=start_anom_q50),shape=21,fill="white",color="black",size=3)+
  geom_linerange(aes(x=Cellsize,ymin=peak_anom_q25,ymax=peak_anom_q75),color="grey30",lwd=1.5)+
  geom_linerange(aes(x=Cellsize,ymin=peak_anom_q10,ymax=peak_anom_q90),color="grey30")+
  geom_point(aes(x=Cellsize,y=peak_anom_q50),fill="white",color="grey30", size=3, shape=22)+
  geom_linerange(aes(x=Cellsize,ymin=end_anom_q25,ymax=end_anom_q75),color="grey50",lwd=1.5)+
  geom_linerange(aes(x=Cellsize,ymin=end_anom_q10,ymax=end_anom_q90),color="grey50")+
  geom_point(aes(x=Cellsize,y=end_anom_q50),fill="white",color="grey50",size=3, shape=23)+
  geom_linerange(aes(x=Cellsize,ymin=length_anom_q25,ymax=length_anom_q75),color="grey70",lwd=1.5)+
  geom_linerange(aes(x=Cellsize,ymin=length_anom_q10,ymax=length_anom_q90),color="grey70")+
  geom_point(aes(x=Cellsize,y=length_anom_q50),fill="white",color="grey70",size=3,shape=24)+
  #scale_shape_manual(values=c(21,22,23,24),labels=c("first flower", "peak flower","last flower","season length"))+
  scale_x_log10("Spatial grain (m, log scale)")+
  scale_y_continuous("Measurement Anomaly (Days)")+
  coord_cartesian(ylim=c(-20,20),xlim=c(1,1000), clip="on")+
  annotation_logticks(sides="b", mid = unit(0.1, "cm"),)+
  theme_bw()+
  theme(panel.grid=element_blank())

svg("./figs/spatial_grain_anomanly.svg",width=3.5,height=3)
p1
dev.off()

#### Subsamples at different time intervals to examine the effects of temporal grain.

temp_grains <- c(1,2,3,4,6,8,12,17)
bounds_grain_list <- list()

for(i in 1:length(temp_grains)){
  print(paste("Temporal grain", temp_grains[i], "days"))
  samples_2m <- filter(samples_all, Cellsize==2)
  samples_num <- samples_num[,1:217]
  grain_seq <- seq(min(days),max(days),by=temp_grains[i])
  day_bin <- days %in% grain_seq
  samples_grain <- samples_num[,day_bin]
  
  ##computes start and end of flowering for all samples.
  bounds_mat_grain <- matrix(NA,nrow=nrow(samples_grain),ncol=3)
  
  for(j in 1:nrow(samples_grain)){
    if(i/1000==round(i/1000)){
      print(paste("Sample",j,"of",nrow(samples_grain)))
    }
    bounds_mat_grain[j,] <- obs_intervals(x=grain_seq,y=as.numeric(samples_grain[j,]),threshold=0.001)
  }
  
  
  bounds_df_grain <- as.data.frame(bounds_mat_grain)
  bounds_df_grain$Species <- samples_2m$Species
  bounds_df_grain$Year <- samples_2m$Year
  bounds_df_grain$Replicate <- samples_2m$Replicate
  bounds_df_grain$Cellsize <- samples_2m$Cellsize
  colnames(bounds_df_grain)[1:3] <- c("DOY_start","DOY_max","DOY_end")
  
  ##Removes values at the start and end of the season.
  bounds_df_grain$DOY_max[bounds_df_grain$DOY_max==90] <- NA
  bounds_df_grain$DOY_max[bounds_df_grain$DOY_max==306] <- NA
  bounds_df_grain$DOY_start[bounds_df_grain$DOY_start==90] <- NA
  bounds_df_grain$DOY_start[bounds_df_grain$DOY_start==306] <- NA
  bounds_df_grain$DOY_end[bounds_df_grain$DOY_start==90] <- NA
  bounds_df_grain$DOY_end[bounds_df_grain$DOY_end==306] <- NA
  
  bounds_df_grain$length <- bounds_df_grain$DOY_end - bounds_df_grain$DOY_start
  
  bounds_spp_grain <- bounds_df_grain %>% group_by(Species, Cellsize, Year, Replicate) %>% 
    summarise(mean_DOY_start = mean(DOY_start),
              sd_DOY_start = sd(DOY_start,na.rm=TRUE),
              mean_DOY_peak =  mean(DOY_max,na.rm=TRUE),
              sd_DOY_peak = sd(DOY_max,na.rm=TRUE),
              mean_DOY_end = mean(DOY_end,na.rm=TRUE),
              sd_DOY_end = sd(DOY_end,na.rm=TRUE),
              mean_DOY_length = mean(length,na.rm=TRUE),
              sd_DOY_length = sd(length,na.rm=TRUE))
  bounds_spp_grain$scenario <- paste(bounds_spp_grain$Species,bounds_spp_grain$Year,bounds_spp_grain$Replicate,sep="-")
  
  ##Joins to 2m estimates to compute anomalies.
  bounds_spp_anom_grain <- left_join(bounds_spp_grain,bounds_2m_sum,by=c("scenario","Species","Year","Replicate"))
  
  bounds_spp_anom_grain$DOY_start_anom <- bounds_spp_anom_grain$mean_DOY_start - bounds_spp_anom_grain$mean_DOY_start_2m
  bounds_spp_anom_grain$DOY_peak_anom <- bounds_spp_anom_grain$mean_DOY_peak - bounds_spp_anom_grain$mean_DOY_peak_2m
  bounds_spp_anom_grain$DOY_end_anom <- bounds_spp_anom_grain$mean_DOY_end - bounds_spp_anom_grain$mean_DOY_end_2m
  bounds_spp_anom_grain$length_anom <- bounds_spp_anom_grain$mean_DOY_length - bounds_spp_anom_grain$mean_DOY_length_2m
  bounds_spp_anom_grain$temp_grain <- temp_grains[i]
  bounds_grain_list[[i]] <- bounds_spp_anom_grain
}

bounds_grain_df <- bind_rows(bounds_grain_list)

anom_sum_grain <- bounds_grain_df %>% group_by(temp_grain) %>% summarise(start_anom_q10=quantile(DOY_start_anom,probs=0.1,na.rm=TRUE),
                                                                 start_anom_q25=quantile(DOY_start_anom,probs=0.25,na.rm=TRUE),
                                                                 start_anom_q50=quantile(DOY_start_anom,probs=0.50,na.rm=TRUE),
                                                                 start_anom_q75=quantile(DOY_start_anom,probs=0.75,na.rm=TRUE),
                                                                 start_anom_q90=quantile(DOY_start_anom,probs=0.90,na.rm=TRUE),
                                                                 peak_anom_q10=quantile(DOY_peak_anom,probs=0.1,na.rm=TRUE),
                                                                 peak_anom_q25=quantile(DOY_peak_anom,probs=0.25,na.rm=TRUE),
                                                                 peak_anom_q50=quantile(DOY_peak_anom,probs=0.50,na.rm=TRUE),
                                                                 peak_anom_q75=quantile(DOY_peak_anom,probs=0.75,na.rm=TRUE),
                                                                 peak_anom_q90=quantile(DOY_peak_anom,probs=0.90,na.rm=TRUE),
                                                                 end_anom_q10=quantile(DOY_end_anom,probs=0.1,na.rm=TRUE),
                                                                 end_anom_q25=quantile(DOY_end_anom,probs=0.25,na.rm=TRUE),
                                                                 end_anom_q50=quantile(DOY_end_anom,probs=0.50,na.rm=TRUE),
                                                                 end_anom_q75=quantile(DOY_end_anom,probs=0.75,na.rm=TRUE),
                                                                 end_anom_q90=quantile(DOY_end_anom,probs=0.90,na.rm=TRUE),
                                                                 length_anom_q10=quantile(length_anom,probs=0.1,na.rm=TRUE),
                                                                 length_anom_q25=quantile(length_anom,probs=0.25,na.rm=TRUE),
                                                                 length_anom_q50=quantile(length_anom,probs=0.50,na.rm=TRUE),
                                                                 length_anom_q75=quantile(length_anom,probs=0.75,na.rm=TRUE),
                                                                 length_anom_q90=quantile(length_anom,probs=0.90,na.rm=TRUE))

p2 <- ggplot(anom_sum_grain)+
  geom_abline(aes(intercept=0,slope=0),linetype="dotted")+
  geom_linerange(aes(x=temp_grain,ymin=start_anom_q25,ymax=start_anom_q75),color="black",lwd=1.5)+
  geom_linerange(aes(x=temp_grain,ymin=start_anom_q10,ymax=start_anom_q90),color="black")+
  geom_point(aes(x=temp_grain,y=start_anom_q50),shape=21,fill="white",color="black",size=3)+
  geom_linerange(aes(x=temp_grain,ymin=peak_anom_q25,ymax=peak_anom_q75),color="grey30",lwd=1.5)+
  geom_linerange(aes(x=temp_grain,ymin=peak_anom_q10,ymax=peak_anom_q90),color="grey30")+
  geom_point(aes(x=temp_grain,y=peak_anom_q50),fill="white",color="grey30", size=3, shape=22)+
  geom_linerange(aes(x=temp_grain,ymin=end_anom_q25,ymax=end_anom_q75),color="grey50",lwd=1.5)+
  geom_linerange(aes(x=temp_grain,ymin=end_anom_q10,ymax=end_anom_q90),color="grey50")+
  geom_point(aes(x=temp_grain,y=end_anom_q50),fill="white",color="grey50",size=3, shape=23)+
  geom_linerange(aes(x=temp_grain,ymin=length_anom_q25,ymax=length_anom_q75),color="grey70",lwd=1.5)+
  geom_linerange(aes(x=temp_grain,ymin=length_anom_q10,ymax=length_anom_q90),color="grey70")+
  geom_point(aes(x=temp_grain,y=length_anom_q50),fill="white",color="grey70",size=3,shape=24)+
  #scale_shape_manual(values=c(21,22,23,24),labels=c("first flower", "peak flower","last flower","season length"))+
  scale_x_log10("Temporal Grain (Days, log scale)")+
  scale_y_continuous("Measurement Anomaly (Days)")+
  coord_cartesian(ylim=c(-20,20),xlim=c(1,17),clip="on")+
  annotation_logticks(sides="b",mid = unit(0.1, "cm"),)+
  theme_bw()+
  theme(panel.grid=element_blank())

svg("./figs/temporal_grain_anomanly.svg",width=4,height=3.5)
p2
dev.off()

####Creates combination plot.
library(gridExtra)
svg("./figs/spatial_temporal_grain_anomaly.svg")
grid.arrange(p1,p2,nrow=1)
dev.off()
