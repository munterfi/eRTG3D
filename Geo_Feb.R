#Loadings#====
library(sf)
library(eRTG3D)
library(move)
library(zoo)
library(plyr)
library(dplyr)
library(lubridate)
library(bcpa)
tryCatch(setwd("/Users/george/Dropbox/gitHub/eRTG3D_2019/"))
# Loading from RData save
load("/Users/george/Dropbox/gitHub/eRTG3D_2019/bird.Rdata") #sigouraki
load('/Users/george/Dropbox/gitHub/eRTG3D_2019/water_bodies.Rdata')
# BIRD_MIGR from eRTG files ====
data <- read.csv("/Users/george/Dropbox/phd/code/online_R_cran_Workspace/stat_RTG_on_R/Final_Exp/eRTG_2018/Input/MPIO white stork lifetime tracking data (2013-2014).csv", stringsAsFactors = F)
# splitting 
data_split <- data %>% group_split(tag.local.identifier)
bird_raw <- data_split[[57]] #52 #53 #54 #55 #56
plot(bird_raw$location.long, bird_raw$location.lat, type = "l")
for (i in 50:69) {
  bird_raw <- data_split[[i]] #52 #53 #54 #55 #56
  plot(bird_raw$location.long, bird_raw$location.lat, type = "l", main = i)
}


plot(bird_raw$location.long, bird_raw$location.lat, type = "l")

steps = 10
preproc_any_traj <- function(bird, steps){
  #Preprocesses the trajectory for Merlin in 3 steps
  # a) adds 3rd dimension
  # b) corrects the column sequence as needed by eRTG3D 
  # c) reprojects the whole thing
  # d) subsamples (thinnens) the traj by the "steps" given. 
  if (class(bird)[1] == "Move") {
    bird$z = 0
    #changing the sequence as requested: 
    bird <- data.frame(x = bird$x, y = bird$y, z = bird$z, time = bird$time)
    # Transform CRS
    # using projected CRS as WGS 84 lat/lon was not working properly 
    t_projected <- transformCRS.3d(bird, fromCRS = 4326, toCRS = 32638)
    t_projected <- track.properties.3d(t_projected)
    t_projected$z <- 0
    ###########=========================
    #add the time
    t_projected$time <- bird$time
    #harmonize time diff
    high_res_train_list <- track.split.3d(track = t_projected, timeLag = as.numeric(diff(t_projected$time))) #consider lag = 5min + tolerance)
    ###########=========================
    #high resolution segments, time regulated
    high_res_train_list_reg <- lapply(X = high_res_train_list, FUN = function(x) {
      n <- seq(1, length(x$x), steps)
      x[n, ]
    })
  }
  if (class(bird)[1] == "tbl_df") {
    bird$z = 0
    #changing the sequence as requested: 
    bird$time <- ymd_hms(bird$timestamp)
    bird <- data.frame(x = bird$location.long, y = bird$location.lat, z = bird$z, time = bird$time)
    bird <- na.omit(bird) #removing the funny reconrds.
    # Transform CRS
    # using projected CRS as WGS 84 lat/lon was not working properly 
    t_projected <- transformCRS.3d(bird, fromCRS = 4326, toCRS = 32638)
    t_projected <- track.properties.3d(t_projected)
    t_projected$z <- 0
    ###########=========================
    #add the time
    t_projected$time <- bird$time
    #harmonize time diff
    high_res_train_list <- track.split.3d(track = t_projected, timeLag = as.numeric(diff(t_projected$time))) #consider lag = 5min + tolerance)
    ###########=========================
    #high resolution segments, time regulated
    high_res_train_list_reg <- lapply(X = high_res_train_list, FUN = function(x) {
      n <- seq(1, length(x$x), steps)
      x[n, ]
    })
  }
  plot2d(high_res_train_list_reg, DEM = NULL)
  return(high_res_train_list_reg)
}

preproc_traj_M <- function(bird, steps){
  #Preprocesses any move or dataframe trajectory for Merlin's package in 3 steps
  # a) adds 3rd dimension
  # b) corrects the column sequence as needed by eRTG3D 
  # c) reprojects the whole thing
  # d) subsamples (thinnens) the traj by the "steps" given. 
  bird$z = 0
  #changing the sequence as requested: 
  bird <- data.frame(x = bird$x, y = bird$y, z = bird$z, time = bird$time)
  # Transform CRS
  # using projected CRS as WGS 84 lat/lon was not working properly 
  t_projected <- transformCRS.3d(bird, fromCRS = 4326, toCRS = 32638)
  t_projected <- track.properties.3d(t_projected)
  t_projected$z <- 0
  ###########=========================
  #add the time
  t_projected$time <- bird$time
  #harmonize time diff
  high_res_train_list <- track.split.3d(track = t_projected, timeLag = as.numeric(diff(t_projected$time))) #consider lag = 5min + tolerance)
  ###########=========================
  #high resolution segments, time regulated
  high_res_train_list_reg <- lapply(X = high_res_train_list, FUN = function(x) {
    n <- seq(1, length(x$x), steps)
    x[n, ]
  })
  plot2d(high_res_train_list_reg, DEM = NULL)
  return(high_res_train_list_reg)
}






end_points_sel <- function(t_projected, migration = TRUE){
  if (migration == T){
    start <- Reduce(c, t_projected[[1]][1, 1:3])
    #searching for a point really south
    df_ls <- (do.call(rbind.data.frame, t_projected))
    i <- which(df_ls$y == min(df_ls$y))
    end <- as.numeric(df_ls[i, 1:3])
  }
  else {
    start <- Reduce(c, t_projected[[1]][1, 1:3])
    end <- Reduce(c, t_projected[[length(t_projected)]][nrow(t_projected[[length(t_projected)]]), 1:3])
    
  }
  a <- as.data.frame(rbind(start, end))
  colnames(a) <- c("x", "y", "z")
  
  #ToDo: needs to be corrected. 
  #(lapply(t_projected, function(x) min(x$y))
  #which.min(t_projected$'246'$y)
  #index <- which.min(sapply(t_projected,function(x){min(x$y)}))
  return(a)
}

#SIM
#2D option B 
# P Prob

t_projected <- preproc_traj_M(bird_migr, steps)#high_res_train_list_reg1

# or

#t_projected <- preproc_any_traj(bird_raw, steps)#high_res_train_list_reg
P <- get.section.densities.3d(t_projected, DEM = NULL)
#using it only for gathering starting and ending point
bird <-  bird_raw #bird_migr  #bird_raw  

# Boundary conditions
#sum(sapply(t_projected, nrow)) + as.integer(222/10)
#TODO: match the filtered start/end point with the original traj's ones. 
n.locs <- as.integer(difftime(bird$time[nrow(bird)], bird$time[1], units = "mins" ) /50 )#time 
n.locs <- as.integer(difftime(bird$timestamp[nrow(bird)], bird$timestamp[1], units = "mins" ) /3500 )#timeinterval

start <- end_points_sel(t_projected)[1,]
end <- end_points_sel(t_projected)[2,]


# 
# start <- Reduce(c, t_projected[[1]][1, 1:3])
# #end <- Reduce(c, t_projected[[n.locs]], 1:3])
# #the plotting stops on the first section - needs to be updated
# end <- Reduce(c, t_projected[[length(t_projected)]][nrow(t_projected[[length(t_projected)]]), 1:3])
a0 <- t_projected[[1]]$a[1]
g0 <- t_projected[[1]]$g[1]

#plotting all
plot2d(do.call("rbind", t_projected[start:end], ), BG = water_bodies)

#alt plot
t <- do.call("rbind", t_projected)
plot(water_bodies)
points(t$x, t$y)

# Q Probs: increase steps allowed by 10%
sim.locs <- round(n.locs * 1.1)
#sim.locs <- n.locs
uerw <- sim.uncond.3d(n.locs = sim.locs*1000, a0 = a0, g0 = g0,
                      start = start, densities = P)
Q <- qProb.3d(uerw, sim.locs, multicore = TRUE)

# Simulate CERWs
cerws <- n.sim.cond.3d(n.sim = 500, multicore = T,
                       sim.locs, start = start, end = end, a0 = a0, g0 = g0,
                       DEM = NULL, BG = water_bodies, densities = P, qProbs = Q)

cerws <- filter.dead.ends(cerws)
plot2d(t_projected, cerws)
#plot2d(t_projected)
plot2d(t_projected, cerws, BG = water_bodies)


# start the evaluation part. 



































# eRTG files ====
data <- read.csv("/Users/george/Dropbox/phd/code/online_R_cran_Workspace/stat_RTG_on_R/Final_Exp/eRTG_2018/Input/MPIO white stork lifetime tracking data (2013-2014).csv", stringsAsFactors = F)
# splitting 
data_split <- data %>% group_split(tag.local.identifier)
bird_raw <- data_split[[6]]