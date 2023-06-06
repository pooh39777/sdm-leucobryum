#Get all the libraries
library(maptools)
library(raster) #read shape file
library(tidyverse)
library(rgdal)
library(maps) #map data
library(broom)
library(dismo) #species distribution
library(rJava)
library(randomForest)

#1. Main map
#1.1 Get map data
thai0 <- getData("GADM", country = "THA", level = 0)
thai1 <- getData("GADM", country = "THA", level = 1)

#1.2 ggplot2 map
ggplot() +
  geom_polygon(aes(long, lat, group = group), data = thai0, 
               fill = "#99d8c9", color = "#2ca25f") +
  coord_equal() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#E9FFFF")) +
  labs(x = "Longtitude (°E)", y = "Latitude (°N)") 

#2. Species occurrence
#2.1 Import from GBIF
gbif <- gbif("leucobryum","aduncum")
##clean data
gbif.clean <- gbif %>% 
  filter(fullCountry == "Thailand") %>% 
  select(lon, lat) %>% 
  na.omit()

#2.2 Import from .csv file
csv <- read_csv("Leucobryum-Specimens.csv")
##clean data
csv.clean <- csv %>% 
  filter(Name == "Leucobryum aduncum var. aduncum") %>% 
  select(lon, lat) %>% 
  na.omit()

gbif.csv <- rbind(gbif.clean,csv.clean)

#2.3 Plot species occurrence
ggplot() +
  geom_polygon(aes(long, lat, group = group), data = thai1, fill = "#99d8c9", color = "#2ca25f") +
  geom_point(aes(lon, lat), data = gbif.csv, color = "red") +
  coord_equal() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#E9FFFF")) +
  labs(x = "Longtitude (°E)", y = "Latitude (°N)") 

#3. Species Distribution Modeling
#3.1 Get bioclim data
##Get the available bioclim data from the worldclim dataset. See more at <https://www.worldclim.org/data/bioclim.html>
bioclim <- getData('worldclim', var='bio', res=2.5)
##Get elevation from SRTM 90 m resolution data between -60 and 60 latitude
elevation <- getData('alt', country = "THA", mask = FALSE)

##Crop bioclim data to Thailand
bioclim.thai.crop <- crop(bioclim, thai0)
bioclim.thai <- mask(bioclim.thai.crop, thai0)
plot(bioclim.thai[[1]])

##Crop elevation data to Thailand
elevation.thai.crop <- crop(elevation, thai0)
elevation.thai.0 <- mask(elevation.thai.crop, thai0)

#extent(elevation.thai) <- extent(bioclim.thai) #Didn't use
##Set bioclim.thai and elevation.thai be equal column number (stack preparation)
elevation.thai <- resample(elevation.thai.0, bioclim.thai, method = "bilinear")
plot(elevation.thai[[1]])

#3.2 Build the models
##select only lon and lat of species occurrence
gbif.csv.xy <- select(gbif.csv, lon, lat)
gbif.csv.xy

##split into training and testing dataset
numrow <- nrow(gbif.csv.xy)
sample <- sample(numrow, 0.3 * numrow)
gbif.csv.xy.train <- gbif.csv.xy[-sample, ]
gbif.csv.xy.test <- gbif.csv.xy[sample, ]

predictors <- stack(bioclim.thai[[c(1:19)]], elevation.thai) #choose bio1-19 and alt

##create random background points
bg <- randomPoints(predictors, 1000)
numrow.bg <- nrow(bg)
sample.bg <- sample(numrow.bg, 0.3 * numrow.bg)
bg.train <- bg[-sample.bg, ]
bg.test <- bg[sample.bg, ]
ssb <- ssb(gbif.csv.xy.test, bg.test, gbif.csv.xy.train)
ssb[,1] / ssb[,2]
#Spatial sorting bias (SSB), If there is no SSB this value should be 1, close to zero indicating that SSB is very strong.

#Subsample to remove SSB
i <- pwdSample(gbif.csv.xy.test, bg.test, gbif.csv.xy.train, n=1, tr=0.33)
gbif.csv.xy.test.pwd <- gbif.csv.xy.test[!is.na(i[,1]), ]
bg.test.pwd <- bg.test[na.omit(as.vector(i)), ]
sb2 <- ssb(gbif.csv.xy.test.pwd, bg.test.pwd, gbif.csv.xy.train)
sb2[1]/sb2[2]

##build the models
bc.model <- bioclim(x = predictors, p = gbif.csv.xy.train)
bc.model
me.model <- maxent(x = predictors, p = gbif.csv.xy.train)
me.model

colnames(bg.train) <- c("lon","lat")
train <- rbind(gbif.csv.xy.train, bg.train)
pb.train <- c(rep(1, nrow(gbif.csv.xy.train)), rep(0, nrow(bg.train)))
envtrain.0 <- raster::extract(predictors, train)
envtrain <- data.frame(cbind(pa=pb.train, envtrain.0))
glm.model <- glm(pa ~ ., family = gaussian(link = "identity"), data=envtrain)
glm.model
#predict(glm.model)

rf.model <- randomForest(pa ~ ., data = envtrain)
rf.model

#3.3 see the response curves
#3.3.1 Bioclim model
#response(bc.model) #Figure margins too large
#response(bc.model, var = "bio1")
response(bc.model, var = c(1:5))
response(bc.model, var = c(6:10))
response(bc.model, var = c(11:15))
response(bc.model, var = c(16:20))

##3.3.2 MaxEnt model
#response(me.model) #Figure margins too large
#response(me.model, var = "bio1")
response(me.model, var = c(1:5))
response(me.model, var = c(6:10))
response(me.model, var = c(11:15))
response(me.model, var = c(16:20))

#3.4 Plot the predicted data onto the map
##Basic plot
predict.presence.bc <- dismo::predict(object = bc.model, x = predictors)
plot(predict.presence.bc)
points(lat ~ lon, data = gbif.csv.xy, col = "red")

predict.presence.me <- dismo::predict(object = me.model, x = predictors)
plot(predict.presence.me)
points(lat ~ lon, data = gbif.csv.xy, col = "red")

ext <- extent(90, 110, 5, 22)
predict.presence.glm <- predict(predictors, glm.model, ext=ext)
plot(predict.presence.glm)
points(lat ~ lon, data = gbif.csv.xy, col = "red")

predict.presence.rf <- predict(predictors, rf.model, ext=ext)
plot(predict.presence.rf)
points(lat ~ lon, data = gbif.csv.xy, col = "red")

##Plot in ggplot2
###Turn raster into df
predict.presence.bc.df <- as.data.frame(predict.presence.bc, xy = TRUE) %>% 
  mutate(layer2 = ifelse(layer=="NaN",NA,layer)) #Clean layer for transparent background
predict.presence.me.df <- as.data.frame(predict.presence.me, xy = TRUE) %>% 
  mutate(layer2 = ifelse(layer=="NaN",NA,layer)) #Clean layer for transparent background
predict.presence.glm.df <- as.data.frame(predict.presence.glm, xy = TRUE) %>% 
  mutate(layer2 = ifelse(layer=="NaN",NA,layer)) #Clean layer for transparent background
predict.presence.rf.df <- as.data.frame(predict.presence.rf, xy = TRUE) %>% 
  mutate(layer2 = ifelse(layer=="NaN",NA,layer)) #Clean layer for transparent background

###Bioclim model
midpoint <- (max(predict.presence.bc.df$layer2, na.rm = TRUE)
             -min(predict.presence.bc.df$layer2, na.rm = TRUE))/2
ggplot() +
  geom_raster(data = predict.presence.bc.df, aes(x,y, fill = layer2)) +
  geom_polygon(aes(long, lat, group = group), data = thai1, fill = NA, color = "gray80", size = 0.1) +
  geom_point(data = gbif.csv.xy, aes(lon, lat), color = "red", shape = 1) +
  scale_fill_gradient2("probability", 
                      high = "red", mid = "yellow", midpoint = midpoint, 
                      low = "lightblue", na.value = NA) +
  coord_equal() +
  labs(x = "Longtitude (°E)", y = "Latitude (°N)")

###MaxEnt model
ggplot() +
  geom_raster(data = predict.presence.me.df, aes(x,y, fill = layer2)) +
  geom_polygon(aes(long, lat, group = group), data = thai1, fill = NA, color = "gray80", size = 0.1) +
  geom_point(data = gbif.csv.xy, aes(lon, lat), color = "red", shape = 1) +
  scale_fill_gradient("probability", 
                      high = "forestgreen", low = "lightyellow", na.value = NA) +
  coord_equal() +
  theme_minimal() +
  labs(x = "Longtitude (°E)", y = "Latitude (°N)")

###GLM model
midpoint <- (max(predict.presence.glm.df$layer2, na.rm = TRUE)
             -min(predict.presence.glm.df$layer2, na.rm = TRUE))/2
ggplot() +
  geom_raster(data = predict.presence.glm.df, aes(x,y, fill = layer2)) +
  geom_polygon(aes(long, lat, group = group), data = thai1, fill = NA, color = "gray80", size = 0.1) +
  geom_point(data = gbif.csv.xy, aes(lon, lat), color = "red", shape = 1) +
  scale_fill_gradient2("probability", 
                       high = "red", mid = "yellow", midpoint = midpoint, 
                       low = "lightblue", na.value = NA) +
  coord_equal() +
  labs(x = "Longtitude (°E)", y = "Latitude (°N)")

###Random Forest model
midpoint <- (max(predict.presence.rf.df$layer2, na.rm = TRUE)
             -min(predict.presence.rf.df$layer2, na.rm = TRUE))/2
ggplot() +
  geom_raster(data = predict.presence.rf.df, aes(x,y, fill = layer2)) +
  geom_polygon(aes(long, lat, group = group), data = thai1, fill = NA, color = "gray80", size = 0.1) +
  geom_point(data = gbif.csv.xy, aes(lon, lat), color = "red", shape = 1) +
  scale_fill_gradient2("probability", 
                       high = "red", mid = "yellow", midpoint = midpoint, 
                       low = "lightblue", na.value = NA) +
  coord_equal() +
  labs(x = "Longtitude (°E)", y = "Latitude (°N)")

#3.5 Evaluate the model
##Before remove SSB
e.bc.0 <- evaluate(bc.model, p=gbif.csv.xy.test, a=bg.test, x= predictors)
e.bc.0
e.me.0 <- evaluate(me.model, p=gbif.csv.xy.test, a=bg.test, x= predictors)
e.me.0
e.glm.0 <- evaluate(glm.model, p=gbif.csv.xy.test, a=bg.test, x= predictors)
e.glm.0
e.rf.0 <- evaluate(rf.model, p=gbif.csv.xy.test, a=bg.test, x= predictors)
e.rf.0

##Visualize AUC with ROC curve
plot(e.bc.0, "ROC")
plot(e.me.0, "ROC")
plot(e.glm.0, "ROC")
plot(e.rf.0, "ROC")

##After remove SSB
e.bc <- evaluate(bc.model, p=gbif.csv.xy.test.pwd, a=bg.test.pwd, x= predictors)
e.bc
e.me <- evaluate(me.model, p=gbif.csv.xy.test.pwd, a=bg.test.pwd, x= predictors)
e.me
e.glm <- evaluate(glm.model, p=gbif.csv.xy.test.pwd, a=bg.test.pwd, x= predictors)
e.glm

#AUC = Area Under Curve is the area under the curve of True Positive Rate (TPR) vs. False Positive Rate (FPR)
#If AUC is high, it means more likely to get the correct prediction (true positives)

##Visualize AUC with ROC curve
plot(e.bc, "ROC")
plot(e.me, "ROC")
