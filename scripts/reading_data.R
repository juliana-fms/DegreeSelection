###-----------------------------------------------------###
###     File #1: reading and treating data              ###
###     Author: Juliana Freitas de Mello e Silva        ###
###-----------------------------------------------------###

### Removing objects
rm(list = ls())

### Loading necessary packages
require("fda")
data("growth")
head(growth)

### Acceleration of growth - boys
acceleration_growth_boys<- sapply(1:ncol(growth$hgtm),
                                  function(b)(as.numeric(diff(growth$hgtm[,b]))))
acceleration_growth_boys[which(acceleration_growth_boys < 0)]<- 0

### Acceleration of growth - girls
acceleration_growth_girls<- sapply(1:ncol(growth$hgtf),
                                  function(g)(as.numeric(diff(growth$hgtf[,g]))))
acceleration_growth_girls[which(acceleration_growth_girls < 0)]<- 0

### Data in the wide format
data_wide<- as.data.frame(rbind(t(acceleration_growth_boys), t(acceleration_growth_girls)))
data_wide$id<- 1:(ncol(acceleration_growth_boys) + ncol(acceleration_growth_girls))
data_wide$gender<- c(rep(x = 0, times = ncol(acceleration_growth_boys)),
                rep(x = 1, times = ncol(acceleration_growth_girls)))
colnames(data_wide)<- c(paste0("height", 1:(length(growth$age) - 1)), "id", "gender")

### Data in the long format
data_long<- data.frame(c(as.numeric(acceleration_growth_boys), as.numeric(acceleration_growth_girls)))
colnames(data_long)<- "height"
data_long$id<- as.numeric(sapply(1:(ncol(acceleration_growth_boys) + ncol(acceleration_growth_girls)),
                      function(i)(rep(x = i, times = length(growth$age) - 1))))
data_long$age<- rep(x = growth$age[-1],
                    times = ncol(acceleration_growth_boys) + ncol(acceleration_growth_girls))
data_long$gender<- c(rep(x = 0, times = ncol(acceleration_growth_boys)*length(growth$age[-1])),
                     rep(x = 1, times = ncol(acceleration_growth_girls)*length(growth$age[-1])))
data_long<- data_long[,c("id", "age", "gender", "height")]

### Saving data
write.table(x = data_long, file = "accelerated_growth_long.txt")


### Now you can go to file #2