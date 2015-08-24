library("R.matlab")
library("bvarr")
all <- readMat("../data/usa_data.mat")
str(all)

usa_data <- all$usa.data
usa_data
str(usa_data)
#colnames(usa_data)<-c('var1','var2','var3','var4','var5','var6','var7','var8','var9','var10','var11','var12','var13')
#data_with_names=data.frame(colnames usa_data)
#write_csv(data_with_names, "../data/usa_2015_final.csv")
#var1<-usa_data[:,1]
#as.data.frame(usa_data)
data_with_names <- data.frame(var1=usa_data[,1], var2=usa_data[,2], var3=usa_data[,3],var4=usa_data[,4],
                              var5=usa_data[,5], var6=usa_data[,6], var7=usa_data[,7],var8=usa_data[,8],
                              var9=usa_data[,9], var10=usa_data[,10], var11=usa_data[,11],var12=usa_data[,12], var13=usa_data[,13])
write_csv(data_with_names, "../data/usa_carriero.csv")
