# It's too computationally expensive to compute traffic jam duration and mileage for each user. So we apply machine learning algorithm to estimate traffic jam
# duration and mileage using mean and standard deviation of speed for each trip. The overall process constains of three steps:
# 1. Merge overall mean and sd from mean and sd of speed for each trip;
# 2. Cross-Validation for model selection;
# 3. Test the simple model.


rm(list=ls())
root.path <- 'F:/data/江西人保/speed_cal/m_dura_with_sd'
setwd(root.path)
filelist <- list.files(root.path, pattern="*.csv")
datalist <- lapply(filelist, function(x) read.csv(x, header=T, sep=','))
trip.statistics <- do.call("rbind", datalist)
trip.statistics <- trip.statistics[trip.statistics$m>0,]

## Valid trip filter.
trip.statistics <- trip.statistics[60*trip.statistics$m/(trip.statistics$dura*trip.statistics$speed.mean)<=1.5,]

## Merge overall sd from trip sd and mean using parallel computing in Windows.
sd.merge <- function(dura.sd){
  ###
  #函数用于计算多段行程速度标准差之间总的标准差
  #输入字段为包含dura,speed.mean,sd三个字段的DF
  ###
  names(dura.sd) <- c('dura','mean','sd')
  n <- nrow(dura.sd)
  if (n == 1){
    sd.cal <- dura.sd
  } else if (n == 2){
    sd.cal <- c(sum(dura.sd$dura, na.rm = T), 
                sum(dura.sd$dura*dura.sd$mean, na.rm = T)/sum(dura.sd$dura, na.rm = T), 
                sqrt((sum((dura.sd$dura-1/60)*dura.sd$sd^2, na.rm = T) + 
                        dura.sd[1,'dura']*dura.sd[2,'dura']*(dura.sd[1,'mean']-dura.sd[2,'mean'])^2/(dura.sd[1,'dura']+dura.sd[2,'dura']))/
                       (sum(dura.sd$dura, na.rm = T)-1/60)))
  } else if (n >= 3){
    n.1 <- floor(n/2)
    sd.1 <- sd.merge(dura.sd[1:n.1,])
    sd.2 <- sd.merge(dura.sd[(n.1+1):n,])
    sd.cal <- sd.merge(as.data.frame(rbind(sd.1,sd.2)))
  }
  return(sd.cal)
}

users <- as.data.frame(unique(trip.statistics$deviceid))

start.time <- proc.time()
library(parallel)
cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("trip.statistics","sd.merge"))
user.data.par <- parRapply(cl = cl, users, 
                           function(x) 
                             matrix(c(x[1], 
                                      colSums(trip.statistics[trip.statistics$deviceid==x[1],c('dura','m','dura.20','m.20')], na.rm = T), 
                                      sd.merge(trip.statistics[trip.statistics$deviceid==x[1],c('dura','speed.mean','speed.sd')])), 
                                    nrow = 1))						  
stopCluster(cl)
print(paste(round((proc.time()-start.time)[3],1),'seconds'))
user.data <- as.data.frame(matrix(unlist(user.data.par), ncol = 8, byrow = T), stringsAsFactors = F)
user.data <- user.data[,-6]
names(user.data) <- c('deviceid','dura','m','dura.20','m.20','speed.mean','speed.sd')

user.data <- user.data[user.data$m>100,]
user.data$m.ratio <- 100*user.data$m.20/user.data$m
user.data$dura.ratio <- 100*user.data$dura.20/user.data$dura
user.data$speed.20 <- 60*user.data$m.20/user.data$dura.20

user.data.train <- user.data[which(user.data$speed.20<20 & user.data$m.ratio<user.data$dura.ratio),]


## Cross-Validation for Model Selection
k <- 10
sample.size <- floor(nrow(user.data.train)/k)
sample.label <- matrix(sample(1:nrow(user.data.train))[1:(k*sample.size)], 
                       nrow = sample.size, ncol = k)
mse <- as.data.frame(matrix(0, nrow=10, ncol=24, dimnames=list(NULL,c(paste('fit.m',1:12,sep='.'),paste('fit.dura',1:12,sep='.')))))
for (i in 1:k){
  require(MASS)
  data.train <- user.data.train[-sample.label[,i],]
  data.test <- user.data.train[sample.label[,i],]
  
  fit.m.1 <- lm(data = data.train, m.ratio~speed.mean+speed.sd+I(speed.mean^2)+I(speed.sd^2))
  bc <- boxcox(fit.m.1); lambda <- bc$x[which.max(bc$y)]
  fit.m.1 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.mean+speed.sd+I(speed.mean^2)+I(speed.sd^2))
  mse[i,'fit.m.1'] <- sum(((predict(fit.m.1, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  fit.m.2 <- lm(data = data.train, m.ratio~speed.mean+speed.sd+I(speed.sd^2))
  bc <- boxcox(fit.m.2); lambda <- bc$x[which.max(bc$y)]
  fit.m.2 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.mean+speed.sd+I(speed.sd^2))
  mse[i,'fit.m.2'] <- sum(((predict(fit.m.2, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  fit.m.3 <- lm(data = data.train, m.ratio~speed.mean+speed.sd+I(speed.mean^2))
  bc <- boxcox(fit.m.3); lambda <- bc$x[which.max(bc$y)]
  fit.m.3 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.mean+speed.sd+I(speed.mean^2))
  mse[i,'fit.m.3'] <- sum(((predict(fit.m.3, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  fit.m.4 <- lm(data = data.train, m.ratio~speed.mean+speed.sd)
  bc <- boxcox(fit.m.4); lambda <- bc$x[which.max(bc$y)]
  fit.m.4 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.mean+speed.sd)
  mse[i,'fit.m.4'] <- sum(((predict(fit.m.4, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  fit.m.5 <- lm(data = data.train, m.ratio~speed.mean+I(speed.mean^2))
  bc <- boxcox(fit.m.5); lambda <- bc$x[which.max(bc$y)]
  fit.m.5 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.mean+I(speed.mean^2))
  mse[i,'fit.m.5'] <- sum(((predict(fit.m.5, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  fit.m.6 <- lm(data = data.train, m.ratio~speed.sd+I(speed.sd^2))
  bc <- boxcox(fit.m.6); lambda <- bc$x[which.max(bc$y)]
  fit.m.6 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.sd+I(speed.sd^2))
  mse[i,'fit.m.6'] <- sum(((predict(fit.m.6, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  fit.m.7 <- lm(data = data.train, m.ratio~speed.mean)
  bc <- boxcox(fit.m.7); lambda <- bc$x[which.max(bc$y)]
  fit.m.7 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.mean)
  mse[i,'fit.m.7'] <- sum(((predict(fit.m.7, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  fit.m.8 <- lm(data = data.train, m.ratio~speed.sd)
  bc <- boxcox(fit.m.8); lambda <- bc$x[which.max(bc$y)]
  fit.m.8 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.sd)
  mse[i,'fit.m.8'] <- sum(((predict(fit.m.8, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  fit.m.9 <- lm(data = data.train, m.ratio~speed.mean+I(speed.mean^2)+I(speed.sd^2))
  bc <- boxcox(fit.m.9); lambda <- bc$x[which.max(bc$y)]
  fit.m.9 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.mean+I(speed.mean^2)+I(speed.sd^2))
  mse[i,'fit.m.9'] <- sum(((predict(fit.m.9, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  fit.m.10 <- lm(data = data.train, m.ratio~speed.sd+I(speed.mean^2)+I(speed.sd^2))
  bc <- boxcox(fit.m.10); lambda <- bc$x[which.max(bc$y)]
  fit.m.10 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.sd+I(speed.mean^2)+I(speed.sd^2))
  mse[i,'fit.m.10'] <- sum(((predict(fit.m.10, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  fit.m.11 <- lm(data = data.train, m.ratio~speed.mean+I(speed.sd^2))
  bc <- boxcox(fit.m.11); lambda <- bc$x[which.max(bc$y)]
  fit.m.11 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.mean+I(speed.sd^2))
  mse[i,'fit.m.11'] <- sum(((predict(fit.m.11, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  fit.m.12 <- lm(data = data.train, m.ratio~speed.sd+I(speed.mean^2))
  bc <- boxcox(fit.m.12); lambda <- bc$x[which.max(bc$y)]
  fit.m.12 <- lm(data = data.train, (m.ratio^lambda-1)/lambda~speed.sd+I(speed.mean^2))
  mse[i,'fit.m.12'] <- sum(((predict(fit.m.12, newdata = data.test)*lambda+1)^(1/lambda)-data.test$m.ratio)^2)/nrow(data.test)
  
  
  
  
  fit.dura.1 <- lm(data = data.train, dura.ratio~speed.mean+speed.sd+I(speed.mean^2)+I(speed.sd^2))
  bc <- boxcox(fit.dura.1); lambda <- bc$x[which.max(bc$y)]
  fit.dura.1 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.mean+speed.sd+I(speed.mean^2)+I(speed.sd^2))
  mse[i,'fit.dura.1'] <- sum(((predict(fit.dura.1, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
  
  fit.dura.2 <- lm(data = data.train, dura.ratio~speed.mean+speed.sd+I(speed.sd^2))
  bc <- boxcox(fit.dura.2); lambda <- bc$x[which.max(bc$y)]
  fit.dura.2 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.mean+speed.sd+I(speed.sd^2))
  mse[i,'fit.dura.2'] <- sum(((predict(fit.dura.2, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
  
  fit.dura.3 <- lm(data = data.train, dura.ratio~speed.mean+speed.sd+I(speed.mean^2))
  bc <- boxcox(fit.dura.3); lambda <- bc$x[which.max(bc$y)]
  fit.dura.3 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.mean+speed.sd+I(speed.mean^2))
  mse[i,'fit.dura.3'] <- sum(((predict(fit.dura.3, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
  
  fit.dura.4 <- lm(data = data.train, dura.ratio~speed.mean+speed.sd)
  bc <- boxcox(fit.dura.4); lambda <- bc$x[which.max(bc$y)]
  fit.dura.4 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.mean+speed.sd)
  mse[i,'fit.dura.4'] <- sum(((predict(fit.dura.4, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
  
  fit.dura.5 <- lm(data = data.train, dura.ratio~speed.mean+I(speed.mean^2))
  bc <- boxcox(fit.dura.5); lambda <- bc$x[which.max(bc$y)]
  fit.dura.5 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.mean+I(speed.mean^2))
  mse[i,'fit.dura.5'] <- sum(((predict(fit.dura.5, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
  
  fit.dura.6 <- lm(data = data.train, dura.ratio~speed.sd+I(speed.sd^2))
  bc <- boxcox(fit.dura.6); lambda <- bc$x[which.max(bc$y)]
  fit.dura.6 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.sd+I(speed.sd^2))
  mse[i,'fit.dura.6'] <- sum(((predict(fit.dura.6, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
  
  fit.dura.7 <- lm(data = data.train, dura.ratio~speed.mean)
  bc <- boxcox(fit.dura.7); lambda <- bc$x[which.max(bc$y)]
  fit.dura.7 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.mean)
  mse[i,'fit.dura.7'] <- sum(((predict(fit.dura.7, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
  
  fit.dura.8 <- lm(data = data.train, dura.ratio~speed.sd)
  bc <- boxcox(fit.dura.8); lambda <- bc$x[which.max(bc$y)]
  fit.dura.8 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.sd)
  mse[i,'fit.dura.8'] <- sum(((predict(fit.dura.8, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
  
  fit.dura.9 <- lm(data = data.train, dura.ratio~speed.mean+I(speed.mean^2)+I(speed.sd^2))
  bc <- boxcox(fit.dura.9); lambda <- bc$x[which.max(bc$y)]
  fit.dura.9 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.mean+I(speed.mean^2)+I(speed.sd^2))
  mse[i,'fit.dura.9'] <- sum(((predict(fit.dura.9, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
  
  fit.dura.10 <- lm(data = data.train, dura.ratio~speed.sd+I(speed.mean^2)+I(speed.sd^2))
  bc <- boxcox(fit.dura.10); lambda <- bc$x[which.max(bc$y)]
  fit.dura.10 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.sd+I(speed.mean^2)+I(speed.sd^2))
  mse[i,'fit.dura.10'] <- sum(((predict(fit.dura.10, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
  
  fit.dura.11 <- lm(data = data.train, dura.ratio~speed.mean+I(speed.sd^2))
  bc <- boxcox(fit.dura.11); lambda <- bc$x[which.max(bc$y)]
  fit.dura.11 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.mean+I(speed.sd^2))
  mse[i,'fit.dura.11'] <- sum(((predict(fit.dura.11, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
  
  fit.dura.12 <- lm(data = data.train, dura.ratio~speed.sd+I(speed.mean^2))
  bc <- boxcox(fit.dura.12); lambda <- bc$x[which.max(bc$y)]
  fit.dura.12 <- lm(data = data.train, (dura.ratio^lambda-1)/lambda~speed.sd+I(speed.mean^2))
  mse[i,'fit.dura.12'] <- sum(((predict(fit.dura.12, newdata = data.test)*lambda+1)^(1/lambda)-data.test$dura.ratio)^2)/nrow(data.test)
}
mse.mean <- colMeans(mse)

#selected model:
#m.ratio~mean+sd+sd^2
#dura.ratio~mean+sd+mean^2+sd^2
fit.m <- lm(data = user.data.train, m.ratio~speed.mean+speed.sd+I(speed.sd^2))
bc <- boxcox(fit.m); lambda.m <- bc$x[which.max(bc$y)]
fit.m <- lm(data = user.data.train, (m.ratio^lambda.m-1)/lambda.m~speed.mean+speed.sd+I(speed.sd^2))

fit.dura <- lm(data = user.data.train, dura.ratio~speed.mean+speed.sd+I(speed.mean^2)+I(speed.sd^2))
bc <- boxcox(fit.dura); lambda.dura <- bc$x[which.max(bc$y)]
fit.dura <- lm(data = user.data.train, (dura.ratio^lambda.dura-1)/lambda.dura~speed.mean+speed.sd+I(speed.mean^2)+I(speed.sd^2))



##prepare the test data
root.path.test <- 'F:/data/江西人保/speed_cal/data_test'
setwd(root.path.test)
filelist.test <- list.files(root.path.test, pattern="*.csv")
datalist.test <- lapply(filelist.test, function(x) read.csv(x, header=T, sep=','))
test.data <- do.call("rbind", datalist.test)
test.data <- test.data[test.data$m>0,]
users.test <- as.data.frame(unique(test.data$deviceid))
start.time <- proc.time()
library(parallel)
cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("test.data","sd.merge"))
user.data.par <- parRapply(cl = cl, users.test, 
                           function(x) 
                             matrix(c(x[1], 
                                      colSums(test.data[test.data$deviceid==x[1],c('dura','m','dura.20','m.20')], na.rm = T), 
                                      sd.merge(test.data[test.data$deviceid==x[1],c('dura','speed.mean','speed.sd')])), 
                                    nrow = 1))						  
stopCluster(cl)
print(paste(round((proc.time()-start.time)[3],1),'seconds'))
user.test <- as.data.frame(matrix(unlist(user.data.par), ncol = 8, byrow = T), stringsAsFactors = F)
user.test <- user.test[,-6]
names(user.test) <- c('deviceid','dura','m','dura.20','m.20','speed.mean','speed.sd')
user.test <- user.test[user.test$m>100,]
user.test$m.ratio <- 100*user.test$m.20/user.test$m
user.test$dura.ratio <- 100*user.test$dura.20/user.test$dura
user.test$speed.20 <- 60*user.test$m.20/user.test$dura.20
user.test <- user.test[user.test$speed.20<20&user.test$m.ratio<=user.test$dura.ratio,]



##cal the mse of the model on test data
user.test$m.predict <- (predict(fit.m, newdata = user.test)*lambda.m+1)^(1/lambda.m)
user.test$dura.predict <- (predict(fit.dura, newdata = user.test)*lambda.dura+1)^(1/lambda.dura)
sum(((predict(fit.m, newdata = user.test)*lambda.m+1)^(1/lambda.m)-user.test$m.ratio)^2)/nrow(user.test)
sum(((predict(fit.dura, newdata = user.test)*lambda.dura+1)^(1/lambda.dura)-user.test$dura.ratio)^2)/nrow(user.test)


## Save the model and parameters.
save(fit.m,fit.dura,lambda.m,lambda.dura, file='F:/data/江西人保/speed_cal/fits_traffic.Rdata')




## Some plot of the data
setwd('../pictures/')
library(ggplot2)
ggplot(user.data.train, aes(m.ratio)) + geom_histogram(binwidth = 1, fill = 'steelblue', color = 'white')
ggsave(filename = '拥堵里程比例分布.png')
ggplot(user.data.train, aes(dura.ratio)) + geom_histogram(binwidth = 1, fill = 'steelblue', color = 'white')
ggsave(filename = '拥堵时长比例分布.png')

library(plotly)
p <- plot_ly(data = user.data.train, x = ~speed.mean, y = ~m.ratio, type = 'scatter', 
             mode = 'markers', name = '拥堵里程比例', alpha = 0.5) %>% 
  add_trace(y = ~dura.ratio, name = '拥堵时间比例', mode = 'markers', alpha = 0.5)%>% 
  layout(legend = list(x = 0.8, y = 0.9), xaxis = list(title = '平均速度'), yaxis = list(title = '比例'))
export(p, "拥堵时长&里程VS速度平均值.png")

p.1 <- plot_ly(data = user.data.train, x = ~speed.sd, y = ~m.ratio, type = 'scatter', 
               mode = 'markers', name = '拥堵里程比例', alpha = 0.5) %>% 
  add_trace(y = ~dura.ratio, name = '拥堵时间比例', mode = 'markers', alpha = 0.5)%>% 
  layout(legend = list(x = 0.8, y = 0.9), xaxis = list(title = '速度标准差'), yaxis = list(title = '比例'))
export(p.1, "拥堵时长&里程VS速度标准差.png")

user.data.train$m.predict <- (predict(fit.m, newdata = user.data.train)*lambda.m+1)^(1/lambda.m)
user.data.train$dura.predict <- (predict(fit.dura, newdata = user.data.train)*lambda.dura+1)^(1/lambda.dura)
p.3 <- plot_ly(data = user.data.train, x = ~dura.ratio, y = ~m.ratio, type = 'scatter', 
               mode = 'markers', name = '训练数据', alpha = 0.5) %>% 
  add_trace(x = ~dura.predict, y = ~m.predict, name = '预测结果', mode = 'markers', alpha = 0.5) %>% 
  layout(legend = list(x = 0.1, y = 0.9), xaxis = list(title = '时长比例'), yaxis = list(title = '里程比例'))
export(p.3, "训练集groundtruth&预测值.png")

p.4 <- plot_ly(data = user.test, x = ~dura.ratio, y = ~m.ratio, type = 'scatter', 
               mode = 'markers', name = '测试数据', alpha = 0.5) %>% 
  add_trace(x = ~dura.predict, y = ~m.predict, name = '预测结果', mode = 'markers', alpha = 0.5)%>% 
  layout(legend = list(x = 0.1, y = 0.9), xaxis = list(title = '时长比例'), yaxis = list(title = '里程比例'))
export(p.4, "测试集groundtruth&预测值.png")


 
p.8 <- plot_ly(data = user.data, alpha = 0.6) %>% 
  add_histogram(x = ~m.ratio, binwidth = 1, name = '拥堵里程比例') %>%
  add_histogram(x = ~dura.ratio, binwidth = 1, name = '拥堵时长比例') %>%
  layout(barmode = "overlay", legend = list(x = 0.8, y = 0.9), xaxis = list(title = '比例'))
export(p.8, "拥堵里程&时长分布.png")

