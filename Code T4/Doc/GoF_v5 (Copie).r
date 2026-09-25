##############################
##############################
# CHARGEMENT DES PACKAGES ####
##############################
##############################

#rm(list=ls())

library(doParallel)
library(foreach)
library(doRNG)

##############################
##############################
# PARALELISATION #############
##############################
##############################

NbCores <- detectCores()
cl <- makeCluster(NbCores)
registerDoParallel(cl)
getDoParWorkers()

##############################
##############################
# PARAMETRES #################
##############################
##############################

#load("tps.Rd")
#tps1= seq(from=250,by=250,to=4000)
#tps2= c(250,500,750,1000,1500,2000,2500,3000,4000)
#tps3= c(250,1250,1750,2250,2750,3000,3500,3750,4000)

#tps= tps1

conf.level= 0.05


a= 0.029
b= 0.071
#a= 8e-3
#b= 0.5



#set.seed(50)

# moyenne= a*b*tps
# variance= a*b^2*tps
# 
dtps= diff(c(0,tps))
# min.dtps= min(dtps)

m= length(dtps)

#n= 15

nbrep= 100000
# nbrep= 100

##############################
##############################
# ESTIMATION DE LA ZONE ######
##############################
##############################

funct.estimation <- function() {
 Y= NULL
 dtps= diff(c(0,tps))
 min.dtps= min(dtps)
 for (i in 1:n) {
  aux= rgamma(m,shape=a*dtps,scale=b)
  #aux= rnorm(m,mean=a*b*dtps,sd=sqrt(a*b^2*dtps))
  Y= cbind(Y,aux)
  Y[Y<0]= NA
 }
 Z= log(Y)
 mat.dtps= matrix(rep(dtps,n),ncol=n)
 mu.hat= sum(Y/mat.dtps,na.rm=TRUE)
 var.hat= mean((Y-mu.hat*mat.dtps)^2/mat.dtps,na.rm=TRUE)
 a.hat= mu.hat^2/var.hat
 b1.hat= var.hat/mu.hat
 b2.hat= mean((Y-a.hat*b1.hat*dtps)*(Z-digamma(a.hat*dtps)+log(b1.hat)),na.rm=TRUE)
 b2.hat= mean((Y-a.hat*b1.hat*dtps)*(Z-digamma(a.hat*dtps)-log(b1.hat)),na.rm=TRUE)
 stat1= b1.hat/b2.hat
 return(c(a.hat,b1.hat,b2.hat,stat1))
}

estimation <- foreach(w=1:nbrep,.combine=cbind,.options.RNG=123) %dorng% funct.estimation()
estimation = t(as.table(estimation))
a.hat <- estimation[,1]
b1.hat <- estimation[,2]
b2.hat <- estimation[,3]
vstat1 <- estimation[,4]

# pdf("kde_h0.pdf")
# moyenne= mean(vstat1)
# ecart.type= sd(vstat1)
# norm.vstat1 = (vstat1-moyenne)/ecart.type
# norm.param= c(moyenne,ecart.type)
# save(x=norm.param,file="norm.param.Rd")
# plot(density(norm.vstat1),lwd=2)
# x= seq(from=-5,by=0.01,to=5)
# y= dnorm(x)
# lines(x,y,lwd=2,lty=3)
# dev.off()
# 
# library(nortest)
# lillie.test(norm.vstat1)
# 
# i.vstat1 <- 1/vstat1
# plot(density(i.vstat1),lwd=2)
# 
# moyenne= mean(i.vstat1)
# ecart.type= sd(i.vstat1)
# norm.i.vstat1 = (i.vstat1-moyenne)/ecart.type
# 
# plot(density(norm.i.vstat1),lwd=2)
# x= seq(from=-2,by=0.01,to=6)
# y= dnorm(x)
# lines(x,y,lwd=2,lty=3)
# lillie.test(norm.i.vstat1)
# 
# l.vstat1 <- log(vstat1)
# norm.l.vstat1 <- scale(l.vstat1)
# plot(density(norm.l.vstat1),lwd=2)
# lillie.test(norm.l.vstat1)
# 
# l.i.vstat1 <- log(i.vstat1)
# norm.l.i.vstat1 <- scale(l.i.vstat1)
# plot(density(norm.l.i.vstat1),lwd=2)
# x= seq(from=-2,by=0.01,to=6)
# y= dnorm(x)
# lines(x,y,lwd=2,lty=3)
# lillie.test(norm.l.i.vstat1)


limits= quantile(vstat1,probs=c(conf.level/2,1-conf.level/2))
#save(x=limits,file="limits.Rd")
cat("\n\n")
cat("*** ZONE D'ACCEPTATION : ***\n")
cat("Limite inférieure :",limits[1],"\n")
cat("Limite supérieure :",limits[2],"\n\n")


##############################
##############################
# CALCUL DE PUISSANCE ########
# PROCESSUS DE GNH    ########
##############################
##############################

a= 0.029/500
b= 0.071
beta <- 1.8

funct.estimation.nhgp <- function() {
  Y= NULL
  dtps= diff(c(0,tps^beta))
  min.dtps= min(dtps)
  for (i in 1:n) {
    aux= rgamma(m,shape=a*dtps,scale=b)
    Y= cbind(Y,aux)
    Y[Y<0]= NA
  }
  Z= log(Y)
  mat.dtps= matrix(rep(dtps,n),ncol=n)
  mu.hat= sum(Y/mat.dtps,na.rm=TRUE)
  var.hat= mean((Y-mu.hat*mat.dtps)^2/mat.dtps,na.rm=TRUE)
  a.hat= mu.hat^2/var.hat
  b1.hat= var.hat/mu.hat
  b2.hat= mean((Y-a.hat*b1.hat*dtps)*(Z-digamma(a.hat*dtps)+log(b1.hat)),na.rm=TRUE)
  b2.hat= mean((Y-a.hat*b1.hat*dtps)*(Z-digamma(a.hat*dtps)-log(b1.hat)),na.rm=TRUE)
  stat1= b1.hat/b2.hat
  return(c(a.hat,b1.hat,b2.hat,stat1))
}

estimation.nhgp <- foreach(w=1:nbrep,.combine=cbind,.options.RNG=123) %dorng% funct.estimation.nhgp()
estimation.nhgp = t(as.table(estimation.nhgp))
a.hat <- estimation.nhgp[,1]
b1.hat <- estimation.nhgp[,2]
b2.hat <- estimation.nhgp[,3]
vstat1 <- estimation.nhgp[,4]

cnt= which((vstat1<limits[1])|(vstat1>limits[2]))
pw.nhgp= length(cnt)/length(vstat1)

norm.vstat1 = (vstat1-norm.param[1])/norm.param[2]
#plot(density(norm.vstat1),lwd=2)
#x= seq(from=-5,by=0.01,to=5)
#y= dnorm(x)
#lines(x,y,lwd=2,lty=3)

cat("*** CALCULS DE PUISSANCE : ***\n")
cat("Processus de Gamma non-homogène :",pw.nhgp,"\n")

cat("\n\n")




##############################
##############################
# ARRET PARALELLISATION ######
##############################
##############################

stopCluster(cl)

