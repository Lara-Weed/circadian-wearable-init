### Sep 08, 2010, updates by Jia
### R code for FPCA
### 

source("func_fpca.r")

### read in activity data and clinical data
### 
activity=read.csv("activity.csv",header=TRUE)
ahi=read.csv("ahi.csv",header=TRUE)

### reformat colnames in activity data
colnames(activity)=sub("X","",colnames(activity)) 

### prepare data for smooth 
###
data=vector("list",2)
data[[1]]=as.matrix(activity)
data[[2]]=ahi

names(data)=c("mat","cov")

### Smooth activity with 9 Fourier bases functions
### smooth family can be Fourier or bspline
### fda.smoothdata returns a list of smoothed functional data and clinical data
### 
FD=fda.smoothdata(data,nbasis=9,basistype="Fourier") 

### FPCA for the smoothed result 
### smooth family can be Fourier or bspline
### fda.fpca returns a list of fpca result and clinical data
### more information in memo
### 
FPCA=fda.fpca(FD,nharm=4,nbasis=9,basistype="bspline",norder=4)

### plot FPCA result
### mean activity is same in all the four components plot
### "+" shows add the component to the mean
### "-" shows subtrac the component to the mean
###
par(mfrow=c(2,2))
plot(FPCA$fpca,cex.main=.8)

### colored plot based on PHQ9
cov=FPCA$cov
fpca=FPCA$fpca

### color code for plot
npt=length(cov) # number of patients
color=rep(0,npt)
color[cov<5]=2
color[cov>30]=4

harms = fpca$harmonics # harmonics
meanfd = fpca$meanfd # mean activity
L=meanfd$basis$rangeval[2] # get the number of samples
score=fpca$scores # scores matrix, row for each patients and column for each component

par(mfrow=c(1,1))
j=1### look at the 1st component

plot(0,0,xlim=c(0,L),ylim=c(0,800),xlab='Time',ylab='Activity',type='n',main=paste(j,"component",sep=" "),xaxt="n")      
for (i in 1:npt)
      lines(score[i,j]*harms[j]+meanfd,col=color[i],lwd=1,xaxt="n") 

axis(1,at=seq(0,1441,360),labels=c("Midnight","6AM","Noon","6PM","Midnight"))
legend("topleft",c("AHI normal", "AHI severe"),col=c(2,4),lty=1)


