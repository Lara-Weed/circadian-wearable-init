library(fda)


fda.smoothdata=function(data,basistype="Fourier",nbasis=9,norder=4) {
  mat=data$mat
  cov=data$cov
  L=dim(mat)[1] # get the number of data samples
  
  # Create basis functions
  if(basistype=="Fourier"){
    fbase = create.fourier.basis(rangeval=c(0,L), nbasis)
  } else if(basistype=="bspline") {
    fbase = create.bspline.basis(rangeval=c(0,L), nbasis, norder)
  } else {
    stop("Basis type not specified correctly!")
  }
  ### smooth the data
  fpar = fdPar(fbase) 
  fd=smooth.basis(c(1:L), mat, fpar)
  FD=list(fd=fd,cov=cov)
  return(FD)
}



## Functional principle component analysis

fda.fpca=function(FD,basistype="Fourier",nbasis=9,norder=4,nharm=4) {
  cov=FD$cov[[2]] # get covariates
  fd=FD$fd   # obtain smoothed data 
  L=dim(fd$y)[1] # get the number of data samples.

  # Create basis functions
  if(basistype=="Fourier") {
    fbase = create.fourier.basis(rangeval=c(0,L), nbasis)
  } else if(basistype=="bspline") {
    fbase = create.bspline.basis(rangeval=c(0,L), nbasis, norder)
  } else {
    stop("Basis type not specified correctly!")
  }
  fpar = fdPar(fbase) 
   
  fpca=pca.fd(fd$fd,nharm,fpar) # functional principle component analysis
  FPCA=list(fpca=fpca,cov=cov)
  return(FPCA)
}




