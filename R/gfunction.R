gfun <- function(pars_obj, intercept=F){
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  #if (length(gfun_index)!=1) stop("Something when wrong with the g function. There are either too many or none of them in the pars_obj.")
  out = pars_obj[[gfun_index]]$mat  %*% pars_obj[[gfun_index]]$pars
  if (intercept){
    fix_index = which(sapply(pars_obj, function(x)x$type)=="fix.eff")
    out <- out + pars_obj[[fix_index]]$pars[1]
  }
  out
}

dgfun <- function(pars_obj){
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  #if (length(gfun_index)!=1) stop("Something when wrong with the g function. There are either too many or none of them in the pars_obj.")
  pars_obj[[gfun_index]]$mat1 %*% pars_obj[[gfun_index]]$pars
}

gfun_inv <- function(W, pars_obj){
  gfun_index = which(sapply(pars_obj, function(x)x$type)=="gfun")
  #fix_index = which(sapply(pars_obj, function(x)x$type)=="fix.eff")
  
  ## To find vstar = gfun_inv(W), we compare W with W0=gfun(v0), with v0 a 1000-array in [0,1]
  ## v <- pars_obj[[gfun_index]]$v
  ## Instead of the prev line, need to recompute pars_obj[[gfun_index]]$mat, according to the original knots and a new v <- seq(0,1,len=1000)
  v0 <- seq(0,1,len=1000)
  bases  <- basis2_mpl(v0,pars_obj[[gfun_index]]$knots,order=pars_obj[[gfun_index]]$order,which=c(1,2,3), splines="isplines")
  mat=bases$Psi
  mat1=bases$psi
  mat2=bases$psi2
  vars=apply(mat,2,var)
  idrop=which.min(vars)
  mat=mat[,-idrop]
  mat1=mat1[,-idrop]
  mat2=mat2[,-idrop]
  R <- formR2(v0, mat2)
  pars_obj[[gfun_index]]$v <- v0
  pars_obj[[gfun_index]]$mat <- mat
  pars_obj[[gfun_index]]$mat1 <- mat1
  pars_obj[[gfun_index]]$mat2 <- mat2
  pars_obj[[gfun_index]]$R <- R
  W0 <- gfun(pars_obj, intercept = TRUE)
  ##
  ind <- findInterval(W, W0)
  v0=c(0,v0,1)
  vstar=(v0[ind+1]+v0[ind+2])/2
  return(vstar)
}



