# PCA using frequency components of multivariate time series: FC - PCA #

######################################################################################################
# LIBRARIES

library(boot)
library(np)
library(vars)
library(astsa)
library(igraph)
library(Matrix)
library(cmvnorm)
library(pracma)
library(ggplot2)
library(forecast)

####################################################################################################

# Additional ancillary functions
source("ancillary_functions.R")


# Bartlett-Priestley Kernel  K(.) #
bp.kern = function(locc,hin){
  
  if(abs(locc*hin) <= hin ) kval = 3/4*(1-locc^2)/hin
  if(abs(locc*hin)> hin ) kval = 0
  
  return(kval)
  
} # end function bp.kern


# Kernel spectral estimator #
kernel_spec = function( x , uvec , h ){
  
  p = dim(x)[2]
  n = dim(x)[1] 
  
  nu = length(uvec)

  # Periodogram
  ss <- mvspec(x,plot=FALSE)$fxx
  ss <- ss/2/pi
  pgm <- list() 
  
  # periodogram list of size n
  #pgm <- lapply( 1:n , pgram.gen , ss1 = ss , n1 = n )
  pgm <- lapply( 1:(n/2) , pgram.gen , ss1 = ss , n1 = n )
  pgm <- c( pgm , rev(pgm) )
  
  aug.frq = c( -rev((1:n)/n) , (1:n)/n , rev( (1:n)/n )  )
  
  kspec <- list()
  
  kspec <- lapply(uvec , kspec_list_of_mats , pgm3 = pgm , h3=h , n3=n , 
                  nu3=nu , aug.frq3=aug.frq )
  
  return(kspec)
  
} # end function kernel_spec


# Inputs a list of periodograms and returns kernel spectral estimate list
kernel_spec_using_pgram = function(pgm.in , uvec , h.in){
  
  p = dim(pgm.in[[1]])[1]  #dim(pgm.in)[1]
  n = length(pgm.in) #dim(pgm.in)[3] 
  
  nu = length(uvec)

  aug.frq <- c( -rev((1:n)/n) , (1:n)/n , rev( (1:n)/n )  )
  
  kspec.pgm.smooth <- lapply(uvec , kspec_list_of_mats , pgm3=pgm.in , h3=h.in , 
                             n3=n , nu3=nu , aug.frq3=aug.frq )
  
  return(kspec.pgm.smooth)
} # end function kernel_spec_using_pgram


# Eigendecomposition of sum of real part of spectral matrices
pca_freq_eigenvectors <- function( x , frq.vec , h ){
  
  p = dim(x)[2]
  n = dim(x)[1]
  
  n.frq = length(frq.vec)
  
  kern.sp <- kernel_spec(x,frq.vec,h)
  
  sum.spec.mat <- list() #matrix(0,p,p)
  
  rknmats <- lapply(1:n.frq , re_kern_mat , ksp4 = kern.sp ,freqs = frq.vec )
  sum.spec.mat <- Reduce('+' , rknmats)
  
  ega = eigen(sum.spec.mat)
  evec = ega$vectors
  
  return(list(evec, t(t(evec)%*%t(x))  ))
  
} # end function pca_freq_eigenvectors


# function to evaluate squared cross spectrum of a bivariate series z
sq_cross_spectrum = function(z,h.in,uvec){
  
  n = dim(z)[1]
  p = 2
  
  nu = length(uvec)
  
  ksp.z = kernel_spec(z,uvec,h.in)
  
  sums <- list()
  sums <- lapply( 1:nu , cross_spec_coh , ksp5 = ksp.z  )
  sum1 <- sum(unlist(sums))
  
  return( sum1*sqrt(h.in)*n/nu) 

} # end function sq_cross_spectrum

 

# function to find the group of components
find_groups = function(emat){
  
  p = dim(emat)[1]
  
  grps = list()
  size.grps = NULL
  g=1
  
  grph = graph_from_adjacency_matrix(emat)
  comps = components(grph)
  n.grp = comps$no
  
  g = 1
  for(i in 1:n.grp){
    grps[[g]] = which(comps$membership==i)
    size.grps = c(size.grps , length(grps[[g]]))
    g = g + 1
  }
  
  cc.out = cbind(comps$membership , 1:p)

  # Naming the largest group as group no. 1
  grp.mat = cbind(1:n.grp,size.grps)
  grp.mat = grp.mat[ order(grp.mat[,2] , decreasing = T) , ]
  
  if(n.grp>1){
  for(r in 1:n.grp){
    group.no = grp.mat[r,1]
    indx = which(comps$membership==group.no)
    cc.out[ indx , 1 ] = r
  }}

  size.grps = NULL
  for(k in 1:n.grp) size.grps = c(size.grps , length(which(cc.out[,1]==k)))
  
  return(list(cc.out , n.grp , size.grps)) 
} # end function find_groups


# testing for zero cross-spectrum among every pair of components
test_pairwise_components = function(x , frq.vec , siglevel , h){
  
  p = dim(x)[2]
  n = dim(x)[1]
  
  n.frq = length(frq.vec)
  
  pcs = pca_freq_eigenvectors(x,frq.vec,h)
  
  y = pcs[[2]] # transformed series 
  L = pcs[[1]] # eigenvectors
  
  e.mat = matrix(0,p,p) # adjacency matrix
  pval.mat = matrix(0,p,p) # matrix of p-values
  
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      
      ytemp = cbind( y[,i] , y[,j]  )
      
      tstat.actual = sq_cross_spectrum(ytemp , h.in = h,frq.vec)#cross.spec.info[[1]]
      
      # Critical values based on Bartlett-Priestley kernel
      mu.tstat = 6/5
      sig.tstat = 2*pi*2672*pi/385
      
      pval = 1 - pnorm(tstat.actual , mean = mu.tstat/sqrt(h) , sd  = sqrt(sig.tstat) )

      if(pval<siglevel) e.mat[i,j] = e.mat[j,i] = 1
      if(pval>siglevel) e.mat[i,j] = e.mat[j,i] = 0
      
    } # end loop over i
  } # end loop over j 
  
  return(list(y,L,e.mat,pval.mat))
  
} # end function test_pairwise_components

# Wrapper function
freq_pca_components = function(x , h, frq.vec =  seq(0.05,0.95,0.1)
                               , siglevel = 0.01 ){
  
  p = dim(x)[2]
  n = dim(x)[1]
  
  pca_info <- test_pairwise_components(x,frq.vec,siglevel,h) 
  A.pca = pca_info[[2]] # estimated A matrix (L matrix)
  
  cc <- find_groups(pca_info[[3]])
  n.clusters = cc[[2]] # no. of groups/communities
  
  cc.mat = cc[[1]]
  
  cc.mat <- cc.mat[ order(cc.mat[,1]) , ]

  pmat <- t( matrix( as.numeric(  as( cc.mat[,2] , "pMatrix" ) )  , p , p ) ) 
  #print(cc[[3]])
  y.final <- t(pmat%*%t(pca_info[[1]]))

  output.list <- list(y.final , pca_info[[1]] , pca_info[[2]]  , pmat , 
                     cc[[1]][,1] , cc[[3]] , 
                     n.clusters , 
                     pca_info[[4]] )
  
  names(output.list) <- c( "y.seg" , "y.pre" , "L" , "P" , "group_number" 
                          , "sizes" , "clusters" , "p-value-matrix" )
  
  return(output.list)
  
} # end function freq_pca_components


# find subspace estimation error
d.tilde.error = function(A1,B1){
   
  r1 <- rankMatrix(A1)[1]
  r2 <- rankMatrix(B1)[1]

  P <- A1%*%solve( t(A1)%*%A1 )%*%t(A1)
  Q <- B1%*%solve( t(B1)%*%B1 )%*%t(B1)
  
return( 1 - 1/min(r1,r2)*sum(diag( P%*%Q  ))  )
} # end function d.tilde.error


# output 'correct segmentation' error from Chang et. al (2018)
subspace_error_count = function(x,A1,q,h){
  
  q.count  = 0
  corr.seg.count = 0
  
  d.vec = NULL
  
  p = dim(x)[2]
  n = dim(x)[1]
  
  freq.pca.info <- freq_pca_components(x,h)
  
  A.es = freq.pca.info[[3]]%*%freq.pca.info[[4]] # estimated A, after permuting columns
    
  q.es <- freq.pca.info$clusters # no. of clusters
  c.size <- unname(freq.pca.info$sizes) # cluster sizes 
  
  s.col = rep(1,q.es) # starting column index
  e.col = rep(c.size[1],q.es) # ending column index
  for(j in 2:q.es){
  s.col[j] = sum(c.size[1:(j-1)]) + 1
  e.col[j] = sum(c.size[1:j])
  }
  #if(q==q.es) q.count  = 1
  if(q.es==q & c.size[1]==3  & c.size[2]==2 & c.size[3]==1 ){
  
  q.count = 1
  corr.seg.count = 0
  for(j in 1:q.es){
      
      if(length(s.col[j]:e.col[j])>1) tmp = d.tilde.error(A.es[,s.col[j]:e.col[j]] , A1[,s.col[j]:e.col[j]])
      if(length(s.col[j]:e.col[j])==1) {
        tmp = d.tilde.error(matrix(A.es[,s.col[j]:e.col[j]]) , matrix(A1[,s.col[j]:e.col[j]]))
      }
      
      d.vec = c(d.vec , tmp)
      
      
      tmp.min = tmp
      
      incorr.seg.flag = 0
      
      for(k in  (1:q.es)[-j] ){
        
        if(c.size[j]!=1 & c.size[k]!=1) tmp.inner = d.tilde.error(A.es[,s.col[j]:e.col[j]] , A1[,s.col[k]:e.col[k]])
        if(c.size[j]==1 & c.size[k]!=1) tmp.inner = d.tilde.error(matrix(A.es[,s.col[j]:e.col[j]]) , A1[,s.col[k]:e.col[k]])
        if(c.size[j]!=1 & c.size[k]==1) tmp.inner = d.tilde.error(A.es[,s.col[j]:e.col[j]] , matrix(A1[,s.col[k]:e.col[k]]))
        if(c.size[j]==1 & c.size[k]==1) tmp.inner = d.tilde.error(matrix(A.es[,s.col[j]:e.col[j]]) , matrix(A1[,s.col[k]:e.col[k]]))
        
        
        if(tmp.inner < tmp)  incorr.seg.flag = 1  

      } # end loop over k 
      if(incorr.seg.flag==0)  corr.seg.count = corr.seg.count + 1
    
    } # end loop over j
  
  if(corr.seg.count==q) corr.seg.count = 1
  
} # end big IF loop
  
  return(c(q.count,corr.seg.count,max(d.vec) , mean(d.vec)))
} # end function subspace error count
      


# function to find eigengaps between an inputted list of spectral matrices
eigen_gap = function(spec.list,nfrq){
  
  ns = length(spec.list)/nfrq
  start.indx = seq(1,length(spec.list),nfrq)
  
  eigen.gp = 100 # initialize 
  for(i in 1:(ns-1)){
    for(j in (i+1):ns){
      
      p1 = dim(spec.list[[ start.indx[i] ]])[1]
      sum1 = matrix(0,p1,p1)
      for(k in 1:nfrq) sum1 = sum1 + Re( spec.list[[ start.indx[i] +k-1 ]] )
      t1 = eigen(sum1)$values
      
      p2 = dim(spec.list[[ start.indx[j] ]])[1]
      
      if( sum(is.null(p2))==1 ) p2=1
      
      sum2 = matrix(0,p2,p2)
      for(k in 1:nfrq) sum2 = sum2 + Re( spec.list[[ start.indx[j] + k -1   ]] )
      t2 = eigen(sum2)$values
      
      X = meshgrid(t1,t2)$X
      Y = meshgrid(t1,t2)$Y
      tmp.eigen.gap  = (min(min(abs(X-Y))))
      if(tmp.eigen.gap<eigen.gp)  eigen.gp = tmp.eigen.gap
    }
  }
  return(eigen.gp)
} # end function eigen_gap
