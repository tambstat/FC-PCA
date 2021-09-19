# PCA using Frequency components #

# Ancillary Function Definitions # 


# 1st Ancillary function: Finding periodogram
pgram.gen = function(indxs1 , ss1 , n1) {
  return(ss1[,,indxs1])
}

# 2nd Ancillary function: Find periodogram times the kernel weight at one frequency
pgram_mat_freqs <- function(indxs2 , pgm2 , h2 , n2 , scan.frq2 ){
  dist2 <- scan.frq2 - indxs2
  pgm2[[floor(abs(indxs2)*(n2))]]*bp.kern((dist2)/h2 , h2)
}


# 3rd ancillary function: returns kernel spectral esimate at a one (uvec3) frequency
kspec_list_of_mats = function( uvec3 , pgm3 , h3 , n3, nu3 , aug.frq3  ){
  
  if(uvec3>0.5) uvec3 = 1 - uvec3 
  
  scan.frq3 <- uvec3
  loc.indx3 <- which.min(abs(aug.frq3 - scan.frq3))[1]
  select.frq3 <- aug.frq3[(loc.indx3 - n3/2+1):(loc.indx3 + n3/2)]
  
  pgms.sum3 <- list()
  pgms.sum3  <- lapply(select.frq3 , pgram_mat_freqs , pgm2=pgm3 , h2=h3 , 
                       n2=n3 , scan.frq2 = scan.frq3 )
  
  sum11 <- Reduce('+' , pgms.sum3) 
  
  return(sum11/n3)
} # end function kspec_list_of_mats 


# 4th ancilary function: returns real part of a spectral matrix at freq corresp to indxs4
re_kern_mat <- function(indxs4 , ksp4 , freqs) {
  
  summ = Re(ksp4[[indxs4]])

  return(summ)
}

# 5th ancillary function: returns squared coherence at freq corresp to indxs5
cross_spec_coh <- function(indxs5 , ksp5) Mod(ksp5[[indxs5]][1,2])^2/Re(ksp5[[indxs5]][1,1])/Re(ksp5[[indxs5]][2,2])



# 7th Ancillary function
sum.over.k <- function(indxs7 , ksp.boot7){
  Mod(ksp.boot7[[indxs7]][1,2])^2/Re(ksp.boot7[[indxs7]][1,1])/Re(ksp.boot7[[indxs7]][2,2])
}
