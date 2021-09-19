
# PCA using frequency components: simulation example  #

# Principal Component Analysis using Frequency Components of Multivariate Time Series, Computational Statistics & Data Analysis, 157, 107164, 2021
# Computational Statistics & Data Analysis #

###################################################################################################################################################

# Model 2 of the paper is chosen for this illustration #
# Note: Model 2 is where dimension p=6, no. of sub-groups m=3 #


# read PCA using frequency components functions
source("pca_freq_components_functions.R")

# Note: install needed R libraries listed in lines 6-15 of the above source functions file


set.seed(31)

# series length
n1 = 500

# dimension of the series
p = 6

h.input = n1^(-0.15) # bandwidth choice for the spectral matrix estimate

# randomly generate orthogonal matrix of size pxp
A = randortho(p)

# generate the multivariate time series (Model 1 of the paper)
yt = matrix(0 , n1 , p )


xtemp1 = arima.sim( model = list(ar=c(0.9), ma = c(0.8,-0.2) , sd  = 1 ) ,   n = n1+3) 
xtemp2 = arima.sim( model = list(ar=c(1.25,-0.75,0.3) , sd  = 1 )   , n = n1+5) 
xtemp3 = arima.sim( model = list(ma = c(1,-1,-0.8) ,  sd  = 1 )   , n = n1+2) 

yt[,1] = xtemp1[1:n1] 
yt[,2] = xtemp1[2:(n1+1)]
yt[,3] = xtemp1[3:(n1+2)]

yt[,4] = xtemp2[1:n1]
yt[,5] = xtemp2[2:(n1+1)]

yt[,6] = xtemp3[1:n1]


xt = t( A%*%t(yt)  ) # xt is the n x p multivariate time series 


# plot the acf of the observed series
ggAcf(xt , main = "")


# apply FC-PCA: PCA using frequency components
pcafr = freq_pca_components( xt , h.input )


# plot the acf of the transformed (segmented) series yt
ggAcf(pcafr$y.seg , main="")
# observe the 3 sub-groups from the ACF plot 