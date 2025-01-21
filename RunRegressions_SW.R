rm(list = ls())
set.seed(1)
source("functions.R")
d=5 #minimal ambient dimension of X_ij
dmax=20 #maximal ambient dimension of X_ij
n=100 #number of training Y_i
ntest=200 #number of test Y_i
N=200 #number of samples X_ij per distribution mu_i
alphainf = 0.05 #bounds for the distribution of dependence parameter alpha
alphasup = 0.1
betabar1 = 0.7 #bounds for the distributions of mean parameters beta1 and beta2
betabar2 = 0.7
multScore = 10 #multiplies X_ij1 - X_ij2 to compute P(Y_ij=1|X_ij)
nproj = 100 #number of random projections for sliced Wasserstein
nCV = 10 #number of cross validation splits between training and test
nTrainCV = 80 #size of training set at each cross validation split
vd = c(5,10,15,20) #vector of dimensions for which distribution regression is done 
nHypCV = 10 #number of values to be tested for scale parameter and ridge parameter
nmc = 50 #number of Monte Carlo repetitions to generate boxplots

#Multidimensional arrays of true values Y_i and predicted values Y_i:
TYtest_1 = array(dim=c(nmc,length(vd),ntest)) #for experiment on impact of ambient dimension d
ThatYtest_1 = array(dim=c(nmc,length(vd),ntest))
TYtest_2_Plus = array(dim=c(nmc,d,n)) #for experiment on effect of each variable k
TYtest_2_Minus = array(dim=c(nmc,d,n))
ThatYtest_2_Plus = array(dim=c(nmc,d,n))
ThatYtest_2_Minus = array(dim=c(nmc,d,n))
for (imc in 1:nmc) {  #Loop on Monte Carlo repetitions
  print(paste0("Monte Carlo run ",imc," out of ",nmc))
  generateData(d,dmax,n,ntest,N,alphainf,alphasup,betabar1,
               betabar2,multScore) #Generate inputs and outputs for distribution regression
  #They are stored in folder Data
  #Run distribution regressions and returns predicted values
  res = runRegression_SW(nproj,nCV,nTrainCV,vd,nHypCV) 
  #From res, get vectors of true and predicted Yi
  TYtest_1[imc,,] = res$MYtest_1
  ThatYtest_1[imc,,] = res$MhatYtest_1
  TYtest_2_Plus[imc,,] = res$MYtest_2_Plus
  TYtest_2_Minus[imc,,] = res$MYtest_2_Minus
  ThatYtest_2_Plus[imc,,] = res$MhatYtest_2_Plus
  ThatYtest_2_Minus[imc,,] = res$MhatYtest_2_Minus
}
#Save all results in Folder Data
save(TYtest_1,ThatYtest_1,TYtest_2_Plus,TYtest_2_Minus,
     ThatYtest_2_Plus,ThatYtest_2_Minus,vd = vd,nmc=nmc,file="Data/res_SW.rda")



