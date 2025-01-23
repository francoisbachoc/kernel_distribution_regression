MdistSW = function(X1,X2,Mproj) {
  #Compute the n1 times n2 matrix of sliced Wasserstein distances
  #between two arrays X1 and X2
  #X1 is an array of dimension [n1,N,d] 
  #corresponding to n1 distributions each represented by N d-dimensional vectors
  #Mproj is the matrix of d-dimensional randomly generated directions
  n1 = dim(X1)[1]
  N = dim(X1)[2]
  d = dim(X1)[3]
  n2 = dim(X2)[1]
  nproj = dim(Mproj)[1]
  TprojSort1 = array(dim=c(n1,N,nproj))
  for (i in 1:n1) {
    mi = Mproj %*% t(X1[i,,]) #projection on the random directions
    mis = apply(mi,1,sort) #sort of the N projected values
    TprojSort1[i,,] = mis
  }
  TprojSort2 = array(dim=c(n2,N,nproj))
  for (j in 1:n2) {
    mj = Mproj %*% t(X2[j,,]) #projection on the random directions
    mjs = apply(mj,1,sort) #sort of the N projected values
    TprojSort2[j,,] = mjs
  }
  Mdistt = matrix(nrow=n1,ncol=n2) 
  for (i in 1:n1) {
    for (j in 1:n2) {
      Mdistt[i,j] = mean((TprojSort1[i,,]-TprojSort2[j,,])^2) #average of square difference 
      #over directions and discrete probabilities
    }
  }
  Mdistt
}

MprojSW = function(nproj,d) {
  #Compute a matrix nproj times d of nproj uniformly-distributed directions in dimension d
  Mproj = matrix(nrow=nproj,ncol=d)
  for (a in 1:nproj){
    z = rnorm(d)
    Mproj[a,] = z / sqrt(sum(z^2))
  }
  Mproj
}

hatY = function(lc,lambda,MdistTrain,MdistTrainTest,Ytrain) {
  #Compute the vector of predictions by kernel distribution regression
  #with squared exponential kernel with scaling lc
  #and regularization lambda
  #with matrices of Hilbertian distance MdistTrain for the training data
  #and MdistTrainTest between training and test data
  #and with vector of training data Ytrain
  n = dim(MdistTrain)[1]
  Rtrain = exp(- MdistTrain^2 / lc^2)
  RtrainTest = exp(- MdistTrainTest^2 / lc^2)
  hatYtest =  t(RtrainTest)%*%solve(Rtrain+n*lambda*diag(n))%*%matrix(nrow=n,ncol=1,data=Ytrain)
  hatYtest
}


SelCV = function(vlc,vlambda,Mdist,nCV,nTrainCV,Ytrain) {
  #Select scaling parameter lc among a vector vlc
  #and regularization parameter lambda among a vector vlambda
  #for the matrix Mdist of Hilbertian distances for the training data
  #Selection is done by minimizing the sum of square error with cross validation
  #with nCV splits between a training subsample of size nTrainCV 
  #and a test subsample with the remaining data
  #The scalar predictands are in the vector Ytrain
  #Returns a list with the optimal lc and lambda and the 
  #matrix of sum of square errors
  n = dim(Mdist)[1]
  MMSE = matrix(nrow = length(vlc),ncol = length(vlambda))
  for (ilc in 1:length(vlc)) {    #double loop on lc and lambda
    for (ilambda in 1:length(vlambda)) {
      lc = vlc[ilc]
      lambda = vlambda[ilambda]
      SMSE = 0
      for (k in 1:nCV) {    #loop on the training and test subsamples
        perm = sample(1:n)
        MdistTrain = Mdist[perm[1:nTrainCV],perm[1:nTrainCV]]    #sub matrix of training Hilbertian distances
        MdistTrainTest = Mdist[perm[1:nTrainCV],perm[(nTrainCV+1):n]] #sub matrix of train-test Hilbertian distances
        YtrainCV = Ytrain[perm[1:nTrainCV]] #subvector of training predictands
        YtestCV = Ytrain[perm[(nTrainCV+1):n]] #subvector of test predictands
        hatYCV =  hatY(lc,lambda,MdistTrain,MdistTrainTest,YtrainCV) #kernel regression prediction  
        SMSE = SMSE + sum((YtestCV-hatYCV)^2)  #sum of square error computation
      }
      MMSE[ilc,ilambda] = SMSE
    }
  }
  ilc = which(MMSE == min(MMSE), arr.ind = TRUE)[1]
  lc = vlc[ilc]   #optimal lc
  ilambda = which(MMSE == min(MMSE), arr.ind = TRUE)[2]
  lambda = vlambda[ilambda]  #optimal lambda
  res = list(MMSE=MMSE,lc=lc,lambda=lambda)
  return(res)
}

generateData = function(d,dmax,n,ntest,N,alphainf,alphasup,betabar1,
                        betabar2,multScore) {
  #Generate data (X_ij,Y_i) that are saved in file Data as .rda files
  #d, dmax: smallest and largest ambient dimension
  #n, ntest: number of training and test observations Y_i
  #N: number of samples per Y_i
  #alphainf, alphasup: bounds for the distribution of dependence parameter alpha
  #betabar1, betabar2: bounds for the distributions of mean parameters beta1 and beta2
  #multScore: multiplies X_ij1 - X_ij2 to compute P(Y_ij=1|X_ij)
  ntot = n+ntest
  #generation of ntot parameters for distributions mu_i of samples X_ij
  valpha = runif(n=ntot,min=alphainf,max=alphasup) 
  vbeta1 = runif(n=ntot,min=-betabar1,max=betabar1)
  vbeta2 = runif(n=ntot,min=-betabar2,max=betabar2)
  #
  #Samples stored as arrays with dimensions ntot times N times dmax
  #since X_ij is at most a dmax-dimensional vector
  X = array(dim=c(ntot,N,dmax)) 
  mY = matrix(nrow=ntot,ncol=N) #matrix of the Y_ij
  for (i in 1:(ntot)) { #loop on i in X_ij
    alpha = valpha[i] #parameters of mu_i
    beta1 = vbeta1[i]
    beta2 = vbeta2[i]
    vmean = c(beta1,beta2,rep(0,dmax-2))  
    for (j in 1:N) { #loop on j in X_ij
      X[i,j,] = runif(n=1,min=-alpha,max=alpha)*rep(1,dmax) +
        rnorm(n=dmax,mean=vmean,sd=1/2) #sampling X_ij in dimension dmax
      score = multScore*(X[i,j,1]-X[i,j,2])
      p = exp(score)/(1+exp(score))  #P(Y_ij=1|X_ij)
      mY[i,j] = as.integer(runif(n=1) < p)
    }
  }
  Y = apply(mY,1,mean) #Computing Y_i as average over j of Y_ij
  #Selection and storage of the train and test outputs Y_i:
  Ytrain = Y[1:n] 
  save(Ytrain=Ytrain,file="Data/1_OutputTrain.rda")
  Ytest = Y[(n+1):ntot] 
  save(Ytest=Ytest,file="Data/1_OutputTest.rda")
  #Loop to create distribution regression problems in increasing dimension:
  for (dd in c(5,10,15,20)) {  
    #Selection and storage of the train and test inputs X_ij:
    Xtrain = X[1:n,,1:dd]
    save(Xtrain,file=paste0("Data/1_InputTrain_d",dd,".rda"))
    Xtest = X[(n+1):ntot,,1:dd]
    save(Xtest,file=paste0("Data/1_InputTest_d",dd,".rda"))
  }
  #For k=1...d, we generate data (X_ij,Y_i) to study 
  #impact of variable k
  #Hence we need an array with 4 entries k,i,j,l where l is the component of X_ij
  XmodPlus = array(dim=c(d,n,N,d)) #Plus: larger values of X_ijk
  XmodMinus = array(dim=c(d,n,N,d)) #Minus: larger values of X_ijk
  YmodPlus = matrix(nrow=d,ncol=n)
  YmodMinus = matrix(nrow=d,ncol=n)
  for (k in 1:d) { #Loop on index of studied variable for impact
    for (i in 1:n) {  #Loop on i in X_ij
      vInd = order(X[i,,k],decreasing=FALSE) #Sort along variable k in X_ij
      #Compute arrays of X_ij and Y_j
      XmodMinus[k,i,1:(N/2),] = X[i,vInd[1:(N/2)],1:d]
      XmodMinus[k,i,(N/2+1):N,] = X[i,vInd[1:(N/2)],1:d]
      YmodMinus[k,i] = mean(mY[i,vInd[1:(N/2)]]) #Y_i is the average of Y_ij for j with smaller component k
      XmodPlus[k,i,1:(N/2),] = X[i,vInd[(N/2+1):N],1:d]
      XmodPlus[k,i,(N/2+1):N,] = X[i,vInd[(N/2+1):N],1:d]
      YmodPlus[k,i] = mean(mY[i,vInd[(N/2+1):N]]) #Y_i is the average of Y_ij for j with larger component k
    }
    #Storing arrays of inputs X_ij and ouputs Y_i
    Xtest = XmodMinus[k,,,]
    save(Xtest,file=paste0("Data/2_InputTestMinus_k",k,".rda"))
    Xtest = XmodPlus[k,,,]
    save(Xtest,file=paste0("Data/2_InputTestPlus_k",k,".rda"))
    Ytest = YmodMinus[k,]
    save(Ytest,file=paste0("Data/2_OutputTestMinus_k",k,".rda"))
    Ytest = YmodPlus[k,]
    save(Ytest,file=paste0("Data/2_OutputTestPlus_k",k,".rda"))
  }
}



runRegression_SW = function(nproj,nCV,nTrainCV,vd,nHypCV) {
  #Once data are created from function generateData
  #this function runs all the sliced Wasserstein regressions in dimensions in vd
  #nproj: number of random projections for sliced Wasserstein
  #nCV: number of cross validation splits between training and test
  #nTrainCV: size of training set
  #vd: vector of dimensions for which distribution regression is done 
  #nHypCV: number of values to be tested for scale parameter and ridge parameter
  #
  #Get outputs for regressions with increasing dimension d:
  load(paste0("Data/1_OutputTrain.rda"))
  load(paste0("Data/1_OutputTest.rda"))
  ntest = length(Ytest)
  n = length(Ytrain)
  MYtest_1 = matrix(nrow=length(vd),ncol=ntest) #For each d, store vector of true Y_i
  MhatYtest_1 = matrix(nrow=length(vd),ncol=ntest) #For each d, store vector of predicted Y_i
  #
  for (id in length(vd):1) { #Loop on dimension d of X_ij
    d = vd[id]
    load(paste0("Data/1_InputTrain_d",d,".rda")) #get corresponding generated array of X_ij
    load(paste0("Data/1_InputTest_d",d,".rda"))
    MYtest_1[id,] = Ytest
    Mproj = MprojSW(nproj,d) #generate random projections for sliced Wasserstein distances
    #Compute n times n matrix of sliced Wasserstein distances for training data:
    MdistTrain = MdistSW(Xtrain,Xtrain,Mproj)  
    #Compute n times ntest matrix of sliced Wasserstein distances between training and test data:
    MdistTrainTest = MdistSW(Xtrain,Xtest,Mproj)
    #vectors of tested values of scale parameter ell and ridge parameter lambda:
    vlc = exp(seq(from=log(median(MdistTrain)/100),to=log(median(MdistTrain)*100),length=nHypCV))
    vlambda = exp(seq(from = log(1e-5),to = 1,length=nHypCV))/n
    #Optimization over the regular grid by cross validation
    res = SelCV(vlc,vlambda,MdistTrain,nCV,nTrainCV,Ytrain)
    #Extracting optimal ell and lambda and final regression
    lc = res$lc
    lambda = res$lambda
    MhatYtest_1[id,] = hatY(lc,lambda,MdistTrain,MdistTrainTest,Ytrain)
  }
  #For each variable k which effect is studied, matrices of true and predicted Y_i
  MhatYtest_2_Minus = matrix(nrow=d,ncol=n)
  MYtest_2_Minus = matrix(nrow=d,ncol=n)
  MhatYtest_2_Plus = matrix(nrow=d,ncol=n)
  MYtest_2_Plus = matrix(nrow=d,ncol=n)
  for (k in 1:d) {  #Loop on variable k
    #Regressions for X_ij and Y_i with smaller values of X_ijk 
    load(paste0("Data/2_InputTestMinus_k",k,".rda"))
    load(paste0("Data/2_OutputTestMinus_k",k,".rda"))
    #Compute n times ntest matrix of sliced Wasserstein distances between training and test data:
    MdistTrainTest = MdistSW(Xtrain,Xtest,Mproj)
    #Note that the corresponding training matrix has already be computed above
    #Sliced Wasserstein regression:
    MhatYtest_2_Minus[k,] = hatY(lc,lambda,MdistTrain,MdistTrainTest,Ytrain)
    #Note that the corresponding optimal parameters ell and lambda have already been computed above
    MYtest_2_Minus[k,] = Ytest
    #Regressions for X_ij and Y_i with larger values of X_ijk
    #with similar comments as above
    load(paste0("Data/2_InputTestPlus_k",k,".rda"))
    load(paste0("Data/2_OutputTestPlus_k",k,".rda"))
    MdistTrainTest = MdistSW(Xtrain,Xtest,Mproj)
    MhatYtest_2_Plus[k,] = hatY(lc,lambda,MdistTrain,MdistTrainTest,Ytrain)
    MYtest_2_Plus[k,] = Ytest
  }
  #return all results: 
  res = list(MYtest_1 = MYtest_1,MhatYtest_1 = MhatYtest_1,
       MYtest_2_Plus = MYtest_2_Plus,MYtest_2_Minus = MYtest_2_Minus,
       MhatYtest_2_Plus = MhatYtest_2_Plus,MhatYtest_2_Minus = MhatYtest_2_Minus)
  return(res)
}










