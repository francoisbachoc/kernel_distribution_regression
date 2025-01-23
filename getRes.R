rm(list = ls())
source("functions.R")

load("Data/res_SW.rda")
MR2SW = matrix(nrow=length(vd),ncol=nmc) #matrix of explained variance
   #scores for each ambient dimension d and Monte Carlo repetition 
for (id in 1:length(vd)) { #Loop on dimension 
  for (imc in 1:nmc){  #Loop on Monte Carlo repetition
    #Compute explained variance score:
    MR2SW[id,imc] = 1 - mean((TYtest_1[imc,id,]-ThatYtest_1[imc,id,])^2)/mean((TYtest_1[imc,id,]-mean(TYtest_1[imc,id,]))^2)
  }
}

pdf(file="ecological_simulated_EVS.pdf")
par(mar = c(5,5,5,5))
#Boxplot of explained variance score as function of d
boxplot(t(MR2SW),names=vd,xlab="d",
        ylab= "Explained Variance Score (EVS)",
        ylim=c(0.9,1),cex=2,
        cex.axis=2.4,cex.lab=2.4)
dev.off()
#
#Study of the effect of the variables
for (imc in 1:5) { #Loop on the first Monte Carlo repetitions
  pdf(file=paste0("ecological_simulated_influence_",imc,".pdf"))
  par(mar = c(5,5,5,5))
  #Boxplot of prediction difference as a function of index of studied variable
  boxplot(t(ThatYtest_2_Plus[imc,,] - ThatYtest_2_Minus[imc,,]),xlab="k",ylab= "Prediction difference",ylim=c(-0.4,0.4),
          cex=2,cex.axis=2.4,cex.lab=2.4)
  dev.off()
}







