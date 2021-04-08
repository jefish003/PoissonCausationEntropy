library(glasso)
library(Matrix)
library(MASS)
library(EnvStats)
library(rdetools)



ReadFolder = "/Data/"
RData = "Data_ER_0.05_nsamp_1200_"
ReadDate = "08-Apr-2021.csv"

WD = dirname(rstudioapi::getSourceEditorContext()$path)
ReadName = paste(ReadFolder,RData,ReadDate,sep="")

WriteDate = Sys.Date()
WriteFront = "Glasso_Results_"
WriteName = paste(ReadFolder,WriteFront,RData,WriteDate,".csv",sep="")

CSV = read.csv(paste(WD,"/",ReadName,sep=""),head=FALSE)

Mat = as.matrix(CSV,nrows=nrow(CSV),ncols=50)
BC = boxcox(CSV[,1]+1,objective.name="Log-Likelihood", optimize=TRUE)
lambda = BC[1]
lambda = as.double(lambda)
M = as.matrix(CSV[,1])
M = M+1
M1 = (M^lambda)/lambda
LS = logspace(-2,0,1000)

for(i in 2:50){
  print(i)
  BC = boxcox(CSV[,i]+1,objective.name="Log-Likelihood", optimize=TRUE)
  lambda = BC[1]
  lambda = as.double(lambda)
  M = as.matrix(CSV[,i])
  M = M+1
  M1 = cbind(M1,(M^lambda)/lambda)
}

S = var(M1)

BIC = -1e100
for(j in 1:1000){
  alpha = LS[j]
  GL = glasso(S,alpha,nobs=nrow(CSV))
  LogLikeliNum = as.double(GL[3])
  Adj = matrix(unlist(GL[2]),nrow=50,byrow=TRUE)
  Wh = which(Adj!=0)
  Adj[Wh] = 1
  B2 = length(Wh)*log(nrow(CSV))/2
  print(paste(i,BIC))
  #BIC Maximization Here
  if(LogLikeliNum-B2>BIC){
    A = Adj
    BIC = LogLikeliNum-B2
    print(BIC)
  }
}

write.table(A,paste(WD,"/",WriteName,sep=""),sep=",",row.names=FALSE,col.names=FALSE)

