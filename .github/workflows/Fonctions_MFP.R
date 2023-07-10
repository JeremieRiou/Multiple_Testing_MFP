#Load used packages
library(lattice)
library(CPMCGLM)
library(survival)
library(numbers)
library(SurvRegCensCov)


#Likelihood-ratio test 
LRTest <- function(Z=Z,coup,ind.ajust,codage,data,NamesY1){
  
  tobs <- rep(NA,dim(codage)[3])
  
  for(i in 1:dim(codage)[3]){
    if(ind.ajust==0){
      tobs[i] <- abs(2*(coxph(eval(parse(text=NamesY1)) ~codage[,1:coup[i],i])$loglik[2]-coxph(eval(parse(text=NamesY1)) ~codage[,1:coup[i],i])$loglik[1]))
    } else{
      tobs[i] <- abs(2*(coxph(eval(parse(text=NamesY1)) ~codage[,1:coup[i],i]+Z)$loglik[2]-coxph(eval(parse(text=NamesY1))~Z,data=data)$loglik[2]))
    }
  }
  
  
  return(tobs)
  
}

#Generate data for resampling methods
generate <- function(n,alpha,eta,coeff,Z,X,censor,temps,evt,varcod){
  
  namesZ <- colnames(Z)
  Z <- matrix(Z)
  cov <- cbind(Z,X)
  u <- runif(n)
  v <- runif(n)
  d <- censor/(100-censor)
  Tevt <- ( -log(1-u)/(eta*exp(coeff%*%t(cov))) )^(1/alpha)
  Tcens <- (-(log(1-v))/(d*eta*exp(coeff%*%t(cov))))^(1/alpha)
  delta <- matrix(as.numeric(Tevt < Tcens),1,length(Tevt))
  Ttilde <- t(pmin(Tevt,Tcens))
  product <- cbind.data.frame(Ttilde,t(delta),cov)
  names(product)[1] <- temps
  names(product)[2] <- evt
  names(product)[3:(ncol(Z)+2)] <- namesZ
  names(product)[ncol(Z)+3] <- varcod
  return(product)
} 


#Boostrap function with formula = formula under H0
BootTest <- function(pval=pval,data,NbBoot,formula,coup=coup,codage=X1ter,Z,txcensure,var.cod,temps,evt,ind.ajust,NamesY1,varcod){
  
  pmin <- min(pval)
  pminboot <- rep(0,NbBoot)
  Nbpoly <- dim(codage)[3] 
  matpval <- mattest <- matrix(0,NbBoot,Nbpoly)
  
  modelH0 <- survreg(formula,dist='weibull',data) 
  etahat <- ConvertWeibull(modelH0)$vars[1]
  alphahat <- ConvertWeibull(modelH0)$vars[2]
  betahat <- ConvertWeibull(modelH0)$vars[3:(ncol(Z)+2)]
  
  n <- nrow(data)
  
  for(nbBoot in 1:NbBoot){
    dataBoot <- generate(n,alphahat,etahat,c(betahat,0),Z,var.cod,censor=txcensure,temps,evt,varcod)
    mattest[nbBoot,] <- LRTest(Z=Z,coup,ind.ajust,codage,data=dataBoot,NamesY1)
    
    for (i in 1:Nbpoly) {
      matpval[nbBoot,i] <- (1 - pchisq(mattest[nbBoot,i], coup[i]))
    }
    
  }
  
  pminboot <- apply(matpval,1,min)
  pvalboot <- (1/NbBoot)*sum(pminboot<pmin)
  
  pvalboottransfo <- matrix(NA,nrow=NbBoot,ncol=Nbpoly)
  for(i in 1:NbBoot){
    pvalboottransfo[i,] <- matpval[i,]<pval
  }
  pvalB <- apply(pvalboottransfo,2,mean)
  
  
  return(list(pvalboot,pvalB))
  
}

#Permutation function
PermutTest <- function(pval=pval,data,Nbpermut,X1ter,var.cod=var.cod,coup=coup,ind.ajust,NamesY1,Z){
  
  pmin <- min(pval)
  pminpermut <- rep(0,Nbpermut)
  Nbpoly <- dim(X1ter)[3] 
  matpval <- mattest <- matrix(0,Nbpermut,Nbpoly)
  
  for(nbpermut in 1:Nbpermut){
    
    epermut <- PF(FP, sample(var.cod))
    X5permut <- epermut[[1]]
    X1terpermut <- X5permut
    
    mattest[nbpermut,] <- LRTest(Z=Z,coup,ind.ajust,X1terpermut,data=data,NamesY1)
    
    for(nbpoly in 1:Nbpoly){
      matpval[nbpermut,nbpoly] <- (1 - pchisq(mattest[nbpermut,nbpoly], coup[nbpoly]))    
    }
  }
  pminpermut <- apply(matpval,1,min)
  pvalpermut <- (1/Nbpermut)*sum(pminpermut<pmin)
  
  pvalboottransfo <- matrix(NA,nrow=Nbpermut,ncol=Nbpoly)
  for(i in 1:Nbpermut){
    pvalboottransfo[i,] <- matpval[i,]<pval
  }
  pvalP <- apply(pvalboottransfo,2,mean)
  
  return(list(pvalpermut,pvalP))
  
}


resampling <- function (formula, data, varcod, FP, N = 1000,txcensure,method,alpha){
  
  pos <- which(colnames(data)%in% colnames(data.frame(model.matrix(formula, data))[,-1]))
  if(sum(is.na(data[,pos]))!=0){
    data <- data[-mod(which(is.na(data[,pos])),nrow(data)),]
  }
  
  m <- match.call(expand.dots = FALSE)
  m$varcod <- m$FP <- m$N <- NULL
  m[[1]] <- as.name("model.frame")
  
  NamesY <-  paste0(all.names(update(formula,"~1"))[2],"(",all.names(update(formula,"~1"))[3],",", all.names(update(formula,"~1"))[4],")")
  NamesY1 <-  paste0(all.names(update(formula,"~1"))[2],"(data$",all.names(update(formula,"~1"))[3],", data$", all.names(update(formula,"~1"))[4],")")
  Y <- model.response(model.frame(formula=update(formula,"~1"),data=data))
  evt <- all.names(update(formula,"~1"))[4]
  temps <- all.names(update(formula,"~1"))[3]
  
  mat.exp <- data.frame(model.matrix(formula, data))[,-1]
  
  factor.names <- function(x) {
    x <- matrix(x, nrow = 1)
    Names <- apply(x, MARGIN = 2, FUN = function(x) {
      if (length(grep("factor", x)) != 0) {
        pos1 <- grep("\\(", unlist(strsplit(x,
                                            split = ""))) + 1
        pos2 <- grep("\\)", unlist(strsplit(x,
                                            split = ""))) - 1
        compris.factor <- substr(x, start = pos1, stop = pos2)
        after.factor <- substr(x, start = (pos2 + 2),
                               stop = length(unlist(strsplit(x, split = ""))))
        paste(compris.factor, after.factor, sep = ".")
      }
      else {
        x
      }
    })
    return(Names)
  }
  
  if (is.matrix(mat.exp)) {
    colnames(mat.exp) <- factor.names(colnames(mat.exp))
  }
  if (is.vector(mat.exp)) {
    Z <- NULL
    nb <- 0
    namesZ <- NULL
    var.cod <- as.vector(mat.exp)
    ind.ajust <- 0
  }
  if (!is.vector(mat.exp)) {
    if (!(varcod %in% colnames(mat.exp))) 
      stop("varcod argument it is not present in the dataset.")
    ind.ajust <- 1
    if (ncol(mat.exp) != 2) {
      Z <- mat.exp[, -grep(varcod, colnames(mat.exp))]
      nb <- ncol(Z)
      namesZ <- NULL
      names1 <- unlist(colnames(Z))
      for (i in 1:nb) {
        if (i < nb) 
          namesZ <- paste(namesZ, names1[i], "+")
        else namesZ <- paste(namesZ, colnames(Z)[i])
      }
    }
    else {
      Z <- as.matrix(mat.exp[, -grep(varcod, colnames(mat.exp))])
      nb <- 1
      namesZ <- NULL
      names1 <- as.matrix(colnames(mat.exp))
      posi <- grep(varcod, names1)
      namesZ <- names1[-posi, ]
      colnames(Z) <- namesZ
    }
    var.cod <- as.vector(mat.exp[, varcod])
  }
  n <- nrow(data)
  coup <- NULL
  e <- PF(FP, var.cod)
  X5 <- e[[1]]
  coup5 <- e[[2]]
  X1ter <- X5
  coup <- coup5
  
  if (ind.ajust == 1) {
    formula1 <- update(formula, paste("~", namesZ))
  } else {
    formula1 <- update(formula, paste("~", 1))
  }
  
  
  t.obs <- LRTest(Z=Z,coup,ind.ajust,X1ter,data=data,NamesY1)
  
  pval <- p.val <- NULL
  
  nbtransf1 <- dim(X1ter)[3] 
  
  for (i in 1:nbtransf1) {
    pval[i] <- (1 - pchisq(t.obs[i], coup[i]))
    p.val[i] <- min(p.val[i - 1], pval[i])
  } 
  
  coeff <- as.vector(coxph(eval(parse(text=NamesY1)) ~Z+var.cod)$coefficients)
  
  pval.raw <- pval
  bestPF <- which.min(pval.raw)
  
  
  if(method=="Bootstrap"){
    pval.Bootstrap <- BootTest(pval=pval,data,N,formula1,coup=coup,codage=X1ter,Z,txcensure,var.cod,temps,evt,ind.ajust,NamesY1,varcod)
    pval_fb <- ifelse(pval.Bootstrap[[1]]<0.0001,"<0.0001",pval.Bootstrap[[1]])
    return(paste("Based on bootstrap method,  the fractionnal polynomial number",bestPF, "is selected with a p.value:",pval_fb))
  }
  
  if(method=="Permutation"){
    pval.Permut <- PermutTest(pval=pval,data,N,X1ter,var.cod=var.cod,coup=coup,ind.ajust,NamesY1,Z)
    pval_fp <- ifelse(pval.Permut[[1]]<0.0001,"<0.0001",pval.Bootstrap[[1]])
    return(paste("Based on permutation method,  the fractionnal polynomial number",bestPF, "is selected with a p.value:",pval_fp))
  }   
}

