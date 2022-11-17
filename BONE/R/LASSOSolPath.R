#' Compute LASSO solution path
#'
#' Computes the LASSO solution path of the length length(lambda) using multinomial logit-model which is augmented with a LASSO-penalty function.
#' @param Y A (nmb of SNPs)x(nmb of individuals) data matrix with numeric entries 0,1 and 2. No missing values allowed. First "SampleSizeBaselinePop" columns are SNPs of baseline individuals, rest are mixture individuals. 
#' @param lambda A user defined vector of non-zero tuning parameter values in decreasing order.
#' @param intercept Should the intercept be included in the logit-model. Default TRUE.
#' @param SampleSizeBaselinePop How many baseline individuals are in data. A numeric value.
#' @param Baseline Should it also be checked how baseline individuals depend on mixture individual. Default TRUE.
#' @param UseWeights Genotype weights (observation weights of the multinomial logit model). If TRUE, the procedure computes weights so that more weight is given to rare genotypes (computed in baseline and mixture group separately). When FALSE (default) weight is 1 for each genotype.
#' @param verbose Should process be printed. The default value is TRUE.
#' @return Returns an array of dimension [ncol(Y),ncol(Y),length(lambda)] and saves the array in the working directory. ith matrix (a non-symmetric adjacency matrix) in the array corresponds to an adjacency matrix associated with the tuning parameter value lambda_i. 
#' @keywords Network LASSO Baseline Mixture Solution path glmnet
#' @export
#' @examples
#' MixtureY = impute(MixtureData) # Simple marker mode imputation
#' ReferenceY = impute(ReferenceData)
#'
#' # Choose only same set of markers:
#' setdiff(rownames(MixtureY),rownames(ReferenceY))
#' Markers = intersect(rownames(MixtureY),rownames(ReferenceY))
#'
#' MixtureY = MixtureY[Markers,]
#' ReferenceY = ReferenceY[Markers,]
#'
#' SampleSizeBaselinePop = ncol(ReferenceY)
#' Y = cbind(ReferenceY,MixtureY)
#'
#' lambda = seq(0.4,0.02,length.out=40)
#'
#' MBapprox = LASSOSolPath(Y,lambda,intercept=T,SampleSizeBaselinePop,Baseline=T)

LASSOSolPath = function(Y, lambda, intercept=T, SampleSizeBaselinePop, UseWeights = F, Baseline=T, verbose=T){
  
  L = nrow(Y)
  N = ncol(Y)
  
  p = N
  
  MBapprox = array(0,c(p,p,length(lambda)))
  
  k = 1
  
  MixInd = (SampleSizeBaselinePop + 1):N
  BaselineInd = 1:SampleSizeBaselinePop
  
  h = 1
  
  cat("Mixture","\n")
  
  for(i in MixInd){
    
    if(UseWeights){
      
      GenFreq = table(Y[ ,MixInd])/length(Y[ ,MixInd])
      
      weights = 1 - GenFreq[as.character(Y[ ,i])]
      
      weights = as.vector(weights)
     
      fit = glmnet(Y[ ,BaselineInd], Y[ ,i], intercept=intercept, lambda=lambda, family="multinomial", weights = weights)
       
    }else{
      
      fit = glmnet(Y[,BaselineInd], Y[,i],intercept=intercept, lambda=lambda, family="multinomial")
      
    }
    
    for(j in 1:length(lambda)){
      
      BetaNames = names(fit$beta)
      
      sarake = 0
      
      for(d in 1:length(BetaNames)){
        
        sarake = sarake + as.vector(abs(coef(fit)[[BetaNames[d]]][,j]))
        
      }
      
      MBapprox[BaselineInd,i,j] = sarake[-1]
      
    }
    
    if(verbose) cat("\r", round(h/length(MixInd)*100,2), "%")
    
    h = h + 1
    
  }
  
  if(Baseline==F) MBapprox = MBapprox + aperm(MBapprox, c(2,1,3))
  
  # Then, explain baseline with mixture
  
  if(Baseline){
  
	cat("\n")
    cat("Baseline","\n")
  
	for(i in BaselineInd){
    
		if(UseWeights){
      
			GenFreq = table(Y[ ,BaselineInd])/length(Y[ ,BaselineInd])
      
			weights = 1 - GenFreq[as.character(Y[ ,i])]
      
			weights = as.vector(weights)
      
			fit = glmnet(Y[ ,MixInd], Y[ ,i], intercept=intercept, lambda=lambda, family="multinomial", weights = weights)
      
		}else{
      
			fit = glmnet(Y[,MixInd], Y[,i],intercept=intercept,lambda=lambda,family="multinomial") 
      
    }
    
    for(j in 1:length(lambda)){
      
      BetaNames = names(fit$beta)
      
      rivi = 0
      
      for(d in 1:length(BetaNames)){
        
        rivi = rivi + as.vector(abs(coef(fit)[[BetaNames[d]]][,j]))
        
      }
      
      MBapprox[MixInd,i,j] = rivi[-1]
      
    }
    
    if(verbose) cat("\r", round(i/length(BaselineInd)*100,2), "%")
    
  }
  
  }
  
  #########################################################
  
  #########################################################
  
  dir.create("NetworkResults", showWarnings = FALSE)
  
  save(MBapprox,file="NetworkResults/SolutionPathNonSymmetricNet.RData")
  
  return(MBapprox) # Hox! Non-symmetric network if Baseline = FALSE !
  
}