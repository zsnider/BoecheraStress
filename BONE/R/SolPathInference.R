#' Solution path method
#'
#' Compute probability of the origin and mixture proportions using the whole LASSO solution path. Weighted averages are counted.
#' @param MBapprox A nonsymmetric network computed with the function LASSOSolPath.
#' @param ReferenceData A RUBIAS format baseline data file. Missing values are NA and A=1, C=2, G=3, and T=4 for SNPs.
#' @param SampleNames A character vector of sample names, both mixture and baseline: baseline first, then mixture.
#' @param SampleSizeBaselinePop A numeric value: how many baseline individuals there are?
#' @param alpha A alpha percentage quantile: individuals whose probability of the origin differs from the baseline population proportions less than or equal to this are determined as outgroup members. Note that there are always individuals which are determined as outgroup memebers! Exercise caution and more detailed examination.
#' @return ProbofOrigin = a probability of the origin matrix where individuals are at rows and populations are columns. The last column is the MSE between individuals probability of the origin and baseline proportions. This table is also saved in a folder "NetworkResults" as "ProbOfOriginSolPath.txt"
#'
#' MixtureProp = Mixture proportions. This table is also saved in a folder "NetworkResults" as "MixPropSolPath.txt"
#'
#' OutgroupMembers = A vector of the potential outgroup memebers
#'
#' NmbOfNeighbors = An array of dimension [nmb of baseline populations, nmb of mixture individuals, length(lambda)]. Each row represents the number of neighbors the individual has with the corresponding baseline population individuals.
#'
#' @keywords Network LASSO Baseline Mixture Solution path
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
#'
#' SampleNames=colnames(Y)
#'
#' NetworkResultsSolpath = SolPathInference(MBapprox,ReferenceData,SampleNames=SampleNames,SampleSizeBaselinePop,alpha=0.05)

SolPathInference = function(MBapprox,ReferenceData,SampleNames,SampleSizeBaselinePop,alpha=0.05){
  
  p = ncol(MBapprox[,,1]) # Total number of samples
  
  lambdalength = dim(MBapprox)[3]
  
  SampleSizeMixPop = p - SampleSizeBaselinePop
  
  MixInd = (SampleSizeBaselinePop + 1):p
  BaselineInd = 1:SampleSizeBaselinePop
  
  AllMixConnected = matrix(0,nrow=SampleSizeMixPop,ncol=lambdalength)
  
  lambda = NULL
  
  for(i in 1:lambdalength){
    lambda[i] = paste("lambda",i,sep="") # lambda_1 > lambda_2 > ... lambda_lambdalength
  }
  
  colnames(AllMixConnected) = lambda
  rownames(AllMixConnected) = SampleNames[MixInd]
  
  AllpsiLambda = array(0,c(p,p,lambdalength))
  
  SamplesBaseline = SampleNames[BaselineInd]
  SamplesMix = SampleNames[MixInd]
  
  for(index in 1:lambdalength){
    
    NewMBapprox2 = MBapprox[,,index]
    
    colnames(NewMBapprox2) = SampleNames
    rownames(NewMBapprox2) = SampleNames
    
    NewMBapprox2 = as.matrix(NewMBapprox2)
    
    NewMBapprox2 = NewMBapprox2*t(NewMBapprox2)
    
    psiLambda = ifelse(NewMBapprox2 != 0 , 1, 0)
    
    ###############################################################################
    
    AllMixConnected[,index] = colSums(psiLambda[BaselineInd,MixInd])
    
    AllpsiLambda[,,index] = psiLambda
    
  }
  
  # Determine the probability of the origin for each individual from the LASSO solution path
  
  b = length(unique(ReferenceData$repunit)) # nmb of classes
  
  Freq = array(0,c(b,SampleSizeMixPop,lambdalength)) # How many neighbors individual has with baseline populations
                                                     # given lambda
  
  dimnames(Freq)[[2]] = colnames(psiLambda)[MixInd]
  
  for(index in 1:lambdalength){
    
    Freq[,,index] = rowsum(AllpsiLambda[BaselineInd,MixInd,index],ReferenceData$repunit)
    
  }
  
  # Populations are in order: 
  
  #sort(unique(ReferenceData$repunit))
  
  # H_0: {Baseline frequencies}
  
  prob = c(table(ReferenceData$repunit)/SampleSizeBaselinePop)
  
  classes = sort(unique(ReferenceData$repunit))
  
  # How many times the neighborhood was present over different penalizations?
  
  Proportion = matrix(0,SampleSizeMixPop,b)
  
  colnames(Proportion) = classes
  rownames(Proportion) = SamplesMix
  
  # Determine a weighted average: larger lambda values (heavier penalization) have larger weights:
  
  w = lambdalength:1
  w = w/sum(w)
  
  for(i in 1:lambdalength){
    
    Proportion = Proportion + t(w[i]*Freq[,,i]) # (conditional) Weighted mean/frequency (given baseline populations)
    
  }
  
  RelProportion = apply(Proportion,1,function(x) ifelse(is.na(x/sum(x)),0,x/sum(x))) # Standardized weighted frequencies
  
  RelProportion = t(RelProportion)
  
  DifferenceProportion = sweep(RelProportion,2,prob)
  
  DifferenceProportion = apply(DifferenceProportion,1,function(x) sum(x^2)/length(x))
  
  RelProportion = cbind(RelProportion,DifferenceProportion)
  
  Outsiders = which(DifferenceProportion <= quantile(DifferenceProportion,probs = alpha))
  
  # Write probability of the origin results into a table:
  
  dir.create("NetworkResults", showWarnings = FALSE)
  
  write.table(RelProportion,"NetworkResults/ProbOfOriginSolPath.txt",quote = F)
  
  # Determine also the mixing proportions:
  
  EstimatedMixProp = apply(RelProportion[,1:b],2,as.numeric)
  rownames(EstimatedMixProp) = rownames(RelProportion)
  EstimatedMixProp = apply(EstimatedMixProp[,1:b],2,mean)
  
  EstimatedMixProp = matrix(EstimatedMixProp,b,1,dimnames = list(classes,"Pop	Est.Pi"))
  
  # Write mixing proportions into a file:
  
  write.table(EstimatedMixProp,"NetworkResults/MixPropSolPath.txt",quote = F)
  
  return(list("ProbofOrigin"=RelProportion,"MixtureProp"=EstimatedMixProp,
              "OutgroupMembers"=Outsiders,"NmbOfNeighbors"=Freq))
  
}


