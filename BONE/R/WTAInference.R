
#' Winner takes it all method
#'
#' By default, the condition is "associate the mixture individual to baseline population where the first neighbor is found"
#' @param MBapprox A nonsymmetric network computed with the function LASSOSolPath.
#' @param ReferenceData A RUBIAS format baseline data file. Missing values are NA. Alleles are coded as integers: A=1, C=2, G=3, and T=4 for SNPs.
#' @param SampleNames A character vector of sample names, both mixture and baseline: baseline first, then mixture.
#' @param SampleSizeBaselinePop A numeric value: how many baseline individuals there are?
#' @param NmbOfParents How many neighbors should be found before we assing individual(s) to corresponding baseline populations? Default is 1.
#' @return ProbofOrigin = A probability of the origin matrix where individuals are at rows and populations are columns. This table is also saved in a folder "NetworkResults" as "ProbOfOriginWTA.txt"
#'
#' MixtureProp = Mixture proportions. This table is also saved in a folder "NetworkResults" as "MixPropWTA.txt".
#'
#' OutgroupSample = IDs of mixture individuals which were not assigned to any baseline population. 
#' 
#' FirstLambda = The index of the penalty parameter lambda when NmbOfParents neighbors were found. This is NA if no neighbors is found.
#' @keywords Network LASSO Baseline 
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
#' NetworkResultsWTA = WTAInference(MBapprox,ReferenceData,SampleNames=SampleNames,SampleSizeBaselinePop,NmbOfParents=1)

WTAInference=function(MBapprox,ReferenceData,SampleNames,SampleSizeBaselinePop,NmbOfParents=1){

  p = ncol(MBapprox[,,1]) # Total number of samples
  
  SampleSizeMixPop = p - SampleSizeBaselinePop
  
  MixInd = (SampleSizeBaselinePop + 1):p # Assuming, that data is ("Baseline columns","Mixture columns")
  BaselineInd = 1:SampleSizeBaselinePop
  
  ############################################################
  
  lambdalength = dim(MBapprox)[3]
  
  AllMixConnected = matrix(0,nrow=SampleSizeMixPop,ncol=lambdalength)
  
  lambda = NULL
  
  for(i in 1:lambdalength){
    lambda[i] = paste("lambda",i,sep="") 
  }
  
  colnames(AllMixConnected) = lambda
  rownames(AllMixConnected) = SampleNames[MixInd]
  
  ############################################################
  
  
  AllpsiLambda = array(0,c(p,p,lambdalength))
  
  for(index in 1:lambdalength){ # Larger graph index means denser graph
    
    NewMBapprox2 = MBapprox[,,index]
    
    colnames(NewMBapprox2) = SampleNames
    rownames(NewMBapprox2) = SampleNames
    
    NewMBapprox2 = as.matrix(NewMBapprox2)
    
    NewMBapprox2 = NewMBapprox2*t(NewMBapprox2)
    
    psiLambda = ifelse(NewMBapprox2 != 0 , 1, 0)
    
    ###############################################################################

    # Count the number of neighbors with a lambda value "lambda_index". Usually less neighbors
    # with smaller lambda = smaller index
    
    AllMixConnected[,index] = colSums(psiLambda[BaselineInd,MixInd])
    
    AllpsiLambda[,,index] = psiLambda
    
  }
  
  dimnames(AllpsiLambda)[[1]] = rownames(psiLambda)
  dimnames(AllpsiLambda)[[2]] = colnames(psiLambda)
  
  #which(rowSums(AllMixConnected) == 0) # No parents
  
  ############################################################
  
  # Find the first adjacency matrix (psiLambda) which fulfils a condition: 
  
  FirstpsiLambda = data.frame(sampleID=rownames(AllMixConnected),index=0)
  
  rownames(FirstpsiLambda) = rownames(AllMixConnected)
  
  for(i in 1:nrow(AllMixConnected)){
    
    FirstpsiLambda[i,2] =  which(AllMixConnected[i,] >= NmbOfParents)[1] # Try different conditions
    
  }
  
  # Are there individuals with no baseline population?
  
  #sum(is.na(FirstpsiLambda[,2]))
  
  NoFounds = which(is.na(FirstpsiLambda[,2]))
  OutgroupSample = as.character(FirstpsiLambda[NoFounds,1])
  
  ############################################################
  
  # Collect neighbours to a summary matrix:
  
  SummarypsiLambda = matrix(0,p,p)
  
  rownames(SummarypsiLambda) = dimnames(AllpsiLambda)[[1]]
  colnames(SummarypsiLambda) = dimnames(AllpsiLambda)[[2]]
  
  FirstpsiLambda[NoFounds,2] = 1
  
  for(i in 1:nrow(AllMixConnected)){
    
    # Neighbourhood of Baseline:
    
    SummarypsiLambda[,as.character(FirstpsiLambda[i,1])] = 
      AllpsiLambda[,as.character(FirstpsiLambda[i,1]),FirstpsiLambda[i,2]]
    
    # Neighbourhood of Mixture:
    
    SummarypsiLambda[as.character(FirstpsiLambda[i,1]),] = 
      AllpsiLambda[as.character(FirstpsiLambda[i,1]),,FirstpsiLambda[i,2]]
    
  }
  
  write.table(SummarypsiLambda,file="SummarypsiLambda.txt",quote=F)
  
  ########################################################################
  ########################################################################
  
  b = length(unique(ReferenceData$repunit)) # nmb of classes
    
  Freq = rowsum(SummarypsiLambda[BaselineInd,MixInd],ReferenceData$repunit)
  
  dimnames(Freq)[[2]] = colnames(psiLambda)[MixInd]
  
  # Populations are in order: 
  
  #sort(unique(ReferenceData$repunit))
  
  classes = sort(unique(ReferenceData$repunit))
  
  # How many times the neighborhood was present over different penalizations?
  
  RelProportion = apply(Freq,2,function(x) ifelse(is.na(x/sum(x)),0,x/sum(x))) # Rescaled weighted frequencies. If no population is found, set to zero.
  
  RelProportion = t(RelProportion)
  
  # Write probability of the origin results into a table:
  
  dir.create("NetworkResults", showWarnings = FALSE)
  
  write.table(RelProportion,"NetworkResults/ProbOfOriginWTA.txt",quote = F)
  
  # Determine also the mixing proportions:
  
  EstimatedMixProp = apply(RelProportion[,1:b],2,as.numeric)
  rownames(EstimatedMixProp) = rownames(RelProportion)
  EstimatedMixProp = apply(EstimatedMixProp[,1:b],2,mean)
  
  EstimatedMixProp = matrix(EstimatedMixProp,b,1,dimnames = list(classes,"Pop	Est.Pi"))
  
  # Return also the lambda index where the first neighbor was found:
  
  FirstpsiLambda[NoFounds,2] = NA
  
  # Write mixing proportions into a file:
  
  write.table(EstimatedMixProp,"NetworkResults/MixPropWTA.txt",quote = F)
  
  return(list("ProbofOrigin"=RelProportion,"MixtureProp"=EstimatedMixProp,"OutgroupSample"=OutgroupSample,"FirstLambda"=FirstpsiLambda))
  
}