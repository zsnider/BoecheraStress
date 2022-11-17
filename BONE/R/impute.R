#' Compute a genotype matrix and impute missing genotypes
#'
#' Compute a genotype matrix. One can also impute missing genotypes once either with marker mode (row mode) or random selection of gentype based on the genotype frequencies across the entire baseline.
#' @param impute Should missing genotypes be imputed. Default TRUE.
#' @param Data A RUBIAS format data table. Missing values are NA. Alleles are coded as integers: A=1, C=2, G=3, and T=4 for SNPs.
#' @param rmMono Should monomorphic SNPs be removed. Default TRUE.
#' @param rmMissing Shold SNPs with more than alpha percentage missing values be removed. Default TRUE. 
#' @param method Imputation method used. If "mode" impute missing values with marker mode (row mode). If "random" impute missing values by selection of genotype based on the genotype frequencies across all populations. Default set as "mode".
#' @param genotyping The genotyping method used. If set as "Axiom", the number of reference alleles (C/G) is returned. If set as "Fluidigm", the original allele coding is preserved and four "genotypes" are returned (alleles are divided into independent observations at each locus).
#' @param alpha A user defined value: Remove genotypes with more than alpha percentage values missing. Default set as 5 percentage. 
#' @param verbose Tell when it is done and everything went smoothly. Default TRUE.
#' @return Returns a (nmb of SNPs)x(nmb of individuals) numeric valued data matrix. One should impute baseline and mixture data separately.
#' @keywords Impute monomorphic missing values genotypes
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

impute = function(Data, impute=T, rmMono=T, rmMissing=T, method=NULL, genotyping=NULL, alpha=0.05, verbose=T){
  
  if(impute == F) NmbOfRm = NmbOfMono = NmbOfImputed = 0
  
  if(is.null(method)) method = "mode"
  
  if(is.null(genotyping)) genotyping = "Axiom"
  
  IndivNames = Data$indiv
  
  Yraw = Data[,-c(1:4)]
  
  LocusNames = colnames(Yraw)
  
  #dim(Yraw)
  
  Yraw = matrix(as.numeric(unlist(Yraw)),nrow=nrow(Yraw))
  
  rownames(Yraw) = IndivNames
  colnames(Yraw) = LocusNames
  
  if(genotyping == "Axiom"){
    
    # A=1, C=2, G=3, and T=4. NA for missing
    
    # Change: A = 0, C = 1, G = 1, T = 0
    
    Yraw[is.na(Yraw)] <- 9 # Missing values
    Yraw[Yraw == 1] <- 0
    Yraw[Yraw == 4] <- 0 
    Yraw[Yraw == 3] <- 1
    Yraw[Yraw == 2] <- 1
    
    nmbofSNPs = ncol(Yraw)/2
    
    a = 2*(0:nmbofSNPs)
    
    b = a + 1
    
    a = a[-1]; b = b[-length(b)]
    
    YFirst = Yraw[,b]
    YSecond = Yraw[,a]
    
    #table(YFirst)
    #table(YSecond)
    
    Y = YFirst + YSecond
	
	}
	
	if(genotyping == "Fluidigm"){
    
		Yraw[is.na(Yraw)] <- 9 # Missing values
		
		Y = Yraw
	
	}
    
  rm(Yraw)
    
  Y = t(Y) # Individuals at columns
    
  if(impute == T){
  
	##################################################################
    
  # Remove SNPs with more than alpha% genotypes missing:
    
  if(rmMissing){
  
	if(genotyping == "Axiom") fun = function(a){ sum(a == 18)/ncol(Y)}
	if(genotyping == "Fluidigm") fun = function(a){ sum(a == 9)/ncol(Y)}
    
	missing = apply(Y,1,fun)
    
	Ind = which(missing > alpha)
    
	NmbOfRm = length(Ind)
    
	if(length(Ind) > 0){
		Y = Y[-Ind,] 
	}
  
  }
    
  ##################################################################
    
  # Remove monomorphic SNPs:
  
  if(rmMono){
  
	fun = function(a){ length(table(as.vector(a[a != 9]))) }
    
	monomorphic = apply(Y,1,fun)
    
	Ind = which(monomorphic == 1)
    
	NmbOfMono = length(Ind)
    
	if(length(Ind) > 0){
		Y = Y[-Ind,]
	}
  
  }
    
  ###############################################################
    
  # Impute missing values either with marker moder (row mode), with the rarest SNP or 
  # random selection of genotype based on the allele frequencies across the entire baseline:
    
  if(method == "mode"){
      
    if(genotyping == "Axiom") fun = function(a){ length(table(as.vector(a[a == 18]))) }
	  if(genotyping == "Fluidigm") fun = function(a){ length(table(as.vector(a[a == 9]))) }
      
    howmanystillmissing = apply(Y,1,fun)
      
    whichmissing = which(howmanystillmissing == 1)
      
    NmbOfImputed = length(whichmissing)
      
    for(i in whichmissing){
        
      tab = table(as.vector(Y[i,]))
        
      if(genotyping == "Axiom") Y[i,Y[i,] == 18] <- as.numeric(names(tab)[tab == max(tab)])[1]
		  if(genotyping == "Fluidigm") Y[i,Y[i,] == 9] <- as.numeric(names(tab)[tab == max(tab)])[1]
        
    }
  }
    
  if(method == "random"){
      
    if(genotyping == "Axiom") fun = function(a){ length(table(as.vector(a[a == 18]))) }
	  if(genotyping == "Fluidigm") fun = function(a){ length(table(as.vector(a[a == 9]))) }
      
    howmanystillmissing = apply(Y,1,fun)
      
    whichmissing = which(howmanystillmissing == 1)
      
    NmbOfImputed = length(whichmissing)
      
    a = table(Y)
      
	  if(genotyping == "Axiom"){
	  
		  b = a[names(a) == 18]
      
		  a = a[names(a) != 18]
      
		  rho = a/sum(a) 
      
		  Y[Y == 18] = sample(as.numeric(names(a)), b, prob=rho, replace = T)
	  
	   }
	  
	   if(genotyping == "Fluidigm"){
	  
		  b = a[names(a) == 9]
      
		  a = a[names(a) != 9]
      
		  rho = a/sum(a) 
      
		  Y[Y == 9] = sample(as.numeric(names(a)), b, prob=rho, replace = T)
	  
	   }
      
  }
  
  }
    
  ###############################################################
  
  if(verbose){
    cat("Done.")
  }
  
  return(list("Y"=Y,"Removed"=c(NmbOfRm,NmbOfMono,NmbOfImputed)))
  
}