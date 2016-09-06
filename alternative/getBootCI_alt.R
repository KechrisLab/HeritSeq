# Compute variance partition coefficition (VPC) confidence intervals (CI) 
# for one or more features.
# 
# Compute VPC CI based on parametric bootstrap for one or more features.
#
# [INPUT]
# CountMatrix: A G by N count matrix. G is the number of features; N is the 
#   total number of samples.
# Strains: A 1 by N vector of strain labels corresponding to each sample.
# which.features: A 1 by k vector of select feature numbers for which CI is 
#   desired. k < G.
# num.boot: Number of bootstraps.
# method: Which method should be used, "CP-fit", "NB-fit" (default), or "VST".
#   "VST" method bootstraps data under negative binomial mixed models. 
#   Infrequently, the "CP-fit" method might produce inaccurate CI due to the
#   potential memory issue in the cplm pakcage.
# alpha: A numerical value between 0 and 1, indicating the significance level 
#   of the CI. Default value is 0.05.
# optimizer: A character string that determines which optimization routine is 
#   to be used. It is only used for method = "CP-fit". Possible choices are
#   "nlminb" (default), "L-BFGS-B", and "bobyqa".
# 
# [OUTPUT]
# A list of two objects. The first object is a k by 2 matrix containing the CI.
#   The second object consists of a k by num.boot matrix of all bootsrapped VPC
#   values.



getBootCI = function(CountMatrix, Strains, which.features, num.boot,
                     method="NB-fit", alpha=0.05, optimizer = "nlminb"){
  all.vpcs = matrix(NA, nrow = length(which.features), ncol=num.boot)
  vec.num.rep = as.numeric(table(Strains))
  
  if (method=="NB-fit"){
    print("Getting initial point estimates using the NB-fit method")
    fit = fit.NB(CountMatrix[which.features,], Strains)
    
    for (i in 1:length(which.features)){
      print(paste("Bootstraping feature",i))
      boot.data = getReadMatrix.NB(vec.num.rep, 
                                   rep(fit$paras[i,1],num.boot), 
                                   rep(fit$paras[i,2],num.boot),
                                   rep(fit$paras[i,3],num.boot))
      fit.i = fit.NB(boot.data, Strains)
      all.vpcs[i,] = computeVPC.NB(fit.i$paras)
    }
    
    intervals = cbind( apply(all.vpcs, 1, quantile, probs = alpha/2),
                       apply(all.vpcs, 1, quantile, probs = 1-alpha/2))
    
    return(list(intervals = intervals, all.vpcs = all.vpcs))
  }
  
  
  
  
  if (method=="CP-fit"){
    print("Getting initial point estimates using the CP-fit method")
    fit = fit.CP(CountMatrix[which.features,], Strains, optimizer = optimizer)
    
    for (i in 1:length(which.features)){
      print(paste("Bootstraping feature",i))
      boot.data = getReadMatrix.CP(vec.num.rep, 
                                   rep(fit$paras[i,1],num.boot), 
                                   rep(fit$paras[i,2],num.boot),
                                   rep(fit$paras[i,3],num.boot),
                                   rep(fit$paras[i,4],num.boot))
      fit.i = fit.CP_alt(boot.data, Strains, optimizer = optimizer)
      all.vpcs[i,] = computeVPC.CP(fit.i$paras)
    }
    
    intervals = cbind( apply(all.vpcs, 1, quantile, probs = alpha/2),
                       apply(all.vpcs, 1, quantile, probs = 1-alpha/2))
    
    return(list(intervals = intervals, all.vpcs = all.vpcs))
  }
  
  
  
  
  if (method=="VST"){
    print("Getting initial point estimates using the VST method")
    fit = fit.NB(CountMatrix[which.features,], Strains)
    
    for (i in 1:length(which.features)){
      print(paste("Bootstraping feature",i))
      boot.data = getReadMatrix.NB(vec.num.rep, 
                                   rep(fit$paras[i,1],num.boot), 
                                   rep(fit$paras[i,2],num.boot),
                                   rep(fit$paras[i,3],num.boot))
      
      cds <- DESeq2::DESeqDataSetFromMatrix(round(boot.data), 
                                            data.frame(strain = Strains), 
                                            formula(~strain))
      cds <- DESeq2::estimateSizeFactors(cds)
      
      if (sum(is.na(sizeFactors(cds)))>0){   
        cds <- DESeq2::DESeqDataSetFromMatrix(round(1+boot.data), 
                                              data.frame(strain = Strains), 
                                              formula(~strain)) 
        # 1 added to avoid too many low counts in those cases
        # This adds a small variation to the data
        cds <- DESeq2::estimateSizeFactors(cds)
      }
      cds <- DESeq2::estimateDispersions(cds, fitType = "local")
      vsd <- DESeq2::varianceStabilizingTransformation(cds, fitType = "local")
      vsd <- SummarizedExperiment::assay(vsd)
      
      all.vpcs[i,] = fitComputeVPC.lmer(vsd, Strains)$vpcs
    }
    
    intervals = cbind( apply(all.vpcs, 1, quantile, probs = alpha/2),
                       apply(all.vpcs, 1, quantile, probs = 1-alpha/2))
    
    return(list(intervals = intervals, all.vpcs = all.vpcs))
  }
}



