fit.CP_alt <- function(CountMatrix, Strains, test = FALSE, optimizer = "nlminb"){
  # Fit a compound Poisson mixed effect model for one or a list of features and 
  #   output the fit parameters. This is an alternative to the function 
  #   fit.CP() in the HeritSeq package. It includes additional command 
  #   detach(cplm, unload = TRUE) & library(cplm) to avoid potential memory 
  #   issue when fitting multiple features.  
  #
  # [INPUT]
  # CountMatrix: Sequencing count matrix for one or more features. Each 
  #   row is for one feature, and the columns are for samples. 
  # Strains: Strain labels for the samples. 
  # test: TRUE or FALSE (default). Test the presence of heritability 
  #   through examining if the random effect variance is zero.
  # optimizer: A character string that determines which optimization 
  #   routine is to be used. Possible choices are "nlminb" (default), 
  #   "L-BFGS-B", and "bobyqa". 
  #
  # [OUTPUT]
  # A list with two objects. The first object is a matrix indicating the 
  # fitted parameters for each feature. The columns are ordered by intercept, 
  # tweedie parameter, random effect variance, and dispersion. Row names are    
  # feature names. If the argument test is set to
  # be true, the second object of the list consists of p-values for testing 
  # the hypothesis that random effects is zero; otherwise, the second object 
  # is NULL. 
  
  
  if(is.null(dim(CountMatrix))){
    print('Fitting a single feature.')
    CountMatrix <- matrix(CountMatrix, nrow = 1)
  }
  
  paras <- t(pbapply::pbsapply(1:nrow(CountMatrix), function(x){
    library(cplm)
    
    CountVector <- CountMatrix[x, ]
    dat_sub <- data.frame(expr = as.numeric(CountVector), strain = Strains)
    
    fit <- tryCatch({
      fit1 <- cpglmm(expr ~ 1 + (1|strain), data = dat_sub, 
                     optimizer = optimizer)
    }, error=function(err){
      fit1 <- try({cpglmm(expr ~ 1 + (1|strain), data = dat_sub, 
                          optimizer = optimizer)}) 
      return(fit1)
    })
    
    if (class(fit) != "try-error"){
      as <- fit$fixef
      sigma_a2 <- as.numeric(cplm::VarCorr(fit)$strain)
      p <- fit$p
      phi <- fit$phi
      
      para <- c(as, sigma_a2, p, phi)
    }else{
      print(paste("Fitting problem for feature", x, "returning NA"))
      para <- rep(NA, 4)
    }
    
    ### Fitting the reduced model for testing significance of the random effect ###
    if (test){
      detach("package:cplm", unload = TRUE)
      suppressMessages(suppressWarnings(library(cplm)))
      attachNamespace("cplm")
      fit.red <- tryCatch({
        fit2 <- cplm::cpglm(expr ~ 1, data = dat_sub, optimizer = optimizer)
      }, error=function(err){
        fit2 <- try({cplm::cpglm(expr ~ 1, data = dat_sub, optimizer = optimizer)})
        return(fit2)
      })
      
      if (class(fit.red) != "try-error" & class(fit) != "try-error"){
        test.stat <- 2*cplm::logLik(fit)+cplm::AIC(fit.red)-2
        
        if (test.stat<1e-6) {test.stat <- 0}
        pval <- 0.5*pchisq(test.stat, df = 1, lower.tail = FALSE) + 
          0.5*as.numeric(test.stat == 0)
      }else{
        print(paste("Cannot do test for feature", x, "fitting problem."))
        pval <- NA
      }
      
      para <- c(para, pval)
    }
    
    detach("package:cplm", unload = TRUE)
    suppressMessages(suppressWarnings(library(cplm)))
    
    return(para)
  }))
  
  paras1 <- matrix(paras[,1:4], ncol = 4)
  rownames(paras1) <- rownames(CountMatrix)  
  colnames(paras1) <- c("alpha_g", "sigma2_g", "p_g", "phi_g")
  
  if (test){
    return(list(paras = paras1,  pvals = paras[ , 5]))
  }else{
    return(list(paras = paras1,  pvals = NULL))
  }
  
}

