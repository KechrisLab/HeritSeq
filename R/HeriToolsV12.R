
############################ Functions for heritability analysis ##############################



######## Packages that we need ########

# You only need to install once #

#install.packages("lme4", repos="http://cran.r-project.org")
#install.packages("cplm", repos="http://cran.r-project.org")
#install.packages("pbapply", repos="http://cran.r-project.org")
#install.packages("R2admb")
#install.packages("glmmADMB", 
#                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                         getOption("repos")),
#                 type="source")
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq")


requireNamespace("lme4")
requireNamespace("glmmADMB")
requireNamespace("tweedie")
requireNamespace("cplm")
requireNamespace("pbapply")
requireNamespace("DESeq")

if(getRversion() >= "3.1.0") utils::globalVariables(c("cplm.data", "cplm"))

###############################################################################
### Generate negative binomial distributed data matrix 

# (INTERNAL)
getNBReads <- function(vec.num.rep, beta0, sig2.strain, phi){
  # Generate possibly unbalanced counts from a NB mixed effect model:
  # log(E(Yij)) = beta0 + Straini 
  # Straini ~ N(0, sig2.strain)
  #
  # vec.num.rep: a vector of replicate numbers for each strain.
  # phi: dispersion parameter in the NB model; common across 
  #   strains. The greater the value of this parameter, the greater is the variance.
  # sig2.strain: variance of the strain random effect.
  # beta0: intercept term in the GLMM.
  #
  # [OUTPUT]
  # NBcounts: a 1 by N matrix with NB counts. N is the total number of samples.
  #   Column names are sample names of the form "Ss_r", where S stands for 
  #   sample, s is the strain number, r is the replicate number within the 
  #   strain. 
  
  if(sig2.strain == 0){
    warning("No random effect.")
  }
  if(phi <= 0){
    stop("Invalid dispersion value.")
  }
  
  num.strains <- length(vec.num.rep)
  strain.means <- exp(beta0 + rnorm(num.strains, sd = sqrt(sig2.strain)))
  
  NBcounts <- lapply(1:num.strains, function(x){
    counts.x <- rnegbin(n = vec.num.rep[x],
                        mu = strain.means[x],
                        theta = 1/phi)
    return(counts.x)
  })
  NBcounts <- matrix(do.call(c, NBcounts), nrow = 1)
  
  sample.names <- lapply(1:num.strains, function(x){
    sample.x <- paste0("S", x, "_", 1:(vec.num.rep[x]))
    return(sample.x)
  })
  sample.names <- as.vector(do.call(c, sample.names))
  colnames(NBcounts) <- sample.names
  
  return(NBcounts)
}

#' Generate a count matrix from NBMM models.
#' 
#' Generate a (possibly unbalanced) count matrix from
#' negative binomial mixed effect models.
#' Each row corresponds to one feature that follows: \cr
#' \eqn{log(E(Y_{ij})) = {\beta}_0 + Strain_i}{log(E(Yij)) = beta0 + Strain_i} \cr
#' \eqn{Strain_i \sim N(0, sig2.strain)}{Strain_i ~ N(0, sig2.strain)}
#' 
#' @param vec.num.rep A vector of replicate numbers for each strain.
#' @param beta0s Vector for the intercepts in the GLMM's.
#' @param sig2.strains Variance vector of the strain random effect.
#' @param phis Dispersion vector in the NB models; common across strains per feature.
#' @return A \eqn{G \times N}{G x N} matrix with NB counts. \eqn{N} is the total number of 
#'  samples; \eqn{G} is the number of features. Column names are sample names
#'  of the form "Ss_r", where S stands for sample, s is the strain number,
#'  r is the replicate number within the strain. Row names are the feature
#'  names of the form "Gene g", where g is the feature index.
#' @examples
#' ## Generate a sequencing dataset with 5 features and 6 strains. 
#' ## Assign parameter values.
#' rep.num <- c(3, 5, 2, 3, 4, 2)
#' b0 <- c(-1, 1, 2, 5, 10)
#' sig2s <- c(10, 0.2, 0.1, 0.03, 0.01)
#' phis <- c(0.5, 1, 0.05, 0.01, 0.1)
#' 
#' set.seed(1234)
#' ## Generate reads:
#' nbData <- getNBReadMatrix(rep.num, b0, sig2s, phis)
#' @export
getNBReadMatrix <- function(vec.num.rep, beta0s, sig2.strains, phis){

  num.probes <- length(beta0s)
  CountMatrix <- lapply(1:num.probes, function(x){
    probe.x <- getNBReads(vec.num.rep, beta0s[x], 
                          sig2.strains[x], phis[x])
  })
  CountMatrix <- do.call(rbind, CountMatrix)
  rownames(CountMatrix) <- paste0("Gene ", 1:num.probes)
  return(CountMatrix)
}



###############################################################################
### Generate compound Poisson distributed data matrix 

# (INTERNAL)
getCPReads <- function(vec.num.rep, beta0, sig2.strains, p, phi){   
  # Generate possibly unbalanced reads from a CP mixed effect model:
  # log(mu) = beta0 + Straini 
  # Straini ~ N(0, sqrt(sig2.strains))
  #
  # vec.num.rep: a vector of replicate numbers for each strain.
  # beta0: intercept.
  # sig2.strains: random effect variance. 
  # p: Power parameter in CP models.
  # phi: Dispersion parameter in CP models.
  #
  # [OUTPUT]
  # CPcounts: a 1 by N matrix with CP reads. N is the total number of samples.
  #   Column names are sample names of the form "Ss_r", where S stands for 
  #   sample, s is the strain number, r is the replicate number within the 
  #   strain. 
  
  
  if(abs(p - 1.5) >= 0.5){
    stop("The power parameter p needs to satisfy 1<p<2.")
  }
  if(sig2.strains == 0){
    warning("No random effect.")
  }
  if(phi <= 0){
    stop("Invalid dispersion value.")
  }
  
  num.strains <- length(vec.num.rep)
  mus <- exp(beta0 + rnorm(num.strains, sd = sqrt(sig2.strains)))
  
  CPcounts <- lapply(1:length(mus), function(x){
    counts.x <- rtweedie(vec.num.rep[x], xi = p, mu = mus[x], phi = phi)
  })
  CPcounts <- matrix(do.call(c, CPcounts), nrow = 1)
  
  sample.names <- lapply(1:num.strains, function(x){
    sample.x <- paste0("S", x, "_", 1:(vec.num.rep[x]))
    return(sample.x)
  })
  sample.names <- as.vector(do.call(c, sample.names))
  
  colnames(CPcounts) <- sample.names
  
  return(CPcounts)
}


#' Generate a read matrix from a CPMM.
#' 
#' Generate a (possibly unbalanced) read matrix from a CP mixed effect model.
#' Model: \cr
#' \eqn{log(\mu) = {\beta}_0 + Strain_i}{log(mu) = beta0 + Strain_i} \cr
#' \eqn{Strain_i \sim N(0, sig2.strains)}{Strain_i ~ N(0, sig2.strains)}
#' 
#' @param vec.num.rep A vector of replicate numbers for each strain.
#' @param beta0s Intercept vector, \eqn{1 \times \texttt{num.features}}{1 x num.features}.
#' @param sig2.strains Random effect variance vector, \eqn{1 \times \texttt{num.features}}{1 x num.features}.
#' @param ps Power parameter in CP models, a \eqn{1 \times \texttt{num.features}}{1 x num.features} vector.
#' @param phis Dispersion parameter in CP models, a \eqn{1 \times \texttt{num.features}}{1 x num.features} vector.
#' @return CountMatrix: a \eqn{G \times N}{G x N} matrix with CP reads. \eqn{N} is the total number of 
#'   samples; \eqn{G} is the number of features. Column names are sample names 
#'   of the form "Ss_r", where S stands for sample, s is the strain number, 
#'   r is the replicate number within the strain. Row names are the feature 
#'   names of the form "Gene g", where g is the feature index.
#' @examples
#' ## Generate a sequencing dataset with 5 features and 6 strains. 
#' ## Assign parameter values.
#' rep.num <- c(3, 5, 2, 3, 4, 2)
#' b0 <- c(-1, 1, 2, 5, 10)
#' sig2s <- c(10, 0.2, 0.1, 0.03, 0.01)
#' ps <- rep(1.5, 5)
#' phis <- c(1.5, 1, 0.5, 0.1, 0.1)
#' 
#' set.seed(1234)
#' ## Generate reads:
#' cpData <- getCPReadMatrix(rep.num, b0, sig2s, ps, phis)
#' ## Generate strain names:
#' str <- sapply(1:length(rep.num), function(x){
#'   str.x <- paste0("S", x)
#'   return(rep(str.x, rep.num[x]))
#' })
#' str <- do.call(c, str)
#' 
#' ## Visualize sequencing reads for one feature.
#' require(ggplot2)
#' require(reshape2)
#' CountVector <- cpData[1, ]
#' raw_melt <- melt(CountVector)
#' raw_melt$Var <- str
#' names(raw_melt) <- c('read', 'strain')
#' ggplot(raw_melt) +
#'   aes(x = reorder(strain, read, FUN = mean), y = read) +
#'   geom_boxplot(aes(fill = strain)) + 
#'   theme(legend.position="none") +
#'   # theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
#'   theme(text = element_text(size = 20)) +
#'   labs(title = "Gene 1") + 
#'   xlab("strain") 
#' @export
getCPReadMatrix <- function(vec.num.rep, beta0s, sig2.strains, ps, phis){   

  num.probes <- length(beta0s)
  CountMatrix <- lapply(1:num.probes, function(x){
    probe.x <- getCPReads(vec.num.rep, beta0s[x], 
                          sig2.strains[x], ps[x], phis[x])
  })
  CountMatrix <- do.call(rbind, CountMatrix)
  rownames(CountMatrix) <- paste0("Gene ", 1:num.probes)
  return(CountMatrix)
}



###############################################################################
### Functions to fit GLMM and compute NB VPC 

#' Fit a NBMM for a list of features
#' 
#' Fit a NBMM for a list of features and output the fit parameters.
#' 
#' @param CountMatrix Sequencing count matrix for a list of features. Each row 
#' is for one feature, and the columns are for samples.
#' @param Strains Strain labels for the samples.
#' @param test TRUE or FALSE (default). Test the presence of heritability through examining \eqn{\sigma_a^2 = 0}{}
#' @return A list with two members. The first is a \eqn{G \times 3}{G x 3} matrix indicating the fitted parameters 
#' for each feature The columns are ordered by "beta0", "sigma_a2", "phi". 
#' Row names are feature names; the second member of the list  
#' consists of p-values for testing the hypothesis that \eqn{\sigma_a^2 = 0}{sigma_a2 = 0}.
#' @examples
#' ## Compute vpc for each feature under NBMM. This will take a while on the
#' ##  entire dataset. For the purpose of illustration, here we only fit on 
#' ##  the first 10 features.
#' result.nb <- fitNBMM(simData[1:10, ], strains)
#' @export
fitNBMM <- function(CountMatrix, Strains, test = FALSE){

  GeneIDs <- rownames(CountMatrix)
  paras <- t(pbsapply(1:nrow(CountMatrix), function(x){
    CountVector <- CountMatrix[x, ]
    GeneID <- GeneIDs[x]
    dat_sub <- data.frame(expr = as.numeric(CountVector), strain = Strains)
    para <- tryCatch({
      model_sub2 <- glmmadmb(formula = expr ~ 1 + (1|strain), 
                             data = dat_sub, family = 'nbinom', link = 'log')
      sigma_a2 <- as.numeric(model_sub2$S$strain[1])
      beta0 <- as.numeric(model_sub2$b)
      phi <- 1/model_sub2$alpha
      para_sub <- c(beta0, sigma_a2, phi)
      
      if (test){
        model_sub2_red <- glmmadmb(formula = expr ~ 1, data = dat_sub, 
                                   family = 'nbinom', link = 'log')
        test.stat <- 2*logLik(model_sub2) - 2*logLik(model_sub2_red)
        if (test.stat<1e-6) {test.stat <- 0} 
        # test.stat <- round(test.stat, digits = 1e-6)
        pval <- 0.5*pchisq(test.stat, df = 1, lower.tail = FALSE) + 
          0.5*as.numeric(test.stat == 0)
        para_sub <- c(para_sub, pval)
      }
      
      para_sub
    }, error = function(err){
      print(paste('Using alt method for', GeneID))
      model_sub <- try({glmer.nb(formula = expr ~ 1 + (1|strain), 
                                 data = dat_sub, verbose = F)}, silent=T)
      
      if (class(model_sub) != "try-error"){
        sigma_a2 <- as.numeric(unlist(lme4::VarCorr(model_sub)))
        beta0 <- as.numeric(lme4::fixef(model_sub))
        phi <- 1/getME(model_sub, "glmer.nb.theta")
        para_sub <- c(beta0, sigma_a2, phi)
      }else{
        print(paste("Fitting problem for feature", x,"returning NA"))
        para_sub <- rep(NA, 3)
      }
      
      if (test){
        model_sub_red <- try({glm.nb(formula = expr ~ 1, 
                                     data = dat_sub, link = 'log')}, 
                             silent = TRUE)
        if (class(model_sub) != "try-error" & 
            class(model_sub_red) != "try-error"){
          test.stat <- 2*logLik(model_sub) - 2*logLik(model_sub_red)
          if (test.stat<1e-6) {test.stat <- 0} 
          # test.stat <- round(test.stat, digits = 1e-6)
          pval <- 0.5*pchisq(test.stat, df = 1, lower.tail = FALSE) + 
            0.5*as.numeric(test.stat == 0)
        }else{
          print(paste("Cannot do test for feature", x,"fitting problem."))
          pval <- NA
        }
        para_sub <- c(para_sub, pval)
      }
      
      return(para_sub)
    })
    
    return(para)
  }))
  
  paras1 <- paras[ , 1:3]
  rownames(paras1) <- GeneIDs
  colnames(paras1) <- c("beta0", "sigma_a2", "phi")
  
  if (test){
    return(list(paras = paras1, pvals = paras[ , 4]))
  }else{
    return(list(paras = paras1, pvals = NULL))
  }
  
}






# (INTERNAL)
compute1NBVPC <- function(beta0, sigma_a2, phi){
  # Calculate the negative binomial icc 
  #
  # beta0: intercept.
  # sigma_a2: variance of the random factor.
  # phi: dispersion. 
  #
  # [OUTPUT]
  # vpc: a numerical value for variance partition coefficient computed based
  #   on negative binomial mixture model (NBMM).
  
  vpc <- (exp(sigma_a2) - 1) / 
    (exp(sigma_a2) - 1 + exp(sigma_a2)*phi + exp(-beta0 - sigma_a2 / 2))
  return(vpc)
  
}


#' Calculate the negative binomial vpc for one or more features.
#' 
#' Calculate the negative binomial vpc for one or more features.
#' 
#' @param para A \eqn{1 \times 3}{1x3} vector or \eqn{k \times 3}{k x 3} matrix of negative binomial fit parameters.
#' Order of the parameters: \eqn{({\beta}_0, {\sigma_a^2, \phi)}}{(beta0, sigma_a2, phi)}.
#' @return A \eqn{G \times 1}{G x 1} matrix consisting of variance partition coefficients for
#' \eqn{G} features based on negative binomial mixed model (NBMM). Column name is "NBMM_vpc"; row names are the 
#' feature names.
#' @examples
#' ## Compute vpc for each feature under NBMM.
#' vpc.nb <- computeAllNBVPC(para_nb)
#' 
#' ## Visulize the distribution of the vpcs. 
#' hist(vpc.nb, breaks = 50, col = "cyan")
#' 
#' ## Plot sorted vpcs.
#' plot(sort(vpc.nb), ylab = "Heritability (h2)", ylim = c(0,1), main = "Sorted NB VPC scores")
#' abline(h = 0.9, lty = 2, col = "red")
#' text(50, 0.92, "h2 = 0.9", col = "red")
#' @export
computeAllNBVPC <- function(para){

  if(is.null(dim(para))){
    vpcs <- compute1NBVPC(para[1], para[2], para[3])
  }else{
    vpcs <- apply(para, 1, function(x){
      vpc <- compute1NBVPC(x[1], x[2], x[3])
      return(vpc)
    })
  }
  
  vpcs <- matrix(vpcs, ncol = 1)
  rownames(vpcs) <- rownames(para)
  colnames(vpcs) <- 'NBMM_vpc'
  
  return(vpcs)
}





###############################################################################
### Fit compound Poisson models and estimated VPCs.

#' Fit a compound Poisson mixed effect model for a list of features
#' 
#' Fit a CPMM for a list of features.
#' 
#' @param CountMatrix Sequencing count matrix for a list of features. Each row 
#' is for one feature, and the columns are for samples. 
#' @param Strains Strain labels for the samples. 
#' @param test TRUE or FALSE (default). Test the presence of heritability through examining \eqn{\sigma_a^2 = 0}{}
#' @param optimizer A character string that determines which optimization routine is
#   to be used. Possible choices are "nlminb" (default), "L-BFGS-B", and 
#   "bobyqa".
#' @return A list with two members. The first is a \eqn{G \times 4}{G x 4} matrix 
#' indicating the fitted parameters for each feature The columns are ordered by 
#' "beta0", "p", sigma_a2", "phi". Row names are feature names; the second member of the list  
#' consists of p-values for testing the hypothesis that \eqn{\sigma_a^2 = 0}{sigma_a2 = 0}.
#' @export
fitCPMM <- function(CountMatrix, Strains, test = FALSE, optimizer = "nlminb"){
  # Fit a compound Poisson mixed effect model for a list of probes/genes and 
  #   output the fit parameters.
  #
  # CountMatrix: sequencing count matrix for a list of probes/genes. Each row 
  #   is for one probe/gene, and the columns are for samples. 
  # Strains: strain label for the samples. 
  # test: Logical argument indicating whether to do a test for significance of 
  #   random effects.
  # optimizer: A character string that determines which optimization routine is
  #   to be used. Possible choices are "nlminb" (default), "L-BFGS-B", and 
  #   "bobyqa".
  #
  # [OUTPUT]
  # Return a list with two members. The first is a G by 4 matrix 
  #   indicating the fitted parameters for each gene. The columns are ordered by 
  #   "beta0", "p", sigma_a2", "phi". Row names are gene names; the second member of the list  
  #   consists of p-values for testing the hypothesis that sigma_a2 = 0.
  
  suppressMessages(suppressWarnings(requireNamespace("cplm")))
  
  paras <- t(pbsapply(1:nrow(CountMatrix), function(x){
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
      beta0 <- cplm::fixef(fit)
      sigma_a2 <- as.numeric(cplm::VarCorr(fit)$strain)
      p <- fit$p
      phi <- fit$phi
      
      para <- c(beta0, sigma_a2, p, phi)
    }else{
      print(paste("Fitting problem for feature", x, "returning NA"))
      para <- rep(NA, 4)
    }
    
    ### Fitting the reduced model for testing significance of the random effect ###
    if (test){
        #     detach("package:cplm", unload=TRUE)
      unloadNamespace("cplm")
      suppressMessages(suppressWarnings(requireNamespace("cplm")))
      fit.red <- tryCatch({
        fit2 <- cpglm(expr ~ 1, data = dat_sub, optimizer = optimizer)
      }, error=function(err){
        fit2 <- try({cpglm(expr ~ 1, data = dat_sub, optimizer = optimizer)})
        return(fit2)
      })
      
      if (class(fit.red) != "try-error" & class(fit) != "try-error"){
        test.stat <- 2*logLik(fit)+AIC(fit.red)-2
        #-2*sum(log(dtweedie(dat_sub$expr,fit.red$p,exp(fit.red$coefficients),fit.red$phi)))
        if (test.stat<1e-6) {test.stat <- 0}
        pval <- 0.5*pchisq(test.stat, df = 1, lower.tail = FALSE) + 
          0.5*as.numeric(test.stat == 0)
      }else{
        print(paste("Cannot do test for feature", x, "fitting problem."))
        pval <- NA
      }
      
      para <- c(para, pval)
    }
    
    #    detach("package:cplm", unload=TRUE)
    unloadNamespace("cplm")
    suppressMessages(suppressWarnings(requireNamespace("cplm")))
    return(para)
  }))
  
  paras1 <- paras[,1:4]
  rownames(paras1) <- rownames(CountMatrix)  
  colnames(paras1) <- c("beta0", "sigma_a2", "p", "phi")
  
  if (test){
    return(list(paras = paras1,  pvals = paras[ , 5]))
  }else{
    return(list(paras = paras1,  pvals = NULL))
  }
  
}






# (INTERNAL)
compute1CPVPC <- function(beta0, sigma_a2, p, phi){
  # Calculate the compound Poisson variance partition coefficient.
  #
  # beta0: intercept.
  # sigma_a2: variance of the random factor.
  # p: power index.
  # phi: dispersion
  #
  # [OUTPUT]
  # vpc: a numerical value for variance partition coefficient computed based
  #   on compound Poisson mixture model (CPMM).
  
  vpc.numerator <- exp(2 * beta0 + 2 * sigma_a2 ) - exp(2 * beta0 + sigma_a2)
  vpc <- vpc.numerator/
    (vpc.numerator + phi * exp(p * beta0 + p^2 * sigma_a2 / 2))
  return(vpc)
}


#' Calculate the compound Poisson VPC for one or more features.
#' 
#' Calculate the compound Poisson VPC for one or more features.
#' 
#' @param para A \eqn{1 \times 4}{1 x 4} vector or \eqn{k \times 4}{k x 4} matrix of negative binomial fit parameters.
#' Order of the parameters: \eqn{({\beta}_0, {\sigma_a^2, p, \phi)}}{(beta0, sigma_a2, p, phi)}.
#' @return vpcs A \eqn{G \times 1}{G x 1} matrix consisting of variance partition coefficients for
#'   G features based on compound Poisson mixture model (CPMM) Column name is "CPMM_vpc"; row names are the 
#'   feature names.
#' @examples
#' ## Compute vpc for each feature under CPMM. 
#' vpc.cp <- computeAllCPVPC(para_cp) 
#' 
#' ## Visulize the distribution of the vpcs. 
#' hist(vpc.cp, breaks = 50, col = "cyan")
#' 
#' ## Plot sorted vpcs.
#' plot(sort(vpc.cp), ylab = "Heritability (h2)", ylim = c(0,1), main = "Sorted CP VPC scores")
#' abline(h = 0.9, lty = 2, col = "red")
#' text(50, 0.92, "h2 = 0.9", col = "red")
#' @export
computeAllCPVPC <- function(para){

  if(is.null(dim(para))){
    vpcs <- compute1CPVPC(para[1], para[2], para[3], para[4])
  }else{
    vpcs <- apply(para, 1, function(x){
      vpc <- compute1CPVPC(x[1], x[2], x[3], x[4])
      return(vpc)
    })
  }
  
  vpcs <- matrix(vpcs, ncol = 1)
  rownames(vpcs) <- rownames(para)
  colnames(vpcs) <- 'CPMM_vpc'
  
  return(vpcs)
}




###############################################################################
### Fit linear mixed models and compute VPCs.

# (INTERNAL)
compute1lmerVPC <- function(CountVector, Strains, PriorWeight = NULL, 
                            test = FALSE){
  # Compute the VPC value for one feature
  #
  # CountVector: sequencing counts for the feature.
  # Strains: strain labels for each sample.
  # PriorWeight: weights used in the lmer function.
  #
  # [OUTPUT]
  # Return a list with two members. The first is a numerical value indicating the variance partition 
  #   coefficient (vpc) under a linear mixed model (LMM); the second member of the list is the p-values 
  #   from testing the hypothesis that there is no random effect.
  
  dat_sub <- data.frame(expr = CountVector, strain = Strains)
  model_sub <- lmer(formula = expr ~ 1 + (1|strain), data = dat_sub, 
                    weights = PriorWeight)
  residual_var <- sigma(model_sub)^2
  random_intercept_var <- as.numeric(unlist(lme4::VarCorr(model_sub)))
  vpc <- random_intercept_var/(residual_var + random_intercept_var)
  
  if (test){
    model_sub_red <- lm(formula = expr ~ 1, data = dat_sub, 
                        weights = PriorWeight)
    test.stat <- 2*logLik(model_sub) - 2*logLik(model_sub_red)
    if (test.stat<1e-6) {test.stat <- 0}
    # test.stat <- round(test.stat, digits = 1e-6)
    pval <- 0.5*pchisq(test.stat, df = 1, lower.tail = FALSE) + 
      0.5*as.numeric(test.stat == 0)
    return(list(vpc = vpc, pval = pval))
  }else{
    return(vpc)
  }  
}


#' Compute the VPC values for a list of features under a linear mixed model (LMM).
#' 
#' Compute the VPC values for a list of features under a linear mixed model (LMM).
#' 
#' @param CountMatrix Sequencing count matrix for a list of features. Each row 
#'   is for one feature, and the columns are for samples. 
#' @param Strains Strain labels for the samples. 
#' @param PriorWeights Weights used in the lmer function.
#' @param test TRUE or FALSE (default). Test the presence of heritability through examining \eqn{\sigma_a^2 = 0}{}
#' @return A list with two members. The first is a \eqn{1 \times G}{1 x G} vector indicating the significance or random effects
#'  for each feature under linear mixed model (LMM); the second member of the list consists of the 
#'  p-values from testing the hypothesis that there is no random effect.
#' @examples
#' ## Compute vpc for each feature under LMM.
#' 
#' ## Provide normalized data with prior weights:
#' result.voom <- computeAlllmerVPC(simData_voom, strains, PriorWeights = weights_voom)
#' vpc.voom <- result.voom[[1]]
#' 
#' ## Provide normalized data without prior weights and include hypothesis 
#' ##  testing on presence of heritability:
#' result.vst <- computeAlllmerVPC(simData_vst, strains, test = TRUE)
#' ## Extract parameters
#' vpc.vst <- result.vst[[1]]
#' ## Extract p-values
#' pval.vst <- result.vst[[2]]
#' 
#' ## Visulize the distribution of p-values.
#' hist(pval.vst, breaks = 30, col = "cyan")
#' 
#' ## Compare the vpc estimates. 
#' require(psych)
#' bothvpc <- cbind(vpc.voom, vpc.vst)
#' colnames(bothvpc) <- c("voom", "vst")
#' pairs.panels(bothvpc, ellipses = FALSE, breaks = 30, main = "vpc comparison")
#' @export
computeAlllmerVPC <- function(CountMatrix, Strains, PriorWeights = NULL, 
                              test = FALSE){

  VPC <- pbsapply(1:nrow(CountMatrix), function(x){
    vpc.x <- compute1lmerVPC(CountMatrix[x, ], Strains, PriorWeights[x, ],
                             test = test)
    if (test){
      return(c(vpc.x$vpc, vpc.x$pval))
    }else{
      return(vpc.x)
    }
  })
  VPC = as.matrix(VPC)  
  if(test){
    VPC = t(VPC)
    pvals <- matrix(VPC[,2], ncol = 1)
    rownames(pvals) <- rownames(CountMatrix)
    colnames(pvals) <- 'P-value'
  }
  
  vpcs = as.matrix(VPC[,1], ncol = 1)
  
  rownames(vpcs) <- rownames(CountMatrix)
  colnames(vpcs) <- 'LMM_vpc'
  
  if (test){
    return(list(vpcs=vpcs, pvals=pvals))
  }else{
    return(list(vpcs=vpcs, pvals=NULL))
  }
}






#' SHORT DESCRIPTION OF GETBOOTCI
#' 
#' LONG DESCRIPTION OF GETBOOTCI
#' 
#' @param CountMatrix The data matrix; rows are features, columns are samples
#' @param Strains The vector of strains corresponding to each sample, length = ncol(CountMatrix)
#' @param which.features The vector of select feature numbers for which CI is desired
#' @param num.boot Number of bootstraps
#' @param method Which method should be used, "CP", "NB", or "VST". "VST" method bootstraps data from "NB"
#' @param alpha The CI will be \eqn{100*(1-\alpha)}{100*(1-alpha)} percent CI
#' @param optimizer A character string that determines which optimization routine is
#'  to be used. Possible choices are "nlminb" (default), "L-BFGS-B", and "bobyqa".
#' @return (i) intervals: a matrix of dimension length(which.features) x 2 containing the CIs
#'  (ii) all.vpcs: a matrix of dimension length(which.features) x num.boot containing all
#'  bootstrapped VPC values.
#' @export
GetBootCI = function(CountMatrix, Strains, which.features, num.boot,
                     method="NB", alpha=0.05, optimizer = "nlminb"){
  all.vpcs = matrix(NA, nrow = length(which.features), ncol=num.boot)
  vec.num.rep = as.numeric(table(Strains))
  
  if (method=="NB"){
    print("Getting initial point estimates using NB method")
    fit = fitNBMM(CountMatrix[which.features,], Strains)
    
    for (i in 1:length(which.features)){
      print(paste("Bootstraping feature",i))
      boot.data = getNBReadMatrix(vec.num.rep, 
                                  rep(fit$paras[i,1],num.boot), 
                                  rep(fit$paras[i,2],num.boot),
                                  rep(fit$paras[i,3],num.boot))
      fit.i = fitNBMM(boot.data, Strains)
      all.vpcs[i,] = computeAllNBVPC(fit.i$paras)
    }
    
    intervals = cbind( apply(all.vpcs, 1, quantile, probs = alpha/2),
                       apply(all.vpcs, 1, quantile, probs = 1-alpha/2))
    
    return(list(intervals = intervals, all.vpcs = all.vpcs))
  }
  
  
  
  
  if (method=="CP"){
    print("Getting initial point estimates using CP method")
    fit = fitCPMM(CountMatrix[which.features,], Strains, optimizer = optimizer)
    
    for (i in 1:length(which.features)){
      print(paste("Bootstraping feature",i))
      boot.data = getCPReadMatrix(vec.num.rep, 
                                  rep(fit$paras[i,1],num.boot), 
                                  rep(fit$paras[i,2],num.boot),
                                  rep(fit$paras[i,3],num.boot),
                                  rep(fit$paras[i,4],num.boot))
      fit.i = fitCPMM(boot.data, Strains, optimizer = optimizer)
      all.vpcs[i,] = computeAllCPVPC(fit.i$paras)
    }
    
    intervals = cbind( apply(all.vpcs, 1, quantile, probs = alpha/2),
                       apply(all.vpcs, 1, quantile, probs = 1-alpha/2))
    
    return(list(intervals = intervals, all.vpcs = all.vpcs))
  }
  
  
  
  
  if (method=="VST"){
    print("Getting initial point estimates using NB method")
    fit = fitNBMM(CountMatrix[which.features,], Strains)
    
    for (i in 1:length(which.features)){
      print(paste("Bootstraping feature",i))
      boot.data = getNBReadMatrix(vec.num.rep, 
                                  rep(fit$paras[i,1],num.boot), 
                                  rep(fit$paras[i,2],num.boot),
                                  rep(fit$paras[i,3],num.boot))
      
      cds=newCountDataSet(round(boot.data), Strains) 
      cds=estimateSizeFactors(cds)
      
      if (sum(is.na(sizeFactors(cds)))>0){   # In a few simulations, there are NA's due to too many low counts
        cds=newCountDataSet(round(1+cplm.data), Strains) # 1 added to avoid too many low counts in those cases
        # This adds a small variation to the data
        cds=estimateSizeFactors(cds)
      }
      #sizeFactors(cds)
      cds=estimateDispersions(cds, method="pooled",fitType="local")
      vsd=getVarianceStabilizedData(cds)
      all.vpcs[i,] = computeAlllmerVPC(vsd, Strains)$vpcs
    }
    
    intervals = cbind( apply(all.vpcs, 1, quantile, probs = alpha/2),
                       apply(all.vpcs, 1, quantile, probs = 1-alpha/2))
    
    return(list(intervals = intervals, all.vpcs = all.vpcs))
  }
}



