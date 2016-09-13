#####cpglmm utility functions #####
#install.packages("statmod")
requireNamespace("statmod", quietly = TRUE)


cpglm.init <- function(fr, link.power = 0){
  p <- 1.5
  # print(p)
  fit <- cpglm.fit(fr, p, link.power)
  # print(p)
  beta <- as.numeric(fit$coefficients)
  # print(p)
  phi <- sum(fit$weights * fit$residuals^2) / fit$df.residual
  # print(p)
  vbeta <- summary.glm(fit)$cov.scaled
  # print(p)
  list(beta = beta, phi = phi, p = p, vcov = vbeta)
}


cpglm.fit <- function(fr, p = 1.5, link.power = 0) {
  fm <- tweedie(var.power = p, link.power = link.power)
  int <- attr(attr(fr$mf,"terms"), "intercept") > 0L
  suppressWarnings(glm.fit(fr$X, fr$Y, weights = fr$wts, offset = fr$off,
                           family = fm, intercept = int))
}

expand.call <-
  function(call = sys.call(sys.parent(1)))
  {
    #given args:
    ans <- as.list(call)
    # ans1 <- ans[[1]]
    # ans <- lapply(ans[-1], eval, envir = sys.frame(sys.parent(2)))
    # ans <- c(ans1, ans)
    
    #possible args:
    frmls <- formals(safeDeparse(ans[[1]]))
    #remove formal args with no presets:
    frmls <- frmls[!sapply(frmls, is.symbol)]
    
    add <- which(!(names(frmls) %in% names(ans)))
    return(as.call(c(ans, frmls[add])))
  }


make.link.power <- function(link) {
  if (!is.character(link) && !is.numeric(link))
    stop("link.power must be either numeric or character.")
  if (is.character(link)){  
    okLinks <- c("log", "identity", "sqrt","inverse")
    if (link %in% okLinks) 
      switch(link, log = 0, identity = 1, sqrt = 0.5, inverse = -1) else
        stop("invalid link function!")
  } else 
    link  
}


safeDeparse <- function(expr){
  ret <- paste(deparse(expr), collapse="")
  #rm whitespace
  gsub("[[:space:]][[:space:]]+", " ", ret)
}



findbars <- function(term)
  ### Return the pairs of expressions that separated by vertical bars
{
  if (is.name(term) || !is.language(term)) return(NULL)
  if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
  if (!is.call(term)) stop("term must be of class call")
  if (term[[1]] == as.name('|')) return(term)
  if (length(term) == 2) return(findbars(term[[2]]))
  c(findbars(term[[2]]), findbars(term[[3]]))
}





subbars <- function(term)
  ### Substitute the '+' function for the '|' function
{
  if (is.name(term) || !is.language(term)) return(term)
  if (length(term) == 2) {
    term[[2]] <- subbars(term[[2]])
    return(term)
  }
  stopifnot(length(term) >= 3)
  if (is.call(term) && term[[1]] == as.name('|'))
    term[[1]] <- as.name('+')
  for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
  term
}



nobars <- function(term)
  ### Return the formula omitting the pairs of expressions that are
  ### separated by vertical bars
{
  if (!('|' %in% all.names(term))) return(term)
  if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
  if (length(term) == 2) {
    nb <- nobars(term[[2]])
    if (is.null(nb)) return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- nobars(term[[2]])
  nb3 <- nobars(term[[3]])
  if (is.null(nb2)) return(nb3)
  if (is.null(nb3)) return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}


expandSlash <- function(bb)
  ### expand any slashes in the grouping factors returned by findbars
{
  if (!is.list(bb)) return(expandSlash(list(bb)))
  ## I really do mean lapply(unlist(... - unlist returns a
  ## flattened list in this case
  unlist(lapply(bb, function(x) {
    if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
      return(lapply(unlist(makeInteraction(trms)),
                    function(trm) substitute(foo|bar,
                                             list(foo = x[[2]],
                                                  bar = trm))))
    x
  }))
}


dimsNames <- c("nt", "n", "p", "q", "s", "np", "LMM", "REML",
               "fTyp", "lTyp", "vTyp", "nest", "useSc", "nAGQ",
               "verb", "mxit", "mxfn", "cvg")
dimsDefault <- list(s = 1L,             # identity mechanistic model
                    mxit= 300L,         # maximum number of iterations
                    mxfn= 900L, # maximum number of function evaluations
                    verb= 0L,           # no verbose output
                    np= 0L,             # number of parameters in ST
                    LMM= 0L,            # not a linear mixed model
                    REML= 0L,         # glmer and nlmer don't use REML
                    fTyp= 2L,           # default family is "gaussian"
                    lTyp= 5L,           # default link is "identity"
                    vTyp= 1L, # default variance function is "constant"
                    useSc= 1L, # default is to use the scale parameter
                    nAGQ= 1L,                  # default is Laplace
                    cvg = 0L)                  # no optimization yet attempted

devNames <- c("ML", "REML", "ldL2", "ldRX2", "sigmaML",
              "sigmaREML", "pwrss", "disc", "usqr", "wrss",
              "dev", "llik", "NULLdev")



lmerFactorList <- function(formula, fr, rmInt, drop)
{
  mf <- fr$mf
  ## record dimensions and algorithm settings
  
  ## create factor list for the random effects
  bars <- expandSlash(findbars(formula[[3]]))
  if (!length(bars)) stop("No random effects terms specified in formula")
  names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
  fl <- lapply(bars,
               function(x)
               {
                 ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
                                       list(fac = x[[3]])), mf)
                 im <- as(ff, "sparseMatrix") # transpose of indicators
                 ## Could well be that we should rather check earlier .. :
                 if(!isTRUE(validObject(im, test=TRUE)))
                   stop("invalid conditioning factor in random effect: ", format(x[[3]]))
                 
                 mm <- model.matrix(eval(substitute(~ expr, # model matrix
                                                    list(expr = x[[2]]))),
                                    mf)
                 if (rmInt) {
                   if (is.na(icol <- match("(Intercept)", colnames(mm)))) break
                   if (ncol(mm) < 2)
                     stop("lhs of a random-effects term cannot be an intercept only")
                   mm <- mm[ , -icol , drop = FALSE]
                 }
                 ans <- list(f = ff,
                             A = do.call(rBind,
                                         lapply(seq_len(ncol(mm)), function(j) im)),
                             Zt = do.call(rBind,
                                          lapply(seq_len(ncol(mm)),
                                                 function(j) {im@x <- mm[,j]; im})),
                             ST = matrix(0, ncol(mm), ncol(mm),
                                         dimnames = list(colnames(mm), colnames(mm))))
                 if (drop) {
                   ## This is only used for nlmer models.
                   ## Need to do something more complicated for A
                   ## here.  Essentially you need to create a copy
                   ## of im for each column of mm, im@x <- mm[,j],
                   ## create the appropriate number of copies,
                   ## prepend matrices of zeros, then rBind and drop0.
                   ans$A@x <- rep(0, length(ans$A@x))
                   ans$Zt <- drop0(ans$Zt)
                 }
                 ans
               })
  dd <-
    VecFromNames(dimsNames, "integer",
                 c(list(n = nrow(mf), p = ncol(fr$X), nt = length(fl),
                        q = sum(sapply(fl, function(el) nrow(el$Zt)))),
                   dimsDefault))
  ## order terms by decreasing number of levels in the factor but don't
  ## change the order if this is already true
  nlev <- sapply(fl, function(el) length(levels(el$f)))
  ## determine the number of random effects at this point
  if (any(diff(nlev)) > 0) fl <- fl[rev(order(nlev))]
  ## separate the terms from the factor list
  trms <- lapply(fl, "[", -1)
  names(trms) <- NULL
  fl <- lapply(fl, "[[", "f")
  attr(fl, "assign") <- seq_along(fl)
  ## check for repeated factors
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    ## check that the lengths of the number of levels coincide
    fl <- fl[match(ufn, fnms)]
    attr(fl, "assign") <- match(fnms, ufn)
  }
  names(fl) <- ufn
  ## check for nesting of factors
  dd["nest"] <- all(sapply(seq_along(fl)[-1],
                           function(i) isNested(fl[[i-1]], fl[[i]])))
  list(trms = trms, fl = fl, dims = dd)
}




### FIXME: somehow the environment of the mf formula does not have
### .globalEnv in its parent list.  example(Mmmec, package = "mlmRev")
### used to have a formula of ~ offset(log(expected)) + ... and the
### offset function was not found in eval(mf, parent.frame(2))
lmerFrames <- function(mc, formula, contrasts, vnms = character(0))
  ### Create the model frame, X, Y, wts, offset and terms
  
  ### mc - matched call of calling function
  ### formula - two-sided formula
  ### contrasts - contrasts argument
  ### vnms - names of variables to be included in the model frame
{
  mf <- mc
  m <- match(c("data", "subset", "weights", "na.action", "offset"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  
  ## The model formula for evaluation of the model frame.  It looks
  ## like a linear model formula but includes any random effects
  ## terms and any names of parameters used in a nonlinear mixed model.
  frame.form <- subbars(formula)      # substitute `+' for `|'
  if (length(vnms) > 0)               # add the variables names for nlmer
    frame.form[[3]] <-
    substitute(foo + bar,
               list(foo = parse(text = paste(vnms, collapse = ' + '))[[1]],
                    bar = frame.form[[3]]))
  
  ## The model formula for the fixed-effects terms only.
  fixed.form <- nobars(formula)       # remove any terms with `|'
  if (!inherits(fixed.form, "formula"))
    ## RHS is empty - use `y ~ 1'
    fixed.form <- as.formula(substitute(foo ~ 1, list(foo = fixed.form)))
  
  ## attach the correct environment
  environment(fixed.form) <- environment(frame.form) <- environment(formula)
  
  ## evaluate a model frame
  mf$formula <- frame.form
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  fe <- mf                            # save a copy of the call
  mf <- eval(mf, parent.frame(2))
  
  ## evaluate the terms for the fixed-effects only (used in anova)
  fe$formula <- fixed.form
  fe <- eval(fe, parent.frame(2)) # allow model.frame to update them
  
  ## response vector
  Y <- model.response(mf, "any")
  ## avoid problems with 1D arrays, but keep names
  if(length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  mt <- attr(fe, "terms")
  
  ## Extract X checking for a null model. This check shouldn't be
  ## needed because an empty formula is changed to ~ 1 but it can't hurt.
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts) else matrix(,NROW(Y),0)
  storage.mode(X) <- "double"      # when ncol(X) == 0, X is logical
  fixef <- numeric(ncol(X))
  names(fixef) <- colnames(X)
  dimnames(X) <- NULL
  
  ## Extract the weights and offset.  For S4 classes we want the
  ## `not used' condition to be numeric(0) instead of NULL
  wts <- model.weights(mf); if (is.null(wts)) wts <- numeric(0)
  off <- model.offset(mf); if (is.null(off)) off <- numeric(0)
  
  ## check weights and offset
  if (any(wts <= 0))
    stop(gettextf("negative weights or weights of zero are not allowed"))
  if(length(off) && length(off) != NROW(Y))
    stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                  length(off), NROW(Y)))
  
  ## remove the terms attribute from mf
  attr(mf, "terms") <- mt
  list(Y = Y, X = X, wts = as.double(wts), off = as.double(off), mf = mf, fixef = fixef)
}


slashTerms <- function(x)
  ### Return the list of '/'-separated terms in an expression that
  ### contains slashes
{
  if (!("/" %in% all.names(x))) return(x)
  if (x[[1]] != as.name("/"))
    stop("unparseable formula for grouping factor")
  list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}


cplm.control <- function(max.iter = 300L,
                         max.fun = 2000L,               
                         bound.p = c(1.01, 1.99),
                         trace = 0,
                         PQL.init = TRUE){         
  if (!is.numeric(max.iter) || max.iter <= 0) 
    stop("value of 'max.iter' must be > 0")
  if (!is.numeric(max.fun) || max.fun <= 0) 
    stop("value of 'max.fun' must be > 0")
  if (!is.numeric(bound.p) || length(bound.p) != 2)
    stop("'bound.p' must be of length 2")
  if (min(bound.p) < 1 || max(bound.p) > 2)
    stop("invalid bounds in 'bound.p'")          
  if (!is.numeric(trace) && !is.logical(trace))
    stop("'trace' must be logical or numeric")
  
  list(max.iter = as.integer(max.iter),
       max.fun = as.integer(max.fun),
       bound.p = as.numeric(sort(bound.p)),
       trace = as.integer(trace),
       PQL.init = as.logical(PQL.init))
}

VecFromNames <- function(nms, mode = "numeric", defaults = list())
{
  ans <- vector(mode = mode, length = length(nms))
  names(ans) <- nms
  ans[] <- NA
  if ((nd <- length(defaults <- as.list(defaults))) > 0) {
    if (length(dnms <- names(defaults)) < nd)
      stop("defaults must be a named list")
    stopifnot(all(dnms %in% nms))
    ans[dnms] <- as(unlist(defaults), mode)
  }
  ans
}


mkZt <- function(FL, start, s = 1L)
  ### Create the standard versions of flist, Zt, Gp, ST, A, Cm,
  ### Cx, and L. Update dd.
{
  dd <- FL$dims
  fl <- FL$fl
  asgn <- attr(fl, "assign")
  trms <- FL$trms
  ST <- lapply(trms, `[[`, "ST")
  Ztl <- lapply(trms, `[[`, "Zt")
  Zt <- do.call(rBind, Ztl)
  Zt@Dimnames <- vector("list", 2)
  Gp <- c(0L, cumsum(vapply(Ztl, nrow, 1L, USE.NAMES=FALSE)))
  .Call("mer_ST_initialize", ST, Gp, Zt)
  A <- do.call(rBind, lapply(trms, `[[`, "A"))
  rm(Ztl, FL)                         # because they could be large
  nc <- sapply(ST, ncol)         # of columns in els of ST
  Cm <- createCm(A, s)
  L <- .Call("mer_create_L", Cm)
  if (s < 2) Cm <- new("dgCMatrix")
  if (!is.null(start) && checkSTform(ST, start)) ST <- start
  
  nvc <- sapply(nc, function (qi) (qi * (qi + 1))/2) # no. of var. comp.
  ### FIXME: Check number of variance components versus number of
  ### levels in the factor for each term. Warn or stop as appropriate
  
  dd["np"] <- as.integer(sum(nvc))    # number of parameters in optimization
  dev <- VecFromNames(devNames, "numeric")
  fl <- do.call(data.frame, c(fl, check.names = FALSE))
  attr(fl, "assign") <- asgn
  
  list(Gp = Gp, ST = ST, A = A, Cm = Cm, L = L, Zt = Zt,
       dd = dd, dev = dev, flist = fl)
}

famNms <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
            "poisson")
linkNms <- c("logit", "probit", "cauchit", "cloglog", "identity",
             "log", "sqrt", "1/mu^2", "inverse")
varNms <- c("constant", "mu(1-mu)", "mu", "mu^2", "mu^3")

famType <- function(family)
{
  if (!(fTyp <- match(family$family, famNms, nomatch = 0)))
    stop(gettextf("unknown GLM family: %s",
                  sQuote(family$family), domain = "R-lme4"))
  if (!(lTyp <- match(family$link, linkNms, nomatch = 0)))
    stop(gettextf("unknown link: %s",
                  sQuote(family$link), domain = "R-lme4"))
  vNam <- switch(fTyp,
                 "mu(1-mu)",          # binomial
                 "constant",          # gaussian
                 "mu^2",              # Gamma
                 "mu^3",              # inverse.gaussian
                 "mu")                # poisson
  if (!(vTyp <- match(vNam, varNms, nomatch = 0)))
    stop(gettextf("unknown GLM family: %s",
                  sQuote(family$family), domain = "R-lme4"))
  c(fTyp = fTyp, lTyp = lTyp, vTyp = vTyp)
}



createCm <- function(A, s)
  ### Create the nonzero pattern for the sparse matrix Cm from A.
  ### ncol(A) is s * ncol(Cm).  The s groups of ncol(Cm) consecutive
  ### columns in A are overlaid to produce Cm.
{
  stopifnot(is(A, "dgCMatrix"))
  s <- as.integer(s)[1]
  if (s == 1L) return(A)
  if ((nc <- ncol(A)) %% s)
    stop(gettextf("ncol(A) = %d is not a multiple of s = %d",
                  nc, s))
  ncC <- as.integer(nc / s)
  TA <- as(A, "TsparseMatrix")
  as(new("dgTMatrix", Dim = c(nrow(A), ncC),
         i = TA@i, j = as.integer(TA@j %% ncC), x = TA@x),
     "CsparseMatrix")
}








############# Redefine cpglmm ##############



mycpglmm <- function(formula, link = "log", data, weights, offset, subset, 
                 na.action, inits = NULL, contrasts = NULL, control = list(), 
                 basisGenerators = c("tp", "bsp", "sp2d"), optimizer = "nlminb", 
                 doFit = TRUE, nAGQ = 1) 
{
  call <- expand.call(match.call())
  if (missing(data)) 
    data <- environment(formula)
  link.power <- make.link.power(link)
  formula <- eval(call$formula)
  tf <- terms.formula(formula, specials = eval(call$basisGenerators, 
                                               parent.frame(2)))
  n.f <- length(unlist(attr(tf, "specials")))
  if (n.f) {
    call2 <- as.list(call)[-1]
    call2 <- call2[-match(c("link", "inits", "control", "optimizer", 
                            "doFit", "nAGQ"), names(call2), 0L)]
    setup <- do.call(frFL, as.list(call2))
    fr <- setup$m$fr
    FL <- setup$m$FL
  }
  else {
    fr <- lmerFrames(call, formula, contrasts)
    FL <- lmerFactorList(formula, fr, 0L, 0L)
  }
  ctr <- do.call(cplm.control, control)
  FL$dims["mxit"] <- ctr$max.iter
  FL$dims["mxfn"] <- ctr$max.fun
  dm <- mkZt(FL, NULL)
  dm$dd["verb"] <- ctr$trace
  if ((nAGQ <- as.integer(nAGQ)) < 1) 
    nAGQ <- 1L
  if (nAGQ%%2 == 0) 
    nAGQ <- nAGQ + 1L
  dm$dd["nAGQ"] <- as.integer(nAGQ)
  AGQlist <- .Call("cpglmm_ghq", nAGQ)
  M1 <- length(levels(dm$flist[[1]]))
  n <- ncol(dm$Zt)
  q <- dm$dd[["q"]]
  d <- dm$dd[["p"]]
  if (is.null(fr$wts) || length(fr$wts) == 0) 
    fr$wts <- as.double(rep(1, n))
  if (is.null(fr$off) || length(fr$off) == 0) 
    fr$off <- as.double(rep(0, n))
  if (M1 >= n) {
    msg1 <- "Number of levels of a grouping factor for the random effects\n"
    msg3 <- "n, the number of observations"
    if (dm$dd["useSc"]) 
      stop(msg1, "must be less than ", msg3)
    else if (M1 == n) 
      message(msg1, "is *equal* to ", msg3)
  }
  if (!is.null(inits)) {
    check.inits.cpglmm(inits, dm$dd["p"], dm$dd["nt"])
  }
  else {
    inits <- cpglm.init(fr, link.power)
  }
  names(inits$beta) <- names(fr$fixef)
  ans <- new(Class = "cpglmm", env = new.env( ), nlmodel = (~I(x))[[2]], 
             frame = fr$mf, call = call, flist = dm$flist, 
             Zt = dm$Zt, X = fr$X, y = as.numeric(fr$Y), 
             pWt = fr$wts, offset = fr$off,
             Gp = unname(dm$Gp), dims = dm$dd, 
             ST = dm$ST, A = dm$A, Cm = dm$Cm, L = dm$L, 
             Cx = rep(1.0, length((dm$A)@x)),  
             # was Cx = dm$A)@x which broke in R3.0.x. 
             # Referecen to the same value. Need duplicate in C level
             deviance = dm$dev, 
             fixef = as.numeric(inits$beta * 1), 
             ranef = numeric(q), 
             u = numeric(q), eta = numeric(n), 
             mu = numeric(n), resid = numeric(n), 
             muEta = numeric(n), var = numeric(n),
             sqrtXWt = as.matrix(numeric(n)), sqrtrWt = numeric(n), 
             RZX = matrix(0, q, d), RX = matrix(0, d, d), 
             ghx = AGQlist[[1]], ghw = AGQlist[[2]], 
             p = as.double(inits$p * 1),      # force copying object
             phi = as.double(inits$phi * 1), 
             link.power = as.double(link.power * 1), 
             bound.p = ctr$bound.p, formula = formula, contrasts = contrasts,
             model.frame = fr$mf, inits = inits, vcov = matrix(0, d, d), smooths = list())
  if (!doFit) 
    return(ans)
  if (optimizer == "nlminb") {
    invisible(.Call("cpglmm_optimize", ans))
    if (ans@dims[["cvg"]] > 6) 
      warning(convergenceMessage(ans@dims[["cvg"]]))
  }
  else {
    cpglmm_dev <- function(parm) {
      .Call("cpglmm_update_dev", ans, parm)
    }
    parm <- c(.Call("cpglmm_ST_getPars", ans), ans$fixef, 
              log(ans$phi), ans$p)
    parm <- unname(parm)
    n.parm <- length(parm)
    lower <- rep(-Inf, n.parm)
    upper <- rep(Inf, n.parm)
    lower[1:sum(sapply(ans$ST, ncol))] <- 0
    lower[n.parm] <- ans$bound.p[1]
    upper[n.parm] <- ans$bound.p[2]
    rslt <- cplm_optim(parm, cpglmm_dev, lower = lower, upper = upper, 
                       control = ctr, optimizer = optimizer)
    ans@dims[["cvg"]] <- as.integer(rslt$convergence)
    if (rslt$convergence) 
      warning(rslt$message)
    invisible(.Call("cpglmm_update_dev", ans, rslt$par))
  }
  invisible(.Call("cpglmm_update_ranef", ans))
  dev <- ans@deviance
  dev["sigmaML"] <- sqrt(ans@phi)
  ans@deviance <- dev
  invisible(.Call("cpglmm_update_RX", ans))
  ans@vcov <- vcov(ans)
  if (n.f) {
    ans@smooths <- indsF(ans, setup$fct, setup$fctterm)
    vars <- unname(sapply(ans@smooths, function(tt) {
      pos <- grep("^x\\d?$", names(tt), perl = TRUE)
      sapply(pos, function(x) as.character(tt[[x]]))
    }))
    vars <- as.character(vars)
    ans@frame <- data.frame(ans@frame, base::subset(data, 
                                                    select = vars))
  }
  return(ans)
}
