#' Apply the generalized Bayesian two-stage normal-based causal model with instrumental variables.
#'
#' @description The \code{gts.nnormal} function applies the generalized Bayesian two-stage
#' normal-based causal model to the categorical treatment data.
#' The model best suits the normally-distributed outcome data
#' that are complete or ignorably missing (i.e., missing completely at random or missing at random).
#'
#' @param formula An object of class formula: a symbolic description of the model to be fitted.
#' The details of the model specification are given under "Details".
#' @param data A dataframe with the variables to be used in the model.
#' @param advanced Logical; if FALSE (default), the model is specified using the formula argument,
#' if TRUE, self-defined models can be specified using the adv.model argument.
#' @param adv.model Specify the self-defined model. Used when advanced=TRUE.
#' @param b0 The mean hyperparameter of the normal distribution (prior distribution)
#' for the first-stage generalized causal model coefficients, i.e., coefficients for the instrumental variables.
#' This can either be a numerical value or a vector with dimensions equal to the number of coefficients
#' for the instrumental variables. If this takes a numerical value, then that values will
#' serve as the mean hyperparameter for all of the coefficients for the instrumental variables.
#' Default value of 0 is equivalent to a noninformative prior for the normal distributions.
#' Used when advanced=FALSE.
#' @param B0 The precision hyperparameter of the normal distribution (prior distribution)
#' for the first stage generalized causal model coefficients.
#' This can either be a numerical value or a vector with dimensions equal to the number of coefficients
#' for the instrumental variables. If this takes a numerical value, then that values will
#' serve as the precision hyperparameter for all of the coefficients for the instrumental variables.
#' Default value of 1.0E-6 is equivalent to a noninformative prior for the normal distributions.
#' Used when advanced=FALSE.
#' @param g0 The mean hyperparameter of the normal distribution (prior distribution)
#' for the second-stage generalized causal model coefficients,
#' i.e., coefficients for the treatment variable and other regression covariates).
#' This can either be a numerical value if there is only one treatment variable in the model,
#' or a if there is a treatment variable and multiple regression covariates,
#' with dimensions equal to the total number of coefficients for the treatment variable and covariates.
#' Default value of 0 is equivalent to a noninformative prior for the normal distributions.
#' Used when advanced=FALSE.
#' @param G0 The precision hyperparameter of the normal distribution (prior distribution)
#' for the second-stage generalized causal model coefficients.
#' This can either be a numerical value if there is only one treatment variable in the model,
#' or a vector if there is a treatment variable and multiple regression covariates,
#' with dimensions equal to the total number of coefficients for the treatment variable and covariates.
#' Default value of 1.0E-6 is equivalent to a noninformative prior for the normal distributions.
#' Used when advanced=FALSE.
#' @param e0 The location hyperparameter of the inverse Gamma distribution (prior for the variance of the
#' normal distribution on the model residual).
#' Default of 0.001 is equivalent to the noninformative prior for the inverse Gamma distribution.
#' @param E0 The shape hyperparameter of the inverse Gamma distribution (prior for the variance of the
#' normal distribution on the model residual).
#' Default of 0.001 is equivalent to the noninformative prior for the inverse Gamma distribution.
#' @param beta.start The starting values for the first-stage generalized causal model coefficients,
#' i.e., coefficients for the instrumental variables.
#' This can either be a numerical value or a column vector with dimensions
#' equal to the number of first-stage coefficients.
#' The default value of NA will use the IWLS (iteratively reweighted least squares) estimate
#' of first-stage coefficients as the starting value.
#' If this is a numerical value, that value will
#' serve as the starting value mean for all the first-stage beta coefficients.
#' @param gamma.start The starting values for the second-stage generalized causal model coefficients,
#' i.e., coefficients for the treatment variable and the model covariates.
#' This can either be a numerical value or a column vector with dimensions
#' equal to the number of second-stage coefficients.
#' The default value of NA will use the IWLS (iteratively reweighted least squares) estimate
#' of second-stage coefficients as the starting value.
#' If this is a numerical value, that value will
#' serve as the starting value mean for all the second-stage gamma coefficients.
#' @param e.start The starting value for the precision hyperparameter of the inverse gamma distribution
#' (prior for the scale parameter of Student's t distribution of the model residual).
#' The default value of NA will use the inverse of the residual variance from the
#' IWLS (iteratively reweighted least square) estimate of the second-stage model.
#' @param n.chains The number of Markov chains. The default is 1.
#' @param n.iter The number of total iterations per chain (including burnin). The default is 10000.
#' @param n.burnin Length of burn in, i.e., number of iterations to discard at the beginning.
#' Default is n.iter/2, that is, discarding the first half of the simulations.
#' @param n.thin The thinning rate. Must be a positive integer. The default is 1.
#' @param DIC Logical; if TRUE (default), compute deviance, pD, and DIC. The rule pD=Dbar-Dhat is used.
#' @param codaPkg Logical; if FALSE (default), an object is returned; if TRUE,
#' file names of the output are returned.
#'
#' @return
#' If \emph{codaPkg=FALSE}(default), returns an object containing summary statistics of
#' the saved parameters, including
#' \item{s1.intercept}{Estimate of the intercept from the first stage.}
#' \item{s1.slopeP}{Estimate of the pth slope from the first stage. }
#' \item{s2.intercept}{Estimate of the intercept from the second stage.}
#' \item{s2.slopeP}{Estimate of the pth slope from the second stage (the first slope is always
#' the \strong{LATE}).}
#' \item{var.e.s2}{Estimate of the residual variance at the second stage.}
#' \item{DIC}{Deviance Information Criterion.}
#' If \emph{codaPkg=TRUE}, the returned value is the path for the output file
#' containing the Markov chain Monte Carlo output.
#'
#' @details
#' \enumerate{
#' \item{Bayesian two-stage causal models are specified symbolically.
#' A typical model has the form \emph{reponse ~ terms|instrumental_variables},
#' where response is the (numeric) response vector and terms is a series of terms
#' which specifies a linear predictor (i.e., the treatment variable and the covariates) for the response.
#' The first specification in the term is always the treatment variable and
#' the remaining specifications are always the covariates for the response.}
#' \item{DIC is computed as \emph{mean(deviance)+pD}.}
#' \item{Prior distributions used in ALMOND.}
#' \itemize{
#' \item Generalized causal model coefficients at both stages: normal distributions.
#' \item The generalized causal model residual term: normal distribution.
#' }
#' }
#'
#' @references
#' Gelman, A., Carlin, J.B., Stern, H.S., Rubin, D.B. (2003).
#' \emph{Bayesian data analysis}, 2nd edition. Chapman and Hall/CRC Press.
#'
#' Spiegelhalter, D. J., Thomas, A., Best, N. G., Gilks, W., & Lunn, D. (1996).
#' BUGS: Bayesian inference using Gibbs sampling.
#' \href{http://www.mrc-bsu.cam.ac.uk/bugs}, 19.
#'
#' @examples
#' \donttest{
#' # Run the model
#' model1 <- gts.nnormal(outcome~treatment|instrument,data=simCatNormMCAR)
#'
#' # Run the model with the self-defined advanced feature
# my.model<- function(){
#   for (i in 1:N){
#     logit(p[i]) <- beta0 + beta1*z[i]
#     x[i] ~ dbern(p[i])
#     muY[i] <- gamma0 + gamma1*p[i]
#     y[i] ~ dnorm(muY[i], pre.u2)
#   }
#
#   beta0 ~ dnorm(0,1)
#   beta1 ~ dnorm(1, 1)
#   gamma0 ~ dnorm(0, 1)
#   gamma1 ~ dnorm(.5, 1)
#
#   pre.u2 ~ dgamma(.001, .001)
#
#   s1.intercept <- beta0
#   s1.slope1 <- beta1
#   s2.intercept <- gamma0
#   s2.slope1 <- gamma1
#   var.e.s2 <- 1/pre.u2
# }
#'
#' model2 <- gts.nnormal(outcome~treatment|instrument,data=simCatNormMCAR,
#' advanced=TRUE, adv.model=my.model)
#'
#' # Extract the model DIC
#' model1$DIC
#'
#' # Extract the MCMC output
#' model3 <- gts.nnormal(outcome~treatment|instrument,data=simCatNormMCAR,codaPkg=TRUE)
#' }
#'
#' @export
gts.normal<-function(formula,data,advanced=FALSE, adv.model,
                     b0=1,B0=1.0E-6, g0=0, G0=1.0E-6, e0=0.001,E0=0.001,
                     beta.start=NULL, gamma.start=NULL, e.start=NULL,
                     n.chains=1,n.burnin=floor(n.iter/2),n.iter=10000,n.thin=1,DIC,debug=FALSE,
                     codaPkg=FALSE){
  .alReplaceSciNotR <- function(x,digits=5){
    x[abs(x)<1e-3||abs(x)>1e+4] <-
      formatC(x[abs(x)<1e-3||abs(x)>1e+4],
              digits=digits,
              format="E")
    return(x)
  }

  formula1 <- formula(formula)
  y <- model.response(model.frame(as.Formula(formula1),data=data,na.action=NULL))
  x <- as.matrix(model.frame(as.Formula(formula1),data=data,na.action=NULL,rhs=1)[,-1])
  z <- as.matrix(model.frame(as.Formula(formula1),data=data,na.action=NULL,rhs=2)[,-1])
  N <- length(y)
  if (length(levels(factor(x)))!=2){
    stop('The factor level of the treatment variable should be 2')
  }else{
    xList <- lapply(1:ncol(x),function(i){x[,i]})
    zList <- lapply(1:ncol(z),function(i){z[,i]})

    if(length(b0)==1){
      b0 <- rep(b0,length(zList)+1)
    } else if(length(b0)>(length(zList)+1)){
      warning(paste0("Number of priors is greater than number of parameters. Only first ",
                     length(zList)+1," values of b0 will be used."))
      b0 <- b0[1:(length(zList)+1)]
    }

    if(length(B0)==1){
      B0 <- rep(B0,length(zList)+1)
    } else if(length(B0)>(length(zList)+1)){
      warning(paste0("Number of priors is greater than number of parameters. Only first ",
                     length(zList)+1," values of B0 will be used."))
      B0 <- B0[1:(length(zList)+1)]
    }

    if(length(g0)==1){
      g0 <- rep(g0,length(xList)+1)
    } else if(length(g0)>(length(xList)+1)){
      warning(paste0("Number of priors is greater than number of parameters. Only first ",
                     length(xList)+1," values of g0 will be used."))
      g0 <- g0[1:(length(xList)+1)]
    }

    if(length(G0)==1){
      G0 <- rep(G0,length(xList)+1)
    } else if(length(G0)>(length(xList)+1)){
      warning(paste0("Number of priors is greater than number of parameters. Only first ",
                     length(xList)+1," values of G0 will be used."))
      G0 <- G0[1:(length(xList)+1)]
    }

    if (length(xList)==1){
      names(xList) <- "x"
    }else{
      names(xList) <- c("x",paste0("x",1:(length(xList)-1)))
    }
    if (length(zList)==1){
      names(zList) <- "z"
    }else{
      names(zList) <- paste0("z",1:length(zList))
    }
    data = c(list("N"=N,"y"=y),xList,zList)

    lm.data = as.data.frame(cbind(y,x,z))
    if (length(xList)==1){
      xnames = "x"
    }else{
      xnames = c("x",paste0("x",1:(length(xList)-1)))
    }
    if (length(zList)==1){
      znames = "z"
    }else{
      znames = paste0("z",1:length(zList))
    }
    colnames(lm.data) = c("y",xnames,znames)

    if (length(zList)==1){
      coefList1s = as.list(coefficients(summary(glm(x~z,data=data,family='binomial',na.action=na.omit)))[,1])
      names(coefList1s) = c("beta0","beta1")
      xpred <- predict(glm(x~z,data=data,family='binomial',na.action=na.exclude))
    }else{
      coefList1s = as.list(coefficients(summary(lm(formula(paste0("x~",
                                                                  paste0("z",1:length(zList),collapse="+"))),data=data,na.action=na.omit)))[,1])
      names(coefList1s) <- paste0("beta",0:(length(coefList1s)-1))
      xpred <- predict(lm(formula(paste0("x~",paste0("z",1:length(zList),collapse='+'))),data=data,
                          na.action=na.exclude))
    }

    if (length(xList)==1){
      coefList2s = as.list(coefficients(summary(lm(y~xpred,na.action=na.omit)))[,1])
      names(coefList2s) <- c("gamma0","gamma1")
      pre.u2 = 1/(summary(lm(y~xpred,na.action=na.omit))$sigma)
    }else{
      coefList2s = as.list(coefficients(summary(lm(formula(paste0("y~",
                                                                  paste0("xpred+",paste0("x",1:(length(xList)-1),collapse="+")))),
                                                   data=data,na.action=na.omit)))[,1])
      names(coefList2s) <- paste0("gamma",0:(length(coefList2s)-1))
      pre.u2 = 1/summary(lm(formula(paste0("y~",
                                           paste0("xpred+",paste0("x",1:(length(xList)-1),collapse="+")))),
                            data=data,na.action=na.omit))$sigma
    }
    names(pre.u2) = c("pre.u2")

    if (is.null(beta.start)==TRUE){
      beta.start = coefList1s
    } else if (length(beta.start) < length(coefList1s)){
      warning("Number of starting values is fewer than the number of parameters. The first
              element of beta.start ",beta.start[1], " will be used for all betas." )
      beta.start = rep(list(beta.start[1]),length(coefList1s))
    } else if (length(beta.start) > length(coefList1s)){
      warning(paste0("Number of starting values is greater than the number of parameters. Only the first ",
                     length(coefList1s)," values of betas will be used."))
      beta.start = as.list(beta.start[1:(length(coefList1s))])
    } else{
      beta.start = as.list(beta.start)
    }
    names(beta.start) = paste0("beta",0:(length(beta.start)-1))

    if (is.null(gamma.start)==TRUE){
      gamma.start = coefList2s
    } else if (length(gamma.start) < length(coefList2s)){
      warning("Number of starting values is fewer than the number of parameters. The first
              element ",gamma.start[1], "will be used for all gammas." )
      gamma.start = rep(list(gamma.start[1]),length(coefList2s))
    } else if (length(gamma.start) > length(coefList2s)){
      warning(paste0("Number of starting values is greater than the number of parameters. Only the first ",
                     length(coefList2s)," values of gammas will be used."))
      gamma.start = as.list(gamma.start[1:(length(coefList2s))])
    } else {
      gamma.start = as.list(gamma.start)
    }
    names(gamma.start) <- paste0("gamma",0:(length(gamma.start)-1))

    if (is.null(e.start)==TRUE){
      e.start = pre.u2
    } else if (length(e.start)>1){
      warning(paste0("Number of starting values has length > 1.
                     Only the first element will be used."))
      e.start = e.start[1]
    } else {
      e.start = e.start
    }
    names(e.start) = c("pre.u2")

    if (advanced==FALSE){
      L1 <- paste0("model\n{\n\tfor (i in 1:N){","\n")
      if (length(zList)==1){
        L2 <- paste0("\t\tlogit(p[i]) <- beta0 + beta1*z[i]")
      }else{
        L2 <- paste0("\t\tlogit(p[i]) <- beta0+",paste0("beta",1:(length(zList)),"*z",1:(length(zList)),
                                                        "[i]",collapse="+"),"\n")
      }
      L3 <- paste0("\t\tx[i] ~ dbern(p[i])","\n")
      if (length(xList)==1){
        L4 <- paste0("\t\tmuY[i] <- gamma0+gamma1*p[i]")
      }else{
        L4 <- paste0("\t\tmuY[i] <- gamma0+gamma1*p[i]+",paste0("gamma",2:(length(xList)),"*x",
                                                                1:(length(xList)-1),"[i]",collapse="+"),"\n")
      }
      L5 <- paste0("\t\ty[i] ~ dnorm(muY[i], pre.u2)\n\t}\n")
      LFormulasBeta <- do.call(paste0,lapply(1:(length(zList)+1),function(i){
        paste0("\tbeta",i-1," ~ dnorm(",.alReplaceSciNotR(b0[i]),",",
               .alReplaceSciNotR(B0[i]),")","\n")
      }))

      LFormulasGamma <- do.call(paste0,lapply(1:(length(xList)+1),function(i){
        paste0("\tgamma",i-1," ~ dnorm(",.alReplaceSciNotR(g0[i]),
               ",",.alReplaceSciNotR(G0[i]),")","\n")
      }))
      LN <- paste0("\tpre.u2 ~ dgamma(",.alReplaceSciNotR(e0),
                   ",",.alReplaceSciNotR(E0),")\n\n")


      tempParNamesOrig <- c(paste0("beta",0:length(zList)),
                            paste0("gamma",0:(length(xList))),
                            "1/pre.u2")
      tempParNamesNew <- c("s1.intercept",paste0("s1.slope",1:length(zList)),
                           "s2.intercept",
                           paste0("s2.slope",1:length(xList)),
                           "var.e.s2")
      LPara <- do.call(paste0,lapply(1:(length(tempParNamesOrig)),
                                     function(j){
                                       paste0("\t",tempParNamesNew[j]," <- ",tempParNamesOrig[j],"\n")
                                     }
      ))

      tmpModel <- paste0(L1,L2,L3,L4,L5,LFormulasBeta,LFormulasGamma,LN,LPara,"}\n")
    }else{
      tmpModel <- capture.output(adv.model)
      tmpModel[1] <- "model\n{"
      tmpModel <- paste(tmpModel,collapse="\n")
      tempParNamesNew <- c("s1.intercept",paste0("s1.slope",1:length(zList)),
                           "s2.intercept",
                           paste0("s2.slope",1:length(xList)),
                           "var.e.s2")
    }
    tempFileName <- tempfile("model")
    tempFileName <- paste(tempFileName, "txt", sep = ".")
    writeLines(tmpModel,con = tempFileName)
    modelLoc <- gsub("\\\\", "/", tempFileName)

    inits <- function(){
      c(beta.start,gamma.start,e.start)
    }

    parameters<- tempParNamesNew

    output<-bugs(data,inits,parameters,modelLoc,
                 n.chains=as.vector(n.chains),n.thin=as.vector(n.thin),
                 n.burnin=as.integer(n.burnin),n.iter=as.integer(n.iter),
                 DIC=TRUE,debug=as.logical(debug),codaPkg=as.logical(codaPkg))
    return(output)
  }
}


