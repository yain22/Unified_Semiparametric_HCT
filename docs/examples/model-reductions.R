# The population example in the model explorer (balanced covariates and constant drift).
# No fitted observations or add-on packages are required. Each row has its exact probability.
model_example <- function(model = c("unrestricted","transportable","rct","ect"),
                          delta = .2, external_variance = 1) {
  model <- match.arg(model)
  stopifnot(length(delta)==1L,is.finite(delta),length(external_variance)==1L,
            is.finite(external_variance),external_variance>0)
  groups <- switch(model, rct=c("treated","internal"), ect=c("treated","external"),
                   c("treated","internal","external"))
  q <- switch(model,rct=1,ect=.4,.5)
  e <- if (model=="ect") 1 else 2/3
  if (model=="transportable") delta <- 0
  lambda <- switch(model,rct=1,ect=0,q*(1-e)/(q*(1-e)+(1-q)/external_variance))
  obs <- expand.grid(x=0:1, group=groups, atom=-1:1,stringsAsFactors=FALSE)
  with(obs,{
    R <- group!="external"; A <- group=="treated"; I <- group=="internal"
    weight <- .5 * ifelse(A,q*e,ifelse(I,q*(1-e),1-q)) * ifelse(atom==0,2/3,1/6)
    residual <- atom*sqrt(3*ifelse(R,1,external_variance))
    # The trial covariate contrast is constant at these settings, so its centered
    # standardization contribution is zero for all records.
    score <- numeric(length(x))
    score[A] <- residual[A]/(q*e)
    if (any(I)) score[I] <- -lambda*residual[I]/(q*(1-e))
    if (any(!R)) score[!R] <- -(1-lambda)*residual[!R]/(1-q)
    if (model=="unrestricted") {
      correction <- -(R-q)/(1-q)
      correction[I] <- correction[I]+(residual[I]^2-1)/(1-e)
      correction[!R] <- correction[!R]-q/(1-q)*(residual[!R]^2/external_variance-1)
      score <- score+delta*lambda*(1-lambda)*correction/q
    }
    # For the ECT this directly compares the observed treatment and external
    # regressions 1 + 0.4x and 0.4x - delta. Its causal gap is not identified.
    target <- if (model=="ect") 1+delta else 1+(1-lambda)*delta
    stopifnot(abs(sum(weight)-1)<1e-12,abs(sum(weight*score))<1e-12)
    data.frame(model=model,internal_share=lambda,target=target,
               gap=if (model=="ect") NA_real_ else target-1,
               bound_per_record=sum(weight*score^2))
  })
}
results <- do.call(rbind,lapply(c("unrestricted","transportable","rct","ect"),model_example))
# From 300 iid observations, the expected retained counts are 150 (RCT) and 250 (ECT).
# Scale each design's iid bound using these planning counts.
results$N <- c(300,300,150,250)
results$variance <- results$bound_per_record/results$N
results$standard_error <- sqrt(results$variance)
print(results,row.names=FALSE)
stopifnot(isTRUE(all.equal(results$bound_per_record,c(4.528125,4.5,4.5,25/6))))
