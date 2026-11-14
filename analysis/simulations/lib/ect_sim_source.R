# ECT AIPW functions retained from the original implementation.
# This study calls the quadratic OLS outcome model and logistic source model.
# Other branches remain unchanged for provenance; they are not public entry points.

truncate_prob <- function(p, eps = 0.01) {
  p <- as.numeric(p)
  pmin(1 - eps, pmax(eps, p))
}

as_numeric_matrix <- function(df, cols) {
  X <- data.matrix(df[, cols, drop = FALSE])
  storage.mode(X) <- "double"
  X
}

make_quad_features <- function(X) {
  out <- cbind(X, X^2)
  storage.mode(out) <- "double"
  out
}

fit_or_model <- function(D0,
                         or_method = c("ridge_glmnet", "lin_lm", "quad_lm", "nnet", "krr_rbf", "bart"),
                         lambda_or = 1e-3,
                         krr_lambda = 1e-3,
                         krr_ell = NULL,
                         nnet_size = 5,
                         nnet_decay = 1e-3,
                         nnet_maxit = 300,
                         ## --- BART ---
                         bart_ntree = 200,
                         bart_k = 2,
                         bart_nskip = 100,
                         bart_ndpost = 200) {
  or_method <- match.arg(or_method)
  Xnames <- grep("^X\\d+$", names(D0), value = TRUE)
  X0 <- as_numeric_matrix(D0, Xnames)
  y0 <- as.numeric(D0$Y)

  if (or_method == "ridge_glmnet") {
    X0f <- make_quad_features(X0)
    fit <- glmnet::glmnet(
      x = X0f, y = y0, alpha = 0, lambda = lambda_or,
      standardize = TRUE, intercept = TRUE
    )
    return(list(method = or_method, fit = fit, Xnames = Xnames, lambda = lambda_or))
  }

  if (or_method == "lin_lm") {
    fit <- stats::lm(y ~ ., data = data.frame(y = y0, X0))
    return(list(method = or_method, fit = fit, Xnames = Xnames))
  }

  if (or_method == "quad_lm") {
    X0f <- make_quad_features(X0)
    fit <- stats::lm(y ~ ., data = data.frame(y = y0, X0f))
    return(list(method = or_method, fit = fit, Xnames = Xnames))
  }

  if (or_method == "nnet") {
    if (!requireNamespace("nnet", quietly = TRUE)) stop("Install 'nnet' for or_method='nnet'.")
    X0f <- make_quad_features(X0)
    mu <- colMeans(X0f)
    sdv <- apply(X0f, 2, stats::sd)
    sdv[sdv == 0] <- 1
    Xs <- scale(X0f, center = mu, scale = sdv)
    fit <- nnet::nnet(
      x = Xs, y = y0, size = nnet_size, decay = nnet_decay,
      linout = TRUE, maxit = nnet_maxit, trace = FALSE
    )
    return(list(method = or_method, fit = fit, Xnames = Xnames, mu = mu, sdv = sdv))
  }

  if (or_method == "krr_rbf") {
    Xtr <- make_quad_features(X0)
    storage.mode(Xtr) <- "double"
    n <- nrow(Xtr)

    if (is.null(krr_ell)) {
      m <- min(n, 300L)
      Z <- Xtr[sample.int(n, m), , drop = FALSE]
      G <- rowSums(Z^2)
      D2 <- outer(G, G, "+") - 2 * tcrossprod(Z)
      D2[D2 < 0] <- 0
      d <- sqrt(D2[upper.tri(D2)])
      ell_use <- stats::median(d[d > 0])
      if (!is.finite(ell_use) || ell_use <= 0) ell_use <- 1
    } else {
      ell_use <- as.numeric(krr_ell)
      if (!is.finite(ell_use) || ell_use <= 0) stop("krr_ell must be positive.")
    }

    lambda_use <- as.numeric(krr_lambda)
    if (!is.finite(lambda_use) || lambda_use <= 0) stop("krr_lambda must be positive.")

    mhat <- fit_krr_intercept_predictor(Xtr, y0, ell = ell_use, lambda = lambda_use)

    return(list(
      method  = or_method,
      Xnames  = Xnames,
      ell     = ell_use,
      lambda  = lambda_use,
      predict = mhat
    ))
  }

  if (or_method == "bart") {
    fit <- dbarts::bart(
      x.train = X0,
      y.train = y0,
      ntree = as.integer(bart_ntree),
      k = as.numeric(bart_k),
      nskip = as.integer(bart_nskip),
      ndpost = as.integer(bart_ndpost),
      keeptrees = TRUE,
      verbose = FALSE
    )
    return(list(
      method = or_method,
      fit = fit,
      Xnames = Xnames,
      bart_ntree = as.integer(bart_ntree),
      bart_k = as.numeric(bart_k),
      bart_nskip = as.integer(bart_nskip),
      bart_ndpost = as.integer(bart_ndpost)
    ))
  }

  stop("Unknown OR method.")
}

predict_or_model <- function(or_obj, Dnew) {
  X <- as_numeric_matrix(Dnew, or_obj$Xnames)

  if (or_obj$method == "ridge_glmnet") {
    Xf <- make_quad_features(X)
    return(as.numeric(glmnet::predict.glmnet(or_obj$fit, newx = Xf, s = or_obj$lambda)))
  }
  if (or_obj$method == "lin_lm") {
    return(as.numeric(stats::predict(or_obj$fit, newdata = data.frame(X))))
  }
  if (or_obj$method == "quad_lm") {
    Xf <- make_quad_features(X)
    return(as.numeric(stats::predict(or_obj$fit, newdata = data.frame(Xf))))
  }
  if (or_obj$method == "nnet") {
    Xf <- make_quad_features(X)
    Xs <- scale(Xf, center = or_obj$mu, scale = or_obj$sdv)
    return(as.numeric(stats::predict(or_obj$fit, Xs)))
  }
  if (or_obj$method == "krr_rbf") {
    Xte <- make_quad_features(X)
    storage.mode(Xte) <- "double"
    return(as.numeric(or_obj$predict(Xte)))
  }
  if (or_obj$method == "bart") {
    ## dbarts predict returns (ndpost x n) matrix; take posterior mean
    yh <- predict(or_obj$fit, newdata = X)
    return(as.numeric(colMeans(yh)))
  }
  stop("Unknown OR object/method.")
}

fit_ps_model <- function(D1, D0,
                         ps_method = c("glm", "quad_glm", "ridge_glmnet"),
                         lambda_ps = 1e-3) {
  ps_method <- match.arg(ps_method)

  Xnames <- grep("^X\\d+$", names(D1), value = TRUE)
  X1 <- as_numeric_matrix(D1, Xnames)
  X0 <- as_numeric_matrix(D0, Xnames)
  X <- rbind(X1, X0)
  A <- c(rep(1, nrow(D1)), rep(0, nrow(D0)))

  if (ps_method == "glm") {
    fit <- stats::glm(A ~ ., data = data.frame(A = A, X), family = stats::binomial())
    return(list(method = ps_method, fit = fit, Xnames = Xnames))
  }
  if (ps_method == "quad_glm") {
    Xf <- make_quad_features(X)
    fit <- stats::glm(A ~ ., data = data.frame(A = A, Xf), family = stats::binomial())
    return(list(method = ps_method, fit = fit, Xnames = Xnames))
  }
  if (ps_method == "ridge_glmnet") {
    fit <- glmnet::glmnet(
      x = X, y = A, family = "binomial",
      alpha = 0, lambda = lambda_ps,
      standardize = TRUE, intercept = TRUE
    )
    return(list(method = ps_method, fit = fit, Xnames = Xnames, lambda = lambda_ps))
  }

  stop("Unknown PS method.")
}

predict_ps_model <- function(ps_obj, Dnew) {
  X <- as_numeric_matrix(Dnew, ps_obj$Xnames)

  if (ps_obj$method == "glm") {
    return(as.numeric(stats::predict(ps_obj$fit, newdata = data.frame(X), type = "response")))
  }
  if (ps_obj$method == "quad_glm") {
    Xf <- make_quad_features(X)
    return(as.numeric(stats::predict(ps_obj$fit, newdata = data.frame(Xf), type = "response")))
  }
  if (ps_obj$method == "ridge_glmnet") {
    p <- glmnet::predict.glmnet(ps_obj$fit, newx = X, s = ps_obj$lambda, type = "response")
    return(as.numeric(p))
  }

  stop("Unknown PS object/method.")
}

predict_qhat <- function(ps_obj, Dnew, eps = 0.01) {
  pi_hat <- truncate_prob(predict_ps_model(ps_obj, Dnew), eps = eps)
  pi_hat / (1 - pi_hat)
}

att_hat_dr <- function(D1, D0, or_obj, ps_obj, eps = 0.01) {
  f1 <- predict_or_model(or_obj, D1)
  f0 <- predict_or_model(or_obj, D0)
  q0 <- predict_qhat(ps_obj, D0, eps = eps)

  mean(D1$Y - f1) - sum(q0 * (D0$Y - f0)) / nrow(D1)
  # mean(D1$Y - f1) - (nrow(D0) / nrow(D1)) * mean(q0 * (D0$Y - f0)) # Equivalent form
}

eif_wald_ci_dr <- function(D1, D0, or_obj, ps_obj, theta_hat, alpha = 0.05, eps = 0.01) {
  N1 <- nrow(D1)
  N0 <- nrow(D0)
  N <- N1 + N0
  p1_hat <- N1 / N

  f1 <- predict_or_model(or_obj, D1)
  f0 <- predict_or_model(or_obj, D0)
  q0 <- predict_qhat(ps_obj, D0, eps = eps)

  phi1 <- (1 / p1_hat) * (D1$Y - f1 - theta_hat)
  phi0 <- -(1 / p1_hat) * q0 * (D0$Y - f0)
  phi <- c(phi1, phi0)

  Vhat <- mean(phi^2)
  se <- sqrt(Vhat / N)
  z <- stats::qnorm(1 - alpha / 2)
  ci <- c(lower = theta_hat - z * se, upper = theta_hat + z * se)

  list(ci = ci, se = se)
}
