suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(matrixcalc))

conf <- list(inf = 10^308, # Approximating infinity
             epsilon = 1e-10, # Small value to compare parameters (tolerance)
             zero = 1e-8, # Lower bound for sigmag2 and sigmad2 in optim
             lambda = 1e-6, # Shrinkage parameter for the covariance matrix
             parscale_sigmaC2 = 1e-5,
             parscale_sigmac2 = 1e-5,
             parscale_sigmab2 = 1e-5,
             maxit = 5e3,
             ndeps = 1e-8)

##### ADD THE RUNNING TIME OF THE FUNCTION #####
FUN_with_time <- function(FUN, ...) {
    start <- Sys.time()
    par <- FUN(...)
    end <- Sys.time()
    
    par$time <- as.numeric(difftime(end, start, units = "secs"))
    
    return(par)
}

##### APPROXIMATING THE MEDIATOR CORRELATION MATRIX #####
estimate_Sigma <- function(med_hat) {
    k <- ncol(med_hat)
    if (k <= 0)
        stop("Estimating Sigma requires a positive number of mediators.")
    
    Sigma <- diag(k)
    
    if (nrow(med_hat) <= 1)
        warning("The number of observations are <= 1, thus Sigma is taken to be identity.")
    else if (k > 1) {
        Sigma <- cor(med_hat, use = "pairwise.complete.obs")[, , drop = FALSE]
        
        if (rankMatrix(Sigma)[1] != k) {
            warning(paste0("Estimated Sigma is not full rank. Applying shrinkage (lambda = ", conf$lambda, ") to make Sigma invertible."))
            
            Sigma <- conf$lambda * diag(k) + (1-conf$lambda)*Sigma
        }
    }
    
    return(Sigma)
}

##### INTEGRATED LIKELIHOOD APPROACH TO MEDIATION ANALYSIS #####
get_L11 <- function(B_hat, mu_sb2, xi1, k, l) {
    if (l > 1.5*k)
        L11 <- 1/mu_sb2 * (diag(l) - 1/mu_sb2 * B_hat %*% solve(diag(1/xi1, k) + 1/mu_sb2 * t(B_hat) %*% B_hat) %*% t(B_hat))
    else
        L11 <- solve(diag(mu_sb2, l) + xi1 * B_hat %*% t(B_hat))
    
    return(L11)
}

ll_integrated_C <- function(x, beta_hat, C_hat, k, m, nX, nZ, sigmaC2, pleiotropy) {
    sigmag2 <- x[1]
    if (pleiotropy)
        sigmaC2 <- x[2]
    
    sZ2 <- 1/nZ + sigmaC2
    sC2 <- sZ2 + 1/nX*sigmag2
    
    xib <- sigmag2 / (sC2 + sigmag2 * sum(beta_hat^2))
    xi3 <- 1 + sigmag2/sC2 * sum(beta_hat^2)
    
    c(-1/2*(m*k*log(2*pi) + m*k*log(sC2) + k*log(xi3)) - 
          1/(2*sC2) * t(c(C_hat)) %*% c(C_hat - xib * beta_hat %*% t(beta_hat) %*% C_hat))
}

ll_integrated_b <- function(x, B_hat, b_hat, k, l, nZ, nY, sigmab2, pleiotropy) {
    sigmad2 <- x[1]
    if (pleiotropy)
        sigmab2 <- x[2]
    
    mu_sb2 <- 1/nY + sigmab2 + 1/nZ * sigmad2 * k # For identity Sigma, k = trace(Sigma) = sum(diag(Sigma))
    
    S11 <- diag(mu_sb2, l) + sigmad2 * B_hat %*% t(B_hat)
    mu <- rep(0, l)
    
    c(-1/2 * (l * log(2*pi) + determinant(S11, logarithm = TRUE)$modulus[1]) - 
          1/2 * (t(b_hat - mu) %*% solve(S11) %*% (b_hat - mu)))
}

ll_integrated_c <- function(alpha, sigmagd, sigmag2, sigmad2, beta_hat, c_hat, k, m, nX, nY, sigmac2) {
    sigmaw2 <- k * (sigmag2*sigmad2 + sigmagd^2)
    mu_sc2 <- 1/nY + sigmac2 + 1/nX * (sigmaw2 + (alpha + k*sigmagd)^2)
    xic <- sigmaw2 / (mu_sc2 + sigmaw2 * sum(beta_hat^2))
    xiD <- alpha + k*sigmagd
    
    ll <- -1/2 * (m*log(2*pi) + (m-1)*log(mu_sc2) + log(sigmaw2) - log(xic)) +
        xic/(2*mu_sc2) * (sum(c_hat*beta_hat) - xiD*sum(beta_hat^2))^2 -
        1/(2*mu_sc2) * sum((c_hat - xiD*beta_hat)^2)
    
    sum(ll)
}

ll_integrated_bC <- function(sigmagd, sigmag2, sigmad2, beta_hat, C_hat, B_hat, b_hat, k, l, m, nX, nZ, nY, sigmaC2, sigmab2) {
    if (sigmagd^2 > sigmag2*sigmad2)
        return(-sqrt(conf$inf)) # Returning Inf here could result in non-finite finite-difference value error in optim (in Hessian?)
    
    mu_sb2 <- 1/nY + sigmab2 + 1/nZ * sigmad2 * k # For identity Sigma, k = trace(Sigma) = sum(diag(Sigma))
    sZ2 <- 1/nZ + sigmaC2
    sC2 <- sZ2 + 1/nX*sigmag2
    
    xib <- sigmag2 / (sC2 + sigmag2 * sum(beta_hat^2))
    xi1 <- sigmad2 - sigmagd^2/sC2 * sum(beta_hat^2) * (1 - xib * sum(beta_hat^2))
    xi2 <- -sigmagd/sC2 * (1 - xib * sum(beta_hat^2))
    xi3 <- 1 + sigmag2/sC2 * sum(beta_hat^2)
    
    L11 <- get_L11(B_hat, mu_sb2, xi1, k, l)
    lrL11 <- b_hat + xi2 * B_hat %*% t(C_hat) %*% beta_hat
    
    c(-1/2 * (t(lrL11)%*%L11%*%lrL11 + 1/sC2 * t(c(C_hat)) %*% c(C_hat - xib*beta_hat%*%t(beta_hat)%*%C_hat)) -
          1/2 * ((m*k + l)*log(2*pi) + l*log(mu_sb2) + m*k*log(sC2) + k*log(xi3) + determinant(diag(k) + xi1/mu_sb2 * t(B_hat) %*% B_hat, logarithm = TRUE)$modulus[1]))
}

# The integrated method assumes identity Sigma
ll_integrated <- function(x, ll_type, pleiotropy, dots) {
    beta_hat <- dots$beta_hat
    C_hat <- dots$C_hat
    B_hat <- dots$B_hat
    c_hat <- dots$c_hat
    b_hat <- dots$b_hat
    k <- dots$k
    l <- dots$l
    m <- dots$m
    nX <- dots$nX
    nZ <- dots$nZ
    nY <- dots$nY
    
    alpha <- x[1]    
    if (ll_type == "fixed") {
        sigmagd <- x[2]
        sigmag2 <- dots$sigmag2
        sigmad2 <- dots$sigmad2
    }
    else if (ll_type == "just_alpha") {
        sigmagd <- (dots$total - alpha) / k
        sigmag2 <- dots$sigmag2
        sigmad2 <- dots$sigmad2
    }
    else if (ll_type == "free") {
        rhogd <- x[2]
        sigmag2 <- x[3]
        sigmad2 <- x[4]
        sigmagd <- rhogd * sqrt(sigmag2*sigmad2)
    }
    
    if (pleiotropy) {
        sigmaC2 <- x[length(x) - 2]
        sigmac2 <- x[length(x) - 1]
        sigmab2 <- x[length(x)]
    }
    else {
        sigmaC2 <- dots$sigmaC2
        sigmac2 <- dots$sigmac2
        sigmab2 <- dots$sigmab2
    }
    
    ll_value <- ll_integrated_c(alpha, sigmagd, sigmag2, sigmad2, beta_hat, c_hat, k, m, nX, nY, sigmac2) + 
        ll_integrated_bC(sigmagd, sigmag2, sigmad2, beta_hat, C_hat, B_hat, b_hat, k, l, m, nX, nZ, nY, sigmaC2, sigmab2)
    
    # Replace +-Inf return values with approximations
    sign(ll_value) * min(abs(ll_value), conf$inf)
}

# LRT is used to calculate SEs (also tried Hessian but it doesn't work as well)
ll_integrated_SE_LRT <- function(par, ll_type, pleiotropy, dots) {
    ll_complex <- ll_integrated(par, ll_type, pleiotropy, dots)
    
    ll_null <- ll_integrated(x = c(0, par[-1]), ll_type, pleiotropy, dots) # Replace estimated direct effect (alpha) with 0
    if (ll_null >= ll_complex) {
        warning("Estimating SE of the direct effect with LRT but null model has higher likelihood than full model. Returning NaN(s). Consider optimizing the full model with 0 seed.")
        return(list(total = ifelse(ll_type == "just_alpha", NA_real_, NaN), alpha = NaN))
    }
    Z <- sqrt(-2*(ll_null - ll_complex))
    SE_alpha <- abs(par[1] / Z)
    
    if (ll_type == "just_alpha")
        return(list(total = NA_real_, alpha = SE_alpha))
    
    ll_null <- ll_integrated(x = c(par[1], 0, par[-(1:2)]), ll_type, pleiotropy, dots) # Replace estimated indirect effect (sigmagd or rhogd) with 0
    if (ll_null >= ll_complex) {
        warning("Estimating SE of the total indirect effect with LRT but null model has higher likelihood than full model. Returning NaN. Consider optimizing the full model with 0 seed.")
        return(list(total = NaN, alpha = SE_alpha))
    }
    Z <- sqrt(-2*(ll_null - ll_complex))
    indirect <- dots$k * ifelse(ll_type != "free", par[2], par[2] * sqrt(par[3]*par[4]))
    SE_indirect <- abs(indirect / Z)
    
    SE_total <- sqrt(SE_indirect^2 + SE_alpha^2)
    
    list(total = SE_total, alpha = SE_alpha)
}

sigmag2_approach_integrated <- function(beta_hat, C_hat, k, m, nX, nZ, sigmaC2, pleiotropy) {
    par <- runif(1) # sigmag2 seed is random
    lower <- conf$zero
    upper <- 1
    parscale <- 1
    ndeps <- conf$ndeps
    sigmag2_method <- "Brent" # Use Brent when optimizing over one variable: it's more stable, not affected by parscale and produces higher optim values
    
    if (pleiotropy) {
        par <- c(par, sigmaC2)
        lower <- c(lower, 0)
        upper <- c(upper, 1)
        parscale <- c(parscale, conf$parscale_sigmaC2)
        ndeps <- c(ndeps, conf$ndeps)
        sigmag2_method <- "L-BFGS-B"
    }
    
    control <- list(parscale=parscale, ndeps=ndeps, fnscale=-1)
    sigmag2_optim <- optim(par=par, fn=ll_integrated_C, method=sigmag2_method, lower=lower, upper=upper, control=control, 
                           beta_hat=beta_hat, C_hat=C_hat, k=k, m=m, nX=nX, nZ=nZ, sigmaC2=sigmaC2, pleiotropy=pleiotropy)
    
    list(sigmag2=sigmag2_optim$par[1], convergence=sigmag2_optim$convergence)
}

sigmad2_approach_integrated <- function(B_hat, b_hat, k, l, nZ, nY, sigmab2, pleiotropy) {
    par <- runif(1) # sigmad2 seed is random
    lower <- conf$zero
    upper <- 1
    parscale <- 1
    ndeps <- conf$ndeps
    sigmad2_method <- "Brent" # Use Brent when optimizing over one variable: it's more stable, not affected by parscale and produces higher optim values
    
    if (pleiotropy) {
        par <- c(par, sigmab2)
        lower <- c(lower, 0)
        upper <- c(upper, 1)
        parscale <- c(parscale, conf$parscale_sigmab2)
        ndeps <- c(ndeps, conf$ndeps)
        sigmad2_method <- "L-BFGS-B"
    }
    
    control <- list(parscale=parscale, ndeps=ndeps, fnscale=-1)
    sigmad2_optim <- optim(par=par, fn=ll_integrated_b, method=sigmad2_method, lower=lower, upper=upper, control=control,
                           B_hat=B_hat, b_hat=b_hat, k=k, l=l, nZ=nZ, nY=nY, sigmab2=sigmab2, pleiotropy=pleiotropy)
    
    list(sigmad2=sigmad2_optim$par[1], convergence=sigmad2_optim$convergence)
}

ll_approach_integrated <- function(beta_hat, C_hat, B_hat, c_hat, b_hat, nX, nZ, nY, total, alpha, sigmag2 = NULL, 
                                   sigmad2 = NULL, ll_type = c("fixed", "free", "just_alpha"), pleiotropy = FALSE) {
    ll_type <- match.arg(ll_type)
    sigmaC2 <- 0
    sigmac2 <- 0
    sigmab2 <- 0
    k <- ncol(B_hat)
    l <- nrow(B_hat)
    m <- nrow(C_hat)
    sigmag2_convergence <- NA_real_
    sigmad2_convergence <- NA_real_
    
    if (is.null(sigmag2)) {
        sigmag2_integrated <- sigmag2_approach_integrated(beta_hat, C_hat, k, m, nX, nZ, sigmaC2, pleiotropy)
        sigmag2 <- sigmag2_integrated$sigmag2
        sigmag2_convergence <- ifelse(ll_type != "free", sigmag2_integrated$convergence, sigmag2_convergence) # Don't need convergence for the free method
    }
    if (is.null(sigmad2)) {
        sigmad2_integrated <- sigmad2_approach_integrated(B_hat, b_hat, k, l, nZ, nY, sigmab2, pleiotropy)
        sigmad2 <- sigmad2_integrated$sigmad2
        sigmad2_convergence <- ifelse(ll_type != "free", sigmad2_integrated$convergence, sigmad2_convergence) # Don't need convergence for the free method
    }
    
    dots <- list(beta_hat=beta_hat, C_hat=C_hat, B_hat=B_hat, c_hat=c_hat, b_hat=b_hat, k=k, l=l, m=m, nX=nX, nZ=nZ, nY=nY, 
                 sigmaC2=sigmaC2, sigmac2=sigmac2, sigmab2=sigmab2, sigmag2=sigmag2, sigmad2=sigmad2)
    sigmagd <- (total - alpha) / k
    is_seed_OK <- sigmagd^2 <= sigmag2 * sigmad2 # Note that we actually force the seed to be OK below but it's useful to know what was the initial state
    
    if (ll_type == "fixed") {
        max_sigmagd <- sqrt(sigmag2 * sigmad2) # Corresponds to correlation 1
        sigmagd <- ifelse(is_seed_OK, sigmagd, sign(sigmagd) * (max_sigmagd - conf$epsilon))
        
        par <- c(alpha, sigmagd)
        lower <- c(-1, -max_sigmagd)
        upper <- c(1, max_sigmagd)
        parscale = c(1, abs(sigmagd/alpha)) # Could use c(1, 1/k) or c(1, 1) here
    }
    else if (ll_type == "just_alpha") {
        max_indirect <- k*sqrt(sigmag2*sigmad2)
        min_indirect <- -max_indirect
        dots$total <- total
        
        par <- alpha
        lower <- total - max_indirect
        upper <- total - min_indirect
        parscale <- c(1)
    }
    else if (ll_type == "free") {
        rhogd <- ifelse(is_seed_OK, sigmagd / sqrt(sigmag2 * sigmad2), sign(sigmagd) * (1 - conf$epsilon))
        par <- c(alpha, rhogd, sigmag2, sigmad2)
        lower <- c(-1, -1, conf$zero, conf$zero)
        upper <- c(1, 1, 1, 1)
        parscale <- abs(c(1, rhogd/alpha, sigmag2/alpha, sigmad2/alpha))
    }
    
    if (pleiotropy) {
        par <- c(par, sigmaC2, sigmac2, sigmab2)
        lower <- c(lower, 0, 0, 0)
        upper <- c(upper, 1, 1, 1)
        parscale <- c(parscale, conf$parscale_sigmaC2, conf$parscale_sigmac2, conf$parscale_sigmab2)
    }
    
    control <- list(maxit=conf$maxit, parscale=parscale, ndeps=rep(conf$ndeps, length(par)), fnscale=-1)
    optimization <- optim(par=par, fn=ll_integrated, ll_type=ll_type, pleiotropy=pleiotropy, dots=dots, method="L-BFGS-B", lower=lower, upper=upper, control=control)
    SEs <- ll_integrated_SE_LRT(par=optimization$par, ll_type=ll_type, pleiotropy=pleiotropy, dots=dots)
    
    alpha <- optimization$par[1]
    if (ll_type == "fixed")
        sigmagd <- optimization$par[2]
    else if (ll_type == "just_alpha")
        sigmagd <- (total - alpha) / k
    else if (ll_type == "free") {
        sigmag2 <- optimization$par[3]
        sigmad2 <- optimization$par[4]
        sigmagd <- optimization$par[2] * sqrt(sigmag2 * sigmad2)
    }
    total = alpha + k*sigmagd
    sigmag2_bound <- abs(sigmag2 - conf$zero) <= conf$epsilon || abs(sigmag2 - 1) <= conf$epsilon
    sigmad2_bound <- abs(sigmad2 - conf$zero) <= conf$epsilon || abs(sigmad2 - 1) <= conf$epsilon
    
    if (pleiotropy) {
        sigmaC2 <- optimization$par[length(optimization$par) - 2]
        sigmac2 <- optimization$par[length(optimization$par) - 1]
        sigmab2 <- optimization$par[length(optimization$par)]
    }
    
    list(convergence=optimization$convergence, optim_value=optimization$value, is_seed_OK=is_seed_OK,
         sigmag2_convergence=sigmag2_convergence, sigmad2_convergence=sigmad2_convergence, sigmag2_bound=sigmag2_bound, sigmad2_bound=sigmad2_bound,
         total=total, alpha=alpha, SE_total=SEs$total, SE_alpha=SEs$alpha, sigmag2=sigmag2, sigmad2=sigmad2, sigmaC2=sigmaC2, sigmac2=sigmac2, sigmab2=sigmab2)
}

##### ORIGINAL LIKELIHOOD APPROACH TO MEDIATION ANALYSIS #####
ll_original_b <- function(delta, B_hat, b_hat, nZ, nY, sigmab2, Sigma) {
    sb2 <- 1/nY + sigmab2 + 1/nZ * c(t(delta) %*% Sigma %*% delta)
    
    ll <- -1/2*log(2*pi) - 1/2*log(sb2) - 1/(2*sb2) * (b_hat - B_hat %*% delta)^2
    
    sum(ll)
}

ll_original_c <- function(gamma, delta, alpha, beta_hat, c_hat, nX, nY, sigmac2) {
    sc2 <- c(1/nY + sigmac2 + 1/nX * c(alpha + sum(gamma*delta))^2)
    
    ll <- -1/2*log(2*pi) - 1/2*log(sc2) - 1/(2*sc2) * c(c_hat - beta_hat * (alpha + sum(gamma*delta)))^2
    
    sum(ll)
}

ll_original_C <- function(gamma, beta_hat, C_hat, k, nX, nZ, sigmaC2, invSigma, logdetSigma) {
    Chat_invSigma <- C_hat %*% invSigma
    Chat_invSigma_Chat <- rowSums(Chat_invSigma * C_hat)
    sC <- nX * (1/nZ + sigmaC2) + c(t(gamma) %*% invSigma %*% gamma)
    
    ll <- -k/2*log(2*pi) - 1/2*logdetSigma - (k-1)/2*log(1/nZ + sigmaC2) -
        1/2*log(1/nZ + sigmaC2 + 1/nX*c(t(gamma) %*% invSigma %*% gamma)) -
        1/2/sC*nX * (Chat_invSigma_Chat - 2*beta_hat * (Chat_invSigma %*% gamma) + beta_hat^2 * c(t(gamma) %*% invSigma %*% gamma)) -
        1/2/sC/(1/nZ + sigmaC2) * (Chat_invSigma_Chat * c(t(gamma) %*% invSigma %*% gamma) - (Chat_invSigma %*% gamma)^2)
    
    sum(ll)
}

ll_original <- function(x, pleiotropy, dots) {
    beta_hat <- dots$beta_hat
    C_hat <- dots$C_hat
    B_hat <- dots$B_hat
    c_hat <- dots$c_hat
    b_hat <- dots$b_hat
    k <- dots$k
    nX <- dots$nX
    nZ <- dots$nZ
    nY <- dots$nY
    Sigma <- dots$Sigma
    invSigma <- dots$invSigma
    logdetSigma <- dots$logdetSigma
    
    gamma <- x[0:k]
    delta <- x[ifelse(k == 0, 0, k + 1):(2*k)]
    alpha <- x[2*k + 1]
    
    if (pleiotropy) {
        sigmaC2 <- x[2*k + 2]
        sigmac2 <- x[2*k + 3]
        sigmab2 <- x[2*k + 4]
    }
    else {
        sigmaC2 <- dots$sigmaC2
        sigmac2 <- dots$sigmac2
        sigmab2 <- dots$sigmab2
    }
    
    ll_value <- ll_original_b(delta, B_hat, b_hat, nZ, nY, sigmab2, Sigma) + 
        ll_original_c(gamma, delta, alpha, beta_hat, c_hat, nX, nY, sigmac2) + 
        ll_original_C(gamma, beta_hat, C_hat, k, nX, nZ, sigmaC2, invSigma, logdetSigma)
    
    # Replace +-Inf return values with approximations
    sign(ll_value) * min(abs(ll_value), conf$inf)
}

# Optimization must have been on log scale and maximizing
ll_original_SE_hessian <- function(hessian, k) {
    if (is.null(hessian))
        return(list(total = NA_real_, alpha = NA_real_))
    if (nrow(hessian) < 2*k + 1)
        stop("The number of mediators is not consistent with the size of Hessian.")
    if (nrow(hessian) > 2*k + 1)
        warning("Using the first 2k + 1 rows and columns of the Hessian.")
    if (k == 0)
        return(list(total = NA_real_, alpha = sqrt(-1/hessian[1, 1])))
    
    neg_hessian <- -hessian
    if (rankMatrix(neg_hessian)[1] != nrow(neg_hessian))
        neg_hessian <- conf$lambda*diag(nrow(neg_hessian)) + (1-conf$lambda)*neg_hessian
    
    covmat <- solve(neg_hessian)
    if (any(diag(covmat) < 0)) {
        warning("The covariance matrix approximation has negative values on the diagonal. Returning NaNs.")
        return(list(total = NaN, alpha = NaN))
    }
    
    Var_omega <- sum(diag(covmat)[1:k] * diag(covmat)[(k+1):(2*k)] + covmat[cbind(1:k, (k+1):(2*k))]^2)
    SE_alpha <- sqrt(covmat[2*k + 1, 2*k + 1])
    SE_total <- sqrt(SE_alpha**2 + Var_omega)
    
    list(total = SE_total, alpha = SE_alpha)
}

ll_original_SE_LRT <- function(gamma, delta, alpha, sigmac2, dots) {
    ll_complex <- ll_original_c(gamma, delta, alpha, dots$beta_hat, dots$c_hat, dots$nX, dots$nY, sigmac2)
    
    # Replace estimated direct effect (alpha) with 0
    ll_null <- ll_original_c(gamma, delta, 0, dots$beta_hat, dots$c_hat, dots$nX, dots$nY, sigmac2)
    if (ll_null >= ll_complex) {
        warning("Estimating SE of the direct effect with LRT but null model has higher likelihood than full model. Returning NaNs.")
        return(list(total = NaN, alpha = NaN))
    }
    Z <- sqrt(-2*(ll_null - ll_complex))
    SE_alpha <- abs(alpha / Z)
    
    # Replace estimated total indirect effect with 0: ll_original_c only uses the total indirect effect = sum(gamma*delta)
    ll_null <- ll_original_c(rep(0, length(gamma)), rep(0, length(delta)), alpha, dots$beta_hat, dots$c_hat, dots$nX, dots$nY, sigmac2)
    if (ll_null >= ll_complex) {
        warning("Estimating SE of the total indirect effect with LRT but null model has higher likelihood than full model. Returning NaN.")
        return(list(total = NaN, alpha = SE_alpha))
    }
    Z <- sqrt(-2*(ll_null - ll_complex))
    SE_indirect <- abs(sum(gamma * delta) / Z)
    
    SE_total <- sqrt(SE_indirect^2 + SE_alpha^2)
    
    list(total = SE_total, alpha = SE_alpha)
}

ll_approach_original <- function(beta_hat, C_hat, B_hat, c_hat, b_hat, nX, nZ, nY, Sigma, gamma, delta, alpha, pleiotropy = FALSE) {
    k = ncol(B_hat)
    sigmaC2 <- 0
    sigmac2 <- 0
    sigmab2 <- 0
    invSigma <- solve(Sigma)
    logdetSigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]
    
    par <- c(gamma, delta, alpha) # gamma, delta, alpha
    lower <- rep(-1, length(par))
    upper <- rep(1, length(par))
    parscale <- rep(1, length(par))
    ndeps <- rep(conf$ndeps, length(par))
    control <- list(maxit=conf$maxit, ndeps=ndeps, parscale=parscale, fnscale=-1)
    dots <- list(beta_hat=beta_hat, C_hat=C_hat, B_hat=B_hat, c_hat=c_hat, b_hat=b_hat, k=k, nX=nX, nZ=nZ, nY=nY, Sigma=Sigma, invSigma=invSigma, logdetSigma=logdetSigma)
    
    if (pleiotropy) {
        par <- c(par, sigmaC2, sigmac2, sigmab2)
        lower <- c(lower, 0, 0, 0)
        upper <- c(upper, 1, 1, 1)
        control$parscale <- c(control$parscale, conf$parscale_sigmaC2, conf$parscale_sigmac2, conf$parscale_sigmab2)
        control$ndeps <- c(control$ndeps, conf$ndeps, conf$ndeps, conf$ndeps)
        
        optimization <- optim(par=par, fn=ll_original, pleiotropy=pleiotropy, dots=dots, method="L-BFGS-B", lower=lower, upper=upper, control=control)
        
        sigmaC2 <- optimization$par[2*k+2]
        sigmac2 <- optimization$par[2*k+3]
        sigmab2 <- optimization$par[2*k+4]
    }
    else {
        dots <- c(dots, list(sigmaC2=sigmaC2, sigmac2=sigmac2, sigmab2=sigmab2))
        
        optimization <- optim(par=par, fn=ll_original, pleiotropy=pleiotropy, dots=dots, method="L-BFGS-B", lower=lower, upper=upper, control=control)
    }
    
    gamma = optimization$par[1:k]
    delta = optimization$par[(k+1):(2*k)]
    alpha = optimization$par[2*k+1]
    
    # Calculating Hessian can fail separately which is why we do it outside of optim
    # SEs_hessian <- tryCatch(optimHess(par=optimization$par, fn=ll_original, pleiotropy=pleiotropy, dots=dots, control=control), error = function(e) NULL) %>%
    #     ll_original_SE_hessian(k)
    SEs <- ll_original_SE_LRT(gamma, delta, alpha, sigmac2, dots)
    
    list(total=sum(gamma*delta) + alpha, alpha=alpha, gamma=gamma, delta=delta, SE_total=SEs$total, SE_alpha=SEs$alpha, 
         sigmag2=var(gamma), sigmad2=var(delta), sigmaC2=sigmaC2, sigmac2=sigmac2, sigmab2=sigmab2,
         convergence=optimization$convergence, optim_value=optimization$value)
}

# Wald ratio test or IVW MR when beta_X is a vector, otherwise MVMR
ivw_mr <- function(beta_X, beta_y, var_beta_y) {
    if (length(beta_X) == 1) { # Wald ratio, equivalent to TwoSampleMR::mr_wald_ratio
        beta <- beta_y / beta_X
        SE <- sqrt(var_beta_y / beta_X^2) # Delta method
    }
    else { # IVW regression, equivalent to TwoSampleMR::mr_ivw
        ivw_mr_model <- summary(lm(beta_y ~ 0 + beta_X, weights = 1 / var_beta_y))
        beta <- ivw_mr_model$coefficients[, "Estimate"]
        SE <- ivw_mr_model$coefficients[, "Std. Error"] / min(1, ivw_mr_model$sigma) # The residual SE in IVW regression, 1/N*sum(w_i * e_i^2), has expectation 1: 1/N*sum(1/sigma_i^2 * sigma_i^2) = 1
    }
    
    list(beta = unname(beta), SE = unname(SE), P = 2*pnorm(-abs(beta / SE)))
}

##### NAIVE APPROACH TO MEDIATION ANALYSIS: IVW FOR TOTAL AND MVMR FOR DIRECT EFFECT #####
naive_approach <- function(beta_hat, beta_med_hat, C_hat, B_hat, c_hat, b_hat, var_C_hat, var_c_hat, var_b_hat) {
    # IVW MR for the total effect
    ivw_mr_result <- ivw_mr(beta_hat, c_hat, var_c_hat)
    total <- ivw_mr_result$beta
    gamma <- 1:ncol(C_hat) %>%
        sapply(function(i) ivw_mr(beta_hat, C_hat[, i], var_C_hat[, i])$beta)
    
    # IVW MVMR for the direct effect
    G <- cbind(c(beta_hat, beta_med_hat), rbind(C_hat, B_hat))
    g <- c(c_hat, b_hat)
    var_g <- c(var_c_hat, var_b_hat)
    ivw_mvmr_results <- ivw_mr(G, g, var_g)
    alpha <- ivw_mvmr_results$beta[1]
    delta <- ivw_mvmr_results$beta[-1]
    
    list(total=total, alpha=alpha, gamma=gamma, delta=delta, 
         SE_total=ivw_mr_result$SE, SE_alpha=ivw_mvmr_results$SE[1],
         sigmag2=var(gamma), sigmad2=var(delta))
}

Burgess_MVMR = function(beta_hat, beta_med_hat, C_hat, B_hat, c_hat, b_hat, 
                        var_beta_hat, var_beta_med_hat, var_C_hat, var_B_hat, var_c_hat, var_b_hat,
                        max_iter = 100, no_ini = 1){
    mv_norm = function(n, mu, Sigma){
        d = dim(Sigma)[1]
        if (length(mu) != d){
            stop('mu and Sigma must be the same dimension.')
        }
        A = chol(Sigma)
        Z = sapply(1:d, function(j){rnorm(n, 0, 1)})
        M = matrix(rep(mu, n), nrow = n, byrow = TRUE)
        X = drop(M) + Z %*% A
    }
    
    bxhat = cbind(c(beta_hat, beta_med_hat), rbind(C_hat, B_hat)) #as.matrix(mrob@betaX)
    byhat = c(c_hat, b_hat) #mrob@betaY
    sebx = sqrt(cbind(c(var_beta_hat, var_beta_med_hat), rbind(var_C_hat, var_B_hat))) #as.matrix(mrob@betaXse)
    seby = sqrt(c(var_c_hat, var_b_hat)) #mrob@betaYse
    p = length(byhat)
    K = dim(bxhat)[2]
    S = diag(seby^-2)
    SigX = lapply(1:p, function(j){diag(sebx[j, ]^2, length(sebx[j, ]))})
    
    l = matrix(nrow = max_iter, ncol = no_ini)
    thest = matrix(nrow = K, ncol = no_ini)
    for (k in 1:no_ini){
        bxtilde = t(matrix(sapply(1:p, function(j){mv_norm(1, bxhat[j, ], SigX[[j]])}), ncol = p))
        for (i in 1:100){
            thest[, k] = solve(t(bxtilde) %*% S %*% bxtilde, t(bxtilde) %*% S %*% byhat)
            l[i, k] = -0.5 * sum(sapply(1:p, function(j){
                (byhat[j] - t(bxhat[j, ]) %*% thest[, k])^2 / (seby[j]^2 + t(thest[, k]) %*% SigX[[j]] %*% thest[, k])
            }))
            for (j in 1:p){
                bxtilde[j, ] = t(solve(thest[, k] %*% t(thest[, k]) / seby[j]^2 + solve(SigX[[j]]), byhat[j] * thest[, k] / seby[j]^2 + solve(SigX[[j]], bxhat[j, ])))
            }
            if (i > 1){
                if (abs(l[i] - l[(i-1)]) < 1e-4) {break}
            }
        }
    }
    k0 = which.max(apply(as.matrix(l[is.na(l[, 1]) == F, ]), 2, max))
    th = thest[, k0]
    
    v = sapply(1:p, function(j){seby[j]^2 + t(th) %*% SigX[[j]] %*% th})
    e = sapply(1:p, function(j){byhat[j] - bxhat[j, ] %*% th})
    t = sapply(1:p, function(j){e[j] / sqrt(v[j])})
    
    dt = matrix(sapply(1:p, function(j){(-v[j] * bxhat[j, ] - e[j] * SigX[[j]] %*% th) / v[j]^(3/2)}), ncol = p)
    B = dt %*% t(dt)
    
    dt2 = vector(length = p, mode = "list")
    for (j in 1:p){
        dt2[[j]] = matrix(nrow = K, ncol = K)
        S = SigX[[j]] %*% th
        for (k in 1:K){
            for (l in 1:K){
                dt2[[j]][k, l] = v[j]^(-3/2) * (-2 * S[l] * bxhat[j, k] + S[k] * bxhat[j, l] - e[j] * SigX[[j]][k, l]) +
                    3 * v[j]^(-5/2) *(v[j] * bxhat[j, k] + e[j] * S[k]) * S[l]
            }
        }
    }
    a = Reduce('+', lapply(1:p, function(j){c(t[j]) * dt2[[j]]}))
    A = (dt %*% t(dt) + a)
    
    Var = solve(A, B) %*% t(solve(A))
    SE <- sqrt(diag(Var))
    
    list(beta = th, SE = SE, P = 2*pnorm(-abs(th / SE))) #return(list("thest" = th, "l" = l, "Var" = Var))
}

Burgess_approach <- function(beta_hat, beta_med_hat, C_hat, B_hat, c_hat, b_hat, 
                             var_beta_hat, var_beta_med_hat, var_C_hat, var_B_hat, var_c_hat, var_b_hat) {
    # IVW MR for the total effect
    ivw_mr_result <- ivw_mr(beta_hat, c_hat, var_c_hat)
    total <- ivw_mr_result$beta
    gamma <- 1:ncol(C_hat) %>%
        sapply(function(i) ivw_mr(beta_hat, C_hat[, i], var_C_hat[, i])$beta)
    
    # Burgess MVMR for the direct effect
    Burgess_mvmr_results <- Burgess_MVMR(beta_hat, beta_med_hat, C_hat, B_hat, c_hat, b_hat, 
                                         var_beta_hat, var_beta_med_hat, var_C_hat, var_B_hat, var_c_hat, var_b_hat)
    alpha <- Burgess_mvmr_results$beta[1]
    delta <- Burgess_mvmr_results$beta[-1]
    
    list(total=total, alpha=alpha, gamma=gamma, delta=delta, 
         SE_total=ivw_mr_result$SE, SE_alpha=Burgess_mvmr_results$SE[1],
         sigmag2=var(gamma), sigmad2=var(delta))
}

mediation_analysis <- function(B_hat, C_hat, b_hat, c_hat, beta_hat, beta_med_hat, Sigma, var_B_hat, var_C_hat, var_b_hat, 
                               var_c_hat, var_beta_hat, var_beta_med_hat, nX, nY, nZ, k_selected, l, m, pleiotropy,
                               method = c("naive", "naive_zero", "Burgess", "original", "original_diagonal", "original_B", "integrated_free", "integrated_fixed_contained", "integrated_fixed_naive", "integrated_just_alpha", "integrated_fixed_PCA")) {
    method <- match.arg(method)
    
    k_selected <- unique(ncol(B_hat), ncol(C_hat), nrow(Sigma), ncol(Sigma), length(gamma), length(delta), k_selected)
    l <- unique(nrow(B_hat), length(b_hat), length(beta_med_hat))
    m <- unique(nrow(C_hat), length(c_hat), length(beta_hat))
    if (!all(length(k_selected) == 1, length(m) == 1, length(l) == 1))
        stop("The dimensions of input data are not consistent.")
    
    # Initialize the parameter vector of results
    par <- list(total = NA, alpha = NA, SE_total = NA, SE_alpha = NA)
    
    if (method == "naive") {
        par <- FUN_with_time(naive_approach, beta_hat, beta_med_hat, C_hat, B_hat, c_hat, b_hat, var_C_hat, var_c_hat, var_b_hat)
    }
    else if (method == "naive_zero") {
        beta_med_hat <- rep(0, length(beta_med_hat))
        par <- FUN_with_time(naive_approach, beta_hat, beta_med_hat, C_hat, B_hat, c_hat, b_hat, var_C_hat, var_c_hat, var_b_hat)
    }
    else if (method == "Burgess") {
        par <- FUN_with_time(Burgess_approach, beta_hat, beta_med_hat, C_hat, B_hat, c_hat, b_hat, 
                             var_beta_hat, var_beta_med_hat, var_C_hat, var_B_hat, var_c_hat, var_b_hat)
    }
    else {
        naive_par <- naive_approach(beta_hat, beta_med_hat, C_hat, B_hat, c_hat, b_hat, var_C_hat, var_c_hat, var_b_hat) # Atm, "naive" is used as the seed for the likelihood methods
        
        if (grepl("^original", method)) {
            if (method == "original_diagonal")
                Sigma <- diag(k_selected)
            else if (method == "original_B" || is.null(Sigma))
                Sigma <- estimate_Sigma(B_hat)
            
            par <- FUN_with_time(ll_approach_original, beta_hat, C_hat, B_hat, c_hat, b_hat, nX, nZ, nY, Sigma, naive_par$gamma, naive_par$delta, naive_par$alpha, pleiotropy = pleiotropy)
        }
        else if (method == "integrated_free") { # Free sigmag2 and sigmad2
            par <- FUN_with_time(ll_approach_integrated, beta_hat, C_hat, B_hat, c_hat, b_hat, nX, nZ, nY, naive_par$total, naive_par$alpha, sigmag2 = NULL, sigmad2 = NULL, ll_type = "free", pleiotropy = pleiotropy)
        }
        else if (method == "integrated_fixed_contained") { # Fixed sigmag2 and sigmad2, contained estimation
            par <- FUN_with_time(ll_approach_integrated, beta_hat, C_hat, B_hat, c_hat, b_hat, nX, nZ, nY, naive_par$total, naive_par$alpha, sigmag2 = NULL, sigmad2 = NULL, ll_type = "fixed", pleiotropy = pleiotropy)
        }
        else if (method == "integrated_fixed_naive" && !is.na(naive_par$sigmag2) && !is.na(naive_par$sigmad2)) { # Fixed sigmag2 and sigmad2, naive estimation
            par <- FUN_with_time(ll_approach_integrated, beta_hat, C_hat, B_hat, c_hat, b_hat, nX, nZ, nY, naive_par$total, naive_par$alpha, sigmag2 = naive_par$sigmag2, sigmad2 = naive_par$sigmad2, ll_type = "fixed", pleiotropy = pleiotropy)
        }
        else if (method == "integrated_just_alpha") {
            par <- FUN_with_time(ll_approach_integrated, beta_hat, C_hat, B_hat, c_hat, b_hat, nX, nZ, nY, naive_par$total, naive_par$alpha, sigmag2 = NULL, sigmad2 = NULL, ll_type = "just_alpha", pleiotropy = pleiotropy)
            par$SE_total <- naive_par$SE_total
        }
        else if (method == "integrated_fixed_PCA") { # PCA solution with fixed sigmag2 and sigmad2, contained estimation
            PCA <- prcomp(B_hat, scale = FALSE) # Corresponds to eigen(cor(B_hat))
            B_hat <- B_hat %*% PCA$rotation
            C_hat <- C_hat %*% PCA$rotation
            par <- FUN_with_time(ll_approach_integrated, beta_hat, C_hat, B_hat, c_hat, b_hat, nX, nZ, nY, naive_par$total, naive_par$alpha, sigmag2 = NULL, sigmad2 = NULL, ll_type = "fixed", pleiotropy = pleiotropy)
        }
        
        par$SE_total_naive <- naive_par$SE_total
        par$SE_alpha_naive <- naive_par$SE_alpha
    }
    
    c(list(method=method, pleiotropy=pleiotropy, k_selected=k_selected, l=l, m=m, nX=nX, nY=nY, nZ=nZ), par)
}
