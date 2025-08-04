library(testthat)

generate_test_data <- function(k = 2, l = 3, m = 4, nX = 300000, nZ = 10000, nY = 350000,
                               sigmaC2 = 1e-5, sigmac2 = 1e-4, sigmab2 = 1e-3) {
    Sigma <- matrix(rnorm(l*k), nrow = l, ncol = k) %>%
        estimate_Sigma()
    invSigma <- solve(Sigma)
    logdetSigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]
    
    alpha <- rnorm(1)
    gamma <- rnorm(k)
    delta <- rnorm(k)
    C_hat <- mvtnorm::rmvnorm(m, sigma = (1/nZ + sigmaC2) * Sigma)
    B_hat <- mvtnorm::rmvnorm(l, sigma = 1/nZ * Sigma)
    beta_hat <- rnorm(m, sd = sqrt(1/nX))
    beta_med_hat <- rnorm(l, sd = sqrt(1/nX))
    c_hat <- rnorm(m, sd = sqrt(1/nY + sigmac2))
    b_hat <- rnorm(l, sd = sqrt(1/nY + sigmab2))
    
    sb2 <- 1/nY + sigmab2 + 1/nZ * c(t(delta) %*% Sigma %*% delta)
    sc2 <- 1/nY + sigmac2 + 1/nX * (alpha + sum(gamma*delta))^2
    SC <- (1/nZ + sigmaC2) * Sigma + 1/nX * gamma %*% t(gamma)
    
    list(alpha=alpha, gamma=gamma, delta=delta, 
         C_hat=C_hat, B_hat=B_hat, beta_hat=beta_hat, c_hat=c_hat, b_hat=b_hat, 
         Sigma=Sigma, invSigma=invSigma, logdetSigma=logdetSigma, sb2=sb2, sc2=sc2, SC=SC,
         k=k, l=l, m=m, nX=nX, nZ=nZ, nY=nY, sigmaC2=sigmaC2, sigmac2=sigmac2, sigmab2=sigmab2)
}

test_that("estimate_Sigma is full rank", {
    rank_Sigma <- function(n, k) {
        med_hat <- matrix(rnorm(n*k), nrow = n, ncol = k)
        Sigma <- estimate_Sigma(med_hat)
        
        Matrix::rankMatrix(Sigma)[1]
    }
    
    # n > k
    expect_equal(rank_Sigma(3, 2), 2)
    expect_equal(rank_Sigma(4, 2), 2)
    
    # n <= k
    expect_warning(rank_Sigma(2, 2), "Estimated Sigma is not full rank. Applying shrinkage")
    expect_warning(rank_Sigma(2, 3), "Estimated Sigma is not full rank. Applying shrinkage")
    expect_equal(suppressWarnings(rank_Sigma(2, 2)), 2)
    expect_equal(suppressWarnings(rank_Sigma(2, 3)), 3)
    
    # n <= 1
    expect_warning(rank_Sigma(1, 1), "The number of observations are <= 1.")
    expect_warning(rank_Sigma(0, 1), "The number of observations are <= 1.")
    expect_equal(suppressWarnings(rank_Sigma(1, 1)), 1)
    expect_equal(suppressWarnings(rank_Sigma(0, 1)), 1)
    
    # k == 0
    expect_error(rank_Sigma(0, 0), "Estimating Sigma requires a positive number of mediators")
})


test_that("ll_original_b is correctly simplified", {
    test_ll_original_b <- function(sigmab2) {
        dat <- generate_test_data(sigmab2 = sigmab2)
        
        expect_equal(ll_original_b(dat$delta, dat$B_hat, dat$b_hat, dat$nZ, dat$nY, dat$sigmab2, dat$Sigma), 
                     mvtnorm::dmvnorm(dat$b_hat, mean = dat$B_hat %*% dat$delta, sigma = diag(rep(dat$sb2, dat$l)), log = TRUE))
    }
    
    test_ll_original_b(sigmab2 = 1e-3)
    test_ll_original_b(sigmab2 = 0)
})

test_that("ll_original_c is correctly simplified", {
    test_ll_original_c <- function(sigmac2) {
        dat <- generate_test_data(sigmac2 = sigmac2)
        
        expect_equal(ll_original_c(dat$gamma, dat$delta, dat$alpha, dat$beta_hat, dat$c_hat, dat$nX, dat$nY, dat$sigmac2),
                     mvtnorm::dmvnorm(dat$c_hat, mean = dat$beta_hat * (dat$alpha + sum(dat$gamma*dat$delta)), sigma = diag(rep(dat$sc2, dat$m)), log = TRUE))
    }
    
    test_ll_original_c(sigmac2 = 1e-4)
    test_ll_original_c(sigmac2 = 0)
})

test_that("ll_original_C is correctly simplified", {
    test_ll_original_C <- function(sigmaC2) {
        dat <- generate_test_data(sigmaC2 = sigmaC2)
        
        expect_equal(ll_original_C(dat$gamma, dat$beta_hat, dat$C_hat, dat$k, dat$nX, dat$nZ, dat$sigmaC2, dat$invSigma, dat$logdetSigma),
                     mvtnorm::dmvnorm(c(dat$C_hat), mean = c(dat$beta_hat %*% t(dat$gamma)), sigma = dat$SC %x% diag(dat$m), log = TRUE))
    }
    
    test_ll_original_C(sigmaC2 = 1e-5)
    test_ll_original_C(sigmaC2 = 0)
})

test_that("ll_original works correctly", {
    test_ll_original <- function(sigmaC2, sigmac2, sigmab2, pleiotropy) {
        dat <- suppressWarnings(generate_test_data(sigmaC2 = sigmaC2, sigmac2 = sigmac2, sigmab2 = sigmab2))
        
        x = c(dat$gamma, dat$delta, dat$alpha, dat$sigmaC2, dat$sigmac2, dat$sigmab2)
        mean <- c(dat$B_hat %*% dat$delta, dat$beta_hat * (dat$alpha + sum(dat$gamma*dat$delta)), dat$beta_hat %*% t(dat$gamma))
        sigma <- Matrix::bdiag(diag(rep(dat$sb2, dat$l)), diag(rep(dat$sc2, dat$m)), dat$SC %x% diag(dat$m)) %>%
            as.matrix()
        
        expect_equal(ll_original(x=x, pleiotropy=pleiotropy, dots=dat),
                     mvtnorm::dmvnorm(c(dat$b_hat, dat$c_hat, dat$C_hat), mean = mean, sigma = sigma, log = TRUE))
    }
    
    test_ll_original(sigmaC2 = 1e-5, sigmac2 = 1e-4, sigmab2 = 1e-3, pleiotropy = TRUE)
    test_ll_original(sigmaC2 = 0, sigmac2 = 0, sigmab2 = 0, pleiotropy = FALSE)
})

test_that("ll_original_SE_hessian works correctly", {
    generate_Hessian <- function(k, n = 10) {
        neg_hessian <- matrix(rnorm(n*k), nrow = n, ncol = k) %>%
            estimate_Sigma() %>%
            solve()
        
        return(-neg_hessian)
    }
    
    expect_SE_list <- function(SEs) {
        expect_type(SEs, "list")
        expect_vector(unlist(SEs), ptype = numeric(), size = 2)
        expect_named(SEs, c("total", "alpha"))
    }
    
    # Default
    SEs <- ll_original_SE_hessian(generate_Hessian(5), 2)
    expect_SE_list(SEs)
    expect_gt(SEs$total, 0)
    expect_gt(SEs$alpha, 0)
    
    # Hessian is NULL
    SEs <- ll_original_SE_hessian(NULL, 2) 
    expect_SE_list(SEs)
    expect_equal(SEs$total, NA_real_)
    expect_equal(SEs$alpha, NA_real_)

    # No mediators
    SEs <- ll_original_SE_hessian(generate_Hessian(1), 0)
    expect_SE_list(SEs)
    expect_equal(SEs$total, NA_real_)
    expect_gt(SEs$alpha, 0)
        
    # Hessian is not compatible with the number of mediators
    expect_error(ll_original_SE_hessian(generate_Hessian(5), 3))
    expect_warning(ll_original_SE_hessian(generate_Hessian(5), 1))
    
    # The covariance matrix approximation has negative values on the diagonal
    hessian <- matrix(c(-0.03, -0.3, 0.3, -0.3, -99, 2.5, 0.3, 2.5, -2.5), nrow = 3, ncol = 3)
    SEs <- expect_warning(ll_original_SE_hessian(hessian, 1), "The covariance matrix approximation has negative values on the diagonal. Returning NaNs.")
    expect_SE_list(SEs)
    expect_equal(SEs$total, NaN)
    expect_equal(SEs$alpha, NaN)
})

test_that("ll_original_SE_LRT works correctly", {
    ll_original_SE_data <- function(pleiotropy, seed = 1) {
        set.seed(seed)
        dat <- suppressWarnings(generate_test_data())
        
        # Get the ML solution
        opt <- suppressWarnings(ll_approach_original(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, 
                                                     dat$nY, dat$Sigma, dat$gamma, dat$delta, dat$alpha, pleiotropy))
        
        list(opt = opt, dots = dat)
    }
    
    expect_SE_list <- function(SEs) {
        expect_type(SEs, "list")
        expect_vector(unlist(SEs), ptype = numeric(), size = 2)
        expect_named(SEs, c("total", "alpha"))
    }
    
    expect_NaN_alpha <- function(opt, dots) {
        opt$alpha <- -opt$alpha
        SEs <- expect_warning(ll_original_SE_LRT(opt$gamma, opt$delta, opt$alpha, opt$sigmac2, dots), "Estimating SE of the direct effect")
        expect_SE_list(SEs)
        expect_equal(SEs$total, NaN)
        expect_equal(SEs$alpha, NaN)
    }
    
    expect_NaN_total <- function(opt, dots) {
        opt$gamma <- rep(0, length(opt$gamma))
        opt$delta <- rep(0, length(opt$delta))
        SEs <- expect_warning(ll_original_SE_LRT(opt$gamma, opt$delta, opt$alpha, opt$sigmac2, dots), "Estimating SE of the total indirect effect")
        expect_SE_list(SEs)
        expect_equal(SEs$total, NaN)
        expect_gt(SEs$alpha, 0)
    }
    
    expect_SEs <- function(opt, dots) {
        SEs <- ll_original_SE_LRT(opt$gamma, opt$delta, opt$alpha, opt$sigmac2, dots)
        
        expect_SE_list(SEs)
        expect_gt(SEs$alpha, 0)
        expect_gt(SEs$total, 0)
        
        expect_NaN_alpha(opt, dots)
        expect_NaN_total(opt, dots)
    }
    
    dat_SE <- ll_original_SE_data(FALSE)
    expect_SEs(dat_SE$opt, dat_SE$dots)
    
    dat_SE <- ll_original_SE_data(TRUE)
    expect_SEs(dat_SE$opt, dat_SE$dots)
})

test_that("ll_approach_original works correctly", {
    expect_opt_list <- function(opt) {
        expect_type(opt, "list")
        expect_length(opt, 13)
        expect_named(opt, c("total", "alpha", "gamma", "delta", "SE_total", "SE_alpha", "sigmag2", 
                            "sigmad2", "sigmaC2", "sigmac2", "sigmab2", "convergence", "optim_value"))
        expect_vector(unlist(opt), ptype = numeric())
        
        expect_equal(opt$convergence, 0)
        expect_true(!is.na(opt$total))
        expect_true(!is.na(opt$alpha))
        expect_true(all(!is.na(opt$gamma)))
        expect_true(all(!is.na(opt$delta)))
        expect_true(is.nan(opt$SE_total) || opt$SE_total > 0)
        expect_true(is.nan(opt$SE_alpha) || opt$SE_alpha > 0)
        expect_true(!is.na(opt$optim_value))
    }
    
    expect_no_pleiotropy <- function(opt) {
        expect_equal(opt$sigmaC2, 0)
        expect_equal(opt$sigmac2, 0)
        expect_equal(opt$sigmab2, 0)
    }
    
    expect_pleiotropy <- function(opt) {
        expect_gte(opt$sigmaC2, 0)
        expect_gte(opt$sigmac2, 0)
        expect_gte(opt$sigmab2, 0)
    }
    
    expect_no_variances <- function(opt) {
        expect_equal(opt$sigmag2, NA_real_)
        expect_equal(opt$sigmad2, NA_real_)
    }
    
    expect_variances <- function(opt) {
        expect_gte(opt$sigmag2, 0)
        expect_gte(opt$sigmad2, 0)
    }
    
    # With pleiotropy, k > 1
    dat <- suppressWarnings(generate_test_data(sigmaC2 = 1e-5, sigmac2 = 1e-4, sigmab2 = 1e-3))
    opt <- suppressWarnings(ll_approach_original(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, 
                                                 dat$nY, dat$Sigma, dat$gamma, dat$delta, dat$alpha, pleiotropy = TRUE))
    expect_opt_list(opt)
    expect_pleiotropy(opt)
    expect_variances(opt)
    
    # With pleiotropy, k == 1
    dat <- suppressWarnings(generate_test_data(sigmaC2 = 1e-5, sigmac2 = 1e-4, sigmab2 = 1e-3, k = 1))
    opt <- suppressWarnings(ll_approach_original(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, 
                                                 dat$nY, dat$Sigma, dat$gamma, dat$delta, dat$alpha, pleiotropy = TRUE))
    expect_opt_list(opt)
    expect_pleiotropy(opt)
    expect_no_variances(opt)
    
    # Without pleiotropy, k > 1
    dat <- suppressWarnings(generate_test_data(sigmaC2 = 0, sigmac2 = 0, sigmab2 = 0))
    opt <- suppressWarnings(ll_approach_original(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, 
                                                 dat$nY, dat$Sigma, dat$gamma, dat$delta, dat$alpha, pleiotropy = FALSE))
    expect_opt_list(opt)
    expect_no_pleiotropy(opt)
    expect_variances(opt)
    
    # Without pleiotropy, k == 1
    dat <- suppressWarnings(generate_test_data(sigmaC2 = 0, sigmac2 = 0, sigmab2 = 0, k = 1))
    opt <- suppressWarnings(ll_approach_original(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, 
                                                 dat$nY, dat$Sigma, dat$gamma, dat$delta, dat$alpha, pleiotropy = FALSE))
    expect_opt_list(opt)
    expect_no_pleiotropy(opt)
    expect_no_variances(opt)
})
