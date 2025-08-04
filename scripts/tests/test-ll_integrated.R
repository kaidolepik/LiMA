library(testthat)

generate_test_data <- function(k = 2, l = 3, m = 4, nX = 300000, nZ = 10000, nY = 350000, 
                               sigmaC2 = 1e-5, sigmac2 = 1e-4, sigmab2 = 1e-3) {
    Sigma <- matrix(rnorm(l*k), nrow = l, ncol = k) %>%
        estimate_Sigma()
    
    alpha <- rnorm(1)
    gamma <- rnorm(k)
    delta <- rnorm(k)
    C_hat <- mvtnorm::rmvnorm(m, sigma = (1/nZ + sigmaC2) * Sigma)
    B_hat <- mvtnorm::rmvnorm(l, sigma = 1/nZ * Sigma)
    beta_hat <- rnorm(m, sd = sqrt(1/nX))
    beta_med_hat <- rnorm(l, sd = sqrt(1/nX))
    c_hat <- rnorm(m, sd = sqrt(1/nY + sigmac2))
    b_hat <- rnorm(l, sd = sqrt(1/nY + sigmab2))
    
    sc2 <- 1/nY + sigmac2 + 1/nX * (alpha + sum(gamma*delta))^2
    Sigma <- diag(k) # ll_integrated assumes identity Sigma
    
    sigmag2 <- ifelse(k > 1, var(gamma), 0.2)
    sigmad2 <- ifelse(k > 1, var(delta), 0.3)
    sigmagd <- ifelse(k > 1, cov(gamma, delta), 0.1)
    if (sigmagd^2 > sigmag2*sigmad2)
        sigmagd <- sign(sigmagd) * (sqrt(sigmag2 * sigmad2) - 1e-8)
    
    sigmaw2 <- k * (sigmag2*sigmad2 + sigmagd^2)
    sZ2 <- 1/nZ + sigmaC2
    sC2 <- sZ2 + 1/nX*sigmag2
    mu_sb2 <- 1/nY + sigmab2 + 1/nZ * sigmad2 * sum(diag(Sigma))
    mu_sc2 <- 1/nY + sigmac2 + 1/nX * (sigmaw2 + (alpha + k*sigmagd)^2)
    xib <- sigmag2 / (sC2 + sigmag2 * sum(beta_hat^2))
    xi1 <- sigmad2 - sigmagd^2/sC2 * sum(beta_hat^2) * (1 - xib * sum(beta_hat^2))
    
    # For just_alpha method
    total <- alpha + k*sigmagd
    
    list(alpha=alpha, gamma=gamma, delta=delta, C_hat=C_hat, B_hat=B_hat, beta_hat=beta_hat, c_hat=c_hat, b_hat=b_hat, sc2=sc2,
         Sigma=Sigma, sigmag2=sigmag2, sigmad2=sigmad2, sigmagd=sigmagd, sigmaw2=sigmaw2, sZ2=sZ2, sC2=sC2, mu_sb2=mu_sb2, mu_sc2=mu_sc2, xib=xib, xi1=xi1,
         total=total, k=k, l=l, m=m, nX=nX, nZ=nZ, nY=nY, sigmaC2=sigmaC2, sigmac2=sigmac2, sigmab2=sigmab2)
}

bC_Sigma <- function(dat) {
    A <- diag(rep(dat$mu_sb2, dat$l)) + dat$sigmad2 * dat$B_hat %*% t(dat$B_hat)
    B <- dat$sigmagd * (dat$B_hat %x% t(dat$beta_hat))
    C <- dat$sigmagd * (t(dat$B_hat) %x% dat$beta_hat)
    D <- dat$Sigma %x% diag(rep(dat$sZ2, dat$m)) + diag(dat$k) %x% (dat$sigmag2 * (diag(1/dat$nX, dat$m) + dat$beta_hat %*% t(dat$beta_hat)))
    
    Sigma <- rbind(cbind(A, B), cbind(C, D))
    
    return(Sigma)
}

test_that("get_L11 returns the correct matrix", {
    expect_L11 <- function(k, l) {
        dat <- generate_test_data(k = k, l = l)
        
        precision_mat <- bC_Sigma(dat) %>%
            solve()
        
        expect_equal(get_L11(dat$B_hat, dat$mu_sb2, dat$xi1, dat$k, dat$l), precision_mat[1:dat$l, 1:dat$l])
    }
    
    expect_L11(3, 4) # l <= 1.5k
    expect_L11(3, 5) # l > 1.5k
})

test_that("ll_integrated_c is correctly simplified", {
    test_ll_integrated_c <- function(sigmac2) {
        dat <- generate_test_data(sigmac2 = sigmac2)
        
        mean <- dat$beta_hat * (dat$alpha + dat$k*dat$sigmagd)
        sigma <- diag(rep(dat$mu_sc2, dat$m)) + dat$sigmaw2 * dat$beta_hat %*% t(dat$beta_hat)
        
        expect_equal(ll_integrated_c(dat$alpha, dat$sigmagd, dat$sigmag2, dat$sigmad2, dat$beta_hat, dat$c_hat, dat$k, dat$m, dat$nX, dat$nY, dat$sigmac2),
                     mvtnorm::dmvnorm(dat$c_hat, mean = mean, sigma = sigma, log = TRUE))
    }
    
    test_ll_integrated_c(sigmac2 = 1e-4)
    test_ll_integrated_c(sigmac2 = 0)
})

test_that("ll_integrated_bC is correctly simplified", {
    ll_integrated_bC_simplified <- function(dat) {
        ll_integrated_bC(dat$sigmagd, dat$sigmag2, dat$sigmad2, dat$beta_hat, dat$C_hat, dat$B_hat, dat$b_hat, dat$k, dat$l, dat$m, dat$nX, dat$nZ, dat$nY, dat$sigmaC2, dat$sigmab2)
    }
    
    ll_integrated_bC_full <- function(dat) {
        mvtnorm::dmvnorm(c(dat$b_hat, dat$C_hat), mean = rep(0, dat$l + dat$m*dat$k), sigma = bC_Sigma(dat), log = TRUE)
    }
    
    # Pleiotropy
    dat <- generate_test_data(sigmaC2 = 1e-5, sigmab2 = 1e-3)
    expect_equal(ll_integrated_bC_simplified(dat), ll_integrated_bC_full(dat))
    
    dat$sigmagd <- sign(dat$sigmagd) * sqrt(dat$sigmag2*dat$sigmad2 + 1e-10)
    expect_equal(ll_integrated_bC_simplified(dat), -sqrt(conf$inf))
    
    dat$sigmagd <- -dat$sigmagd
    expect_equal(ll_integrated_bC_simplified(dat), -sqrt(conf$inf))
    
    # No pleiotropy
    dat <- generate_test_data(sigmaC2 = 0, sigmab2 = 0)
    expect_equal(ll_integrated_bC_simplified(dat), ll_integrated_bC_full(dat))
    
    dat$sigmagd <- sign(dat$sigmagd) * sqrt(dat$sigmag2*dat$sigmad2 + 1e-10)
    expect_equal(ll_integrated_bC_simplified(dat), -sqrt(conf$inf))
    
    dat$sigmagd <- -dat$sigmagd
    expect_equal(ll_integrated_bC_simplified(dat), -sqrt(conf$inf))
})

test_that("ll_integrated_b is correctly simplified", {
    ll_integrated_b_simplified <- function(dat, x, pleiotropy) {
        ll_integrated_b(x, dat$B_hat, dat$b_hat, dat$k, dat$l, dat$nZ, dat$nY, dat$sigmab2, pleiotropy)
    }
    
    ll_integrated_b_full <- function(dat) {
        mvtnorm::dmvnorm(dat$b_hat, mean = rep(0, dat$l), sigma = bC_Sigma(dat)[1:dat$l, 1:dat$l], log = TRUE)
    }
    
    # Pleiotropy
    dat <- generate_test_data(sigmab2 = 1e-4)
    expect_equal(ll_integrated_b_simplified(dat, c(dat$sigmad2, dat$sigmab2), TRUE), ll_integrated_b_full(dat))
    expect_equal(ll_integrated_b_simplified(dat, c(dat$sigmad2), TRUE), NA_real_)
    expect_equal(ll_integrated_b_simplified(dat, c(dat$sigmad2, 1e-3), FALSE), ll_integrated_b_full(dat))
    expect_true(ll_integrated_b_simplified(dat, c(dat$sigmad2, 1e-3), TRUE) != ll_integrated_b_full(dat))
    
    # No pleiotropy
    dat <- generate_test_data(sigmab2 = 0)
    expect_equal(ll_integrated_b_simplified(dat, c(dat$sigmad2), FALSE), ll_integrated_b_full(dat))
    expect_equal(ll_integrated_b_simplified(dat, c(dat$sigmad2, dat$sigmab2), FALSE), ll_integrated_b_full(dat))
    expect_equal(ll_integrated_b_simplified(dat, c(dat$sigmad2, 1e-3), FALSE), ll_integrated_b_full(dat))
    expect_true(ll_integrated_b_simplified(dat, c(dat$sigmad2, 1e-3), TRUE) != ll_integrated_b_full(dat))
})

test_that("ll_integrated_C is correctly simplified", {
    ll_integrated_C_simplified <- function(dat, x, pleiotropy) {
        ll_integrated_C(x, dat$beta_hat, dat$C_hat, dat$k, dat$m, dat$nX, dat$nZ, dat$sigmaC2, pleiotropy)
    }
    
    ll_integrated_C_full <- function(dat) {
        mvtnorm::dmvnorm(c(dat$C_hat), mean = rep(0, dat$m*dat$k), sigma = bC_Sigma(dat)[(dat$l+1):(dat$l + dat$m*dat$k), (dat$l+1):(dat$l + dat$m*dat$k)], log = TRUE)
    }
    
    # Pleiotropy
    dat <- generate_test_data(sigmaC2 = 1e-5)
    expect_equal(ll_integrated_C_simplified(dat, c(dat$sigmag2, dat$sigmaC2), TRUE), ll_integrated_C_full(dat))
    expect_equal(ll_integrated_C_simplified(dat, c(dat$sigmag2), TRUE), NA_real_)
    expect_equal(ll_integrated_C_simplified(dat, c(dat$sigmag2, 1e-3), FALSE), ll_integrated_C_full(dat))
    expect_false(ll_integrated_C_simplified(dat, c(dat$sigmag2, 1e-3), TRUE) == ll_integrated_C_full(dat))
    
    # No pleiotropy
    dat <- generate_test_data(sigmaC2 = 0)
    expect_equal(ll_integrated_C_simplified(dat, c(dat$sigmag2), FALSE), ll_integrated_C_full(dat))
    expect_equal(ll_integrated_C_simplified(dat, c(dat$sigmag2, dat$sigmaC2), FALSE), ll_integrated_C_full(dat))
    expect_equal(ll_integrated_C_simplified(dat, c(dat$sigmag2, 1e-3), FALSE), ll_integrated_C_full(dat))
    expect_false(ll_integrated_C_simplified(dat, c(dat$sigmag2, 1e-3), TRUE) == ll_integrated_C_full(dat))
})

test_that("ll_integrated works correctly", {
    ll_integrated_simplified <- function(dat, x, ll_type, pleiotropy) {
        ll_integrated(x = x, ll_type = ll_type, pleiotropy = pleiotropy, dots = dat)
    }
    
    ll_integrated_full <- function(dat) {
        mean <- c(dat$beta_hat * (dat$alpha + dat$k*dat$sigmagd), rep(0, dat$l + dat$m*dat$k))
        sigma <- Matrix::bdiag(diag(rep(dat$mu_sc2, dat$m)) + dat$sigmaw2 * dat$beta_hat %*% t(dat$beta_hat), bC_Sigma(dat)) %>%
            as.matrix()
        
        mvtnorm::dmvnorm(c(dat$c_hat, dat$b_hat, dat$C_hat), mean = mean, sigma = sigma, log = TRUE)
    }
    
    expect_correct_flow <- function(dat, ll_type) {
        sigmagd <- ifelse(ll_type == "free", dat$sigmagd / sqrt(dat$sigmag2 * dat$sigmad2), dat$sigmagd)
            
        expect_equal(ll_integrated_simplified(dat, c(dat$alpha, sigmagd, dat$sigmag2, dat$sigmad2, dat$sigmaC2, dat$sigmac2, dat$sigmab2), ll_type, TRUE), ll_integrated_full(dat))
        expect_equal(ll_integrated_simplified(dat, c(dat$alpha, sigmagd, dat$sigmag2, dat$sigmad2, dat$sigmaC2, dat$sigmac2, dat$sigmab2), ll_type, FALSE), ll_integrated_full(dat))
        expect_equal(ll_integrated_simplified(dat, c(dat$alpha, sigmagd, dat$sigmag2, dat$sigmad2, 0.1, 0.1, 0.1), ll_type, FALSE), ll_integrated_full(dat))
        expect_false(ll_integrated_simplified(dat, c(dat$alpha, sigmagd, dat$sigmag2, dat$sigmad2, 0.1, 0.1, 0.1), ll_type, TRUE) == ll_integrated_full(dat))
    }
    
    test_ll_integrated <- function(sigmaC2, sigmac2, sigmab2) {
        dat <- suppressWarnings(generate_test_data(sigmaC2 = sigmaC2, sigmac2 = sigmac2, sigmab2 = sigmab2))
        
        expect_correct_flow(dat, "fixed")
        expect_equal(ll_integrated_simplified(dat, c(dat$alpha, dat$sigmagd, 0, 0, dat$sigmaC2, dat$sigmac2, dat$sigmab2), "fixed", TRUE), ll_integrated_full(dat))
        expect_equal(ll_integrated_simplified(dat, c(dat$alpha, dat$sigmagd, 0, 0, dat$sigmaC2, dat$sigmac2, dat$sigmab2), "fixed", FALSE), ll_integrated_full(dat))
        
        expect_correct_flow(dat, "free")
        expect_false(ll_integrated_simplified(dat, c(dat$alpha, 0.1, 0.1, 0.1, dat$sigmaC2, dat$sigmac2, dat$sigmab2), "free", TRUE) == ll_integrated_full(dat))
        expect_false(ll_integrated_simplified(dat, c(dat$alpha, 0.1, 0.1, 0.1, dat$sigmaC2, dat$sigmac2, dat$sigmab2), "free", FALSE) == ll_integrated_full(dat))
        
        expect_correct_flow(dat, "just_alpha")
        expect_equal(ll_integrated_simplified(dat, c(dat$alpha, 0, 0, 0, dat$sigmaC2, dat$sigmac2, dat$sigmab2), "just_alpha", TRUE), ll_integrated_full(dat))
        expect_equal(ll_integrated_simplified(dat, c(dat$alpha, 0, 0, 0, dat$sigmaC2, dat$sigmac2, dat$sigmab2), "just_alpha", FALSE), ll_integrated_full(dat))
    }
    
    set.seed(1) # Numerical issues can sometimes cause just_alpha method to have sigmagd^2 > sigmag2*sigmad2
    test_ll_integrated(sigmaC2 = 1e-5, sigmac2 = 1e-4, sigmab2 = 1e-3)
    test_ll_integrated(sigmaC2 = 0, sigmac2 = 0, sigmab2 = 0)
})

test_that("ll_integrated_SE_LRT works correctly", {
    ll_integrated_SE_data <- function(ll_type = "fixed", pleiotropy = TRUE, seed = 1) {
        set.seed(seed)
        dat <- suppressWarnings(generate_test_data())
        
        # Get the ML solution
        opt <- suppressWarnings(ll_approach_integrated(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, dat$nY, 
                                                       dat$total, dat$alpha, dat$sigmag2, dat$sigmad2, ll_type, pleiotropy))
        par <- c(opt$alpha, (opt$total - opt$alpha) / dat$k, opt$sigmag2, opt$sigmad2, opt$sigmaC2, opt$sigmac2, opt$sigmab2)
        dat$total <- opt$total
        
        list(par = par, dots = dat)
    }
    
    expect_SE_list <- function(SEs) {
        expect_type(SEs, "list")
        expect_vector(unlist(SEs), ptype = numeric(), size = 2)
        expect_named(SEs, c("total", "alpha"))
    }
    
    expect_NaN_alpha <- function(dat_SE, ll_type, pleiotropy) {
        dat_SE$par[1] <- -dat_SE$par[1]
        SEs <- expect_warning(ll_integrated_SE_LRT(dat_SE$par, ll_type, pleiotropy, dat_SE$dots), "Estimating SE of the direct effect")
        expect_SE_list(SEs)
        expect_equal(SEs$total, NaN)
        expect_equal(SEs$alpha, NaN)
    }
    
    expect_NaN_total <- function(dat_SE, ll_type, pleiotropy) {
        dat_SE$par[2] <- -dat_SE$par[2]
        SEs <- expect_warning(ll_integrated_SE_LRT(dat_SE$par, "fixed", TRUE, dat_SE$dots), "Estimating SE of the total indirect effect")
        expect_SE_list(SEs)
        expect_equal(SEs$total, NaN)
        expect_gt(SEs$alpha, 0)
    }
    
    expect_SEs <- function(dat_SE, ll_type, pleiotropy) {
        SEs <- ll_integrated_SE_LRT(dat_SE$par, ll_type, pleiotropy, dat_SE$dots)
        
        expect_SE_list(SEs)
        expect_gt(SEs$alpha, 0)
        expect_NaN_alpha(dat_SE, ll_type, pleiotropy)
        
        if (ll_type == "just_alpha")
            expect_equal(SEs$total, NA_real_)
        else {
            expect_gt(SEs$total, 0)
            expect_NaN_total(dat_SE, ll_type, pleiotropy)
        }
    }
    
    dat_SE <- ll_integrated_SE_data()
    
    expect_SEs(dat_SE, "fixed", TRUE)
    expect_SEs(dat_SE, "fixed", FALSE)
    expect_SEs(dat_SE, "free", TRUE)
    expect_SEs(dat_SE, "free", FALSE)
    expect_SEs(dat_SE, "just_alpha", TRUE)
    expect_SEs(dat_SE, "just_alpha", FALSE)
})

test_that("sigmag2_approach_integrated works correctly", {
    expect_opt_list <- function(opt) {
        expect_type(opt, "list")
        expect_length(opt, 2)
        expect_named(opt, c("sigmag2", "convergence"))
        
        expect_gt(opt$sigmag2, 0)
        expect_equal(opt$convergence, 0)
    }
    
    set.seed(1)
    dat <- suppressWarnings(generate_test_data())
    
    expect_opt_list(sigmag2_approach_integrated(dat$beta_hat, dat$C_hat, dat$k, dat$m, dat$nX, dat$nZ, dat$sigmaC2, TRUE))
    expect_opt_list(sigmag2_approach_integrated(dat$beta_hat, dat$C_hat, dat$k, dat$m, dat$nX, dat$nZ, dat$sigmaC2, FALSE))
})

test_that("sigmad2_approach_integrated works correctly", {
    expect_opt_list <- function(opt) {
        expect_type(opt, "list")
        expect_length(opt, 2)
        expect_named(opt, c("sigmad2", "convergence"))
        
        expect_gt(opt$sigmad2, 0)
        expect_equal(opt$convergence, 0)
    }
    
    set.seed(2)
    dat <- suppressWarnings(generate_test_data())
    
    expect_opt_list(sigmad2_approach_integrated(dat$B_hat, dat$b_hat, dat$k, dat$l, dat$nZ, dat$nY, dat$sigmab2, TRUE))
    expect_opt_list(sigmad2_approach_integrated(dat$B_hat, dat$b_hat, dat$k, dat$l, dat$nZ, dat$nY, dat$sigmab2, FALSE))
})

test_that("ll_approach_integrated works correctly", {
    expect_opt_list <- function(opt, variances_convergence_value, ll_type) {
        expect_type(opt, "list")
        expect_length(opt, 16)
        expect_named(opt, c("convergence", "optim_value", "is_seed_OK", "sigmag2_convergence", "sigmad2_convergence", "sigmag2_bound", "sigmad2_bound",
                            "total", "alpha", "SE_total", "SE_alpha", "sigmag2", "sigmad2", "sigmaC2", "sigmac2", "sigmab2"))

        expect_equal(opt$convergence, 0)
        expect_true(!is.na(opt$optim_value))
        expect_true(is.logical(opt$is_seed_OK))
        expect_equal(opt$sigmag2_convergence, variances_convergence_value)
        expect_equal(opt$sigmad2_convergence, variances_convergence_value)
        expect_true(is.logical(opt$sigmag2_bound))
        expect_true(is.logical(opt$sigmad2_bound))
        
        expect_true(!is.na(opt$total))
        expect_true(!is.na(opt$alpha))
        expect_true(is.nan(opt$SE_total) || (is.na(opt$SE_total) && ll_type == "just_alpha") || opt$SE_total > 0)
        expect_true(is.nan(opt$SE_alpha) || (is.na(opt$SE_total) && ll_type == "just_alpha") || opt$SE_alpha > 0)
        
        expect_gte(opt$sigmag2, 0)
        expect_gte(opt$sigmad2, 0)
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
    
    set.seed(10)
    dat <- suppressWarnings(generate_test_data(sigmaC2 = 1e-5, sigmac2 = 1e-4, sigmab2 = 1e-3, k = 1))
    
    # With pleiotropy, fixed    
    opt <- suppressWarnings(ll_approach_integrated(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, dat$nY, 
                                                   dat$total, dat$alpha, NULL, NULL, ll_type = "fixed", pleiotropy = TRUE))
    expect_opt_list(opt, 0, "fixed")
    expect_pleiotropy(opt)
    
    opt <- suppressWarnings(ll_approach_integrated(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, dat$nY, 
                                                   dat$total, dat$alpha, dat$sigmag2, dat$sigmad2, ll_type = "fixed", pleiotropy = TRUE))
    expect_opt_list(opt, NA_real_, "fixed")
    expect_pleiotropy(opt)
    
    opt <- suppressWarnings(ll_approach_integrated(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, dat$nY, 
                                                   dat$total, dat$alpha, NULL, NULL, ll_type = "fixed", pleiotropy = FALSE))
    expect_opt_list(opt, 0, "fixed")
    expect_pleiotropy(opt)
    
    opt <- suppressWarnings(ll_approach_integrated(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, dat$nY, 
                                                   dat$total, dat$alpha, dat$sigmag2, dat$sigmad2, ll_type = "fixed", pleiotropy = FALSE))
    expect_opt_list(opt, NA_real_, "fixed")
    expect_pleiotropy(opt)
    
    # With pleiotropy, free
    opt <- suppressWarnings(ll_approach_integrated(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, dat$nY, 
                                                   dat$total, dat$alpha, NULL, NULL, ll_type = "free", pleiotropy = TRUE))
    expect_opt_list(opt, NA_real_, "free")
    expect_pleiotropy(opt)
    
    opt <- suppressWarnings(ll_approach_integrated(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, dat$nY, 
                                                   dat$total, dat$alpha, NULL, NULL, ll_type = "free", pleiotropy = FALSE))
    expect_opt_list(opt, NA_real_, "free")
    expect_pleiotropy(opt)
    
    # With pleiotropy, just_alpha
    opt <- suppressWarnings(ll_approach_integrated(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, dat$nY, 
                                                   dat$total, dat$alpha, NULL, NULL, ll_type = "just_alpha", pleiotropy = TRUE))
    expect_opt_list(opt, 0, "just_alpha")
    expect_pleiotropy(opt)
    
    opt <- suppressWarnings(ll_approach_integrated(dat$beta_hat, dat$C_hat, dat$B_hat, dat$c_hat, dat$b_hat, dat$nX, dat$nZ, dat$nY, 
                                                   dat$total, dat$alpha, NULL, NULL, ll_type = "just_alpha", pleiotropy = FALSE))
    expect_opt_list(opt, 0, "just_alpha")
    expect_pleiotropy(opt)
})
