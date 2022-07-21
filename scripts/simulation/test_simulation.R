# a simple simulation workflow


# set up parameters
n = 100 # variants
m = 10000 # samples
p_dosage = 0.01
h2_beta = 0.00
h2_theta = 0.05
pi_beta = 0.1
pi_theta = 0.1

# generate matrix 
H1 <- matrix(rbinom(n*m, 1, p_dosage), nrow = n, ncol = m)
H2 <- matrix(rbinom(n*m, 1, p_dosage), nrow = n, ncol = m)
G <- H1+H2

# only keep entries with SD > 0
keep <- apply(G, 1, sd) > 0
G <- G[keep, ]
n <- nrow(G)
beta <- make_theta(n, h2_beta, pi = pi_beta)
beta_sign <- sign(beta) # note that this is zero if beta is zero
beta_sign[beta_sign == 0] <- 1

# normalize genotypes
G_norm <- (G - apply(G, 1, mean )) / apply(G, 1, sd)

# Only save genotype matrix when alternate alleles are present
G_alt_var <- G==2
G_norm_alt <- G_norm
G_norm_alt[G_alt_var] <- 0
theta <- abs(make_theta(n, h2_theta, pi = pi_theta)) * beta_sign


# genetic effect contribution
y_no_noise_add <- as.vector(beta %*% G_norm) 
y_no_noise_dom <- as.vector(theta %*% G_norm_alt) 
var(y_no_noise_add)
var(y_no_noise_dom)

beta_theta_ratio <- theta / beta
beta_theta_ratio[beta_theta_ratio == Inf] <- NA
hist(beta_theta_ratio)

# final phenotype with noise
y_no_noise <- y_no_noise_add + y_no_noise_dom
var(y_no_noise)

y <- y_no_noise + y_no_noise_alt + rnorm(m, sd=sqrt(1-h2_beta-h2_theta))
y <- y / sd(y)

# simple linear regression (gwas)
gwas <- do.call(rbind, lapply(1:n, function(i){
  fit <- summary(lm(y ~ G_norm[i,]))
  fit <- as.data.frame(t(as.matrix(fit$coefficients[2,])))
  colnames(fit) <- c("est","std","t","p")
  fit$logp <- round(-log10(fit$p),3)
  fit$ac <- AC[i]
  fit$theta <- theta[i]
  fit$beta <- beta[i]
  fit$ratio <- beta_theta_ratio[i]
  fit$significant <- fit$p < (0.05/n)
  fit$sig <- ifelse(fit$significant,'*','')
  return(fit)
}))

lims <- max(abs(gwas$beta), abs(gwas$theta)) * 1.2

params <- paste0("Beta-h2/P: ", h2_beta, "/", pi_beta, ". Theta-h2/P: ", h2_theta, "/", pi_theta)
ggplot(gwas, aes(x=beta, y=theta, color=significant)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-lims, lims) +
  ylim(-lims, lims) +
  xlab("Beta (Additive effects)") +
  ylab("Theta (Non-additive effects)") +
  labs(color = "Bonf. Sig. Marker") +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
        ) +
  ggtitle("Simulation", params)
