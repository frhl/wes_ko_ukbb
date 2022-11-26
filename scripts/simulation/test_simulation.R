# a simple simulation workflow

require(ggplot2)
require(cowplot)

make_theta <- function(M, h2, pi = NULL){
  pi <- ifelse(!is.null(pi), pi, 1)
  stopifnot((pi > 0 & pi <= 1))
  stopifnot((h2 >= 0 & h2 <= 1))
  causal <- rbinom(M, size = 1, prob = pi)
  effect <- rnorm(M, mean = 0, sd = sqrt(h2/(M*pi)))
  return(causal * effect)
}

rescale_variance <- function(X, X_variance, Y_variance) {
  # Rescale variance
  X_sd = X_variance^(1/2)
  Y_sd = Y_variance^(1/2)
  K = Y_sd / X_sd
  Y = K*X
  return(Y)
}

annotate_effects <- function(beta, theta){
  both_sig <- (abs(beta) > 0) & (abs(theta) > 0)
  beta_sig <- (abs(beta) > 0)
  theta_sig <- (abs(theta) > 0)
  effects <-ifelse(both_sig, "both", ifelse(beta_sig, "beta", ifelse(theta_sig, "theta", NA)))
  return(effects)
}

#""" Make effect sizes for either infintesimal or spike and slab model. """
#M = mt.count()[0]
#pi_temp = 1 if pi == None else pi
#return(hl.rand_bool(pi_temp)*hl.rand_norm(0, hl.sqrt(h2/(M*pi_temp))))

graphics.off()

# set up parameters
n = 100 # variants
m = 10000 # samples
p_dosage = 0.01
h2_beta = 0.5
h2_theta = 0.5
h2 = 0.02
pi_beta = 0.1
pi_theta = 0.1

# generate matrix
H1 <- matrix(rbinom(n*m, 1, p_dosage), nrow = n, ncol = m)
H2 <- matrix(rbinom(n*m, 1, p_dosage), nrow = n, ncol = m)
G <- H1+H2
rownames(G) <- paste0("chr0:",1:nrow(G),"A:A")

# only keep entries with SD > 0
keep <- apply(G, 1, sd) > 0
G <- G[keep, ]
n <- nrow(G)
beta <- make_theta(n, h2_beta, pi = pi_beta)
beta_sign <- sign(beta) # note that this is zero if beta is zero
beta_sign[beta_sign == 0] <- 1
df_beta <- data.frame(rsid=rownames(G), beta = beta)

# normalize genotypes
G_norm <- (G - apply(G, 1, mean )) / apply(G, 1, sd)

# Only save genotype matrix when alternate alleles are present
G_only_one_allele <- G==1
G_alt <- G
G_alt[G_only_one_allele] <- 0

keep <- apply(G_alt, 1, sd) > 0
G_alt <- G_alt[keep, ]
G_norm_alt <- (G_alt - apply(G_alt, 1, mean )) / apply(G_alt, 1, sd)

#G_norm_alt[G_only_one_allele] <- 0
theta <- abs(make_theta(n, h2_theta, pi = pi_theta)) * beta_sign
theta <- theta[keep]
print(sum(theta > 0 ))

# get gene effects in single matrix
df_theta <- data.frame(rsid=rownames(G_alt), theta = theta)
df_effects <- merge(df_theta, df_beta, all = TRUE)

# genetic effect contribution
y_no_noise_add <- as.vector(beta %*% G_norm)
y_no_noise_rec <- as.vector(theta %*% G_norm_alt)

var(y_no_noise_add) # should approximate h2 on avg
var(y_no_noise_rec) # should approximate h2 on avg

# add up and rescale
y_no_noise <- y_no_noise_add + y_no_noise_rec
y_no_noise_rescaled <- rescale_variance(y_no_noise, var(y_no_noise), h2)
var(y_no_noise_rescaled)

# final phenotype with noise
noise <-rnorm(m, sd=sqrt(1-h2))
y <- y_no_noise_rescaled + noise
AC_add <- rowSums(G)
AC_rec <- rowSums(G_alt)

# simple linear regression (gwas)
gwas_add <- do.call(rbind, lapply(1:nrow(G_norm), function(i){
  fit <- summary(lm(y ~ G_norm[i,]))
  fit <- as.data.frame(t(as.matrix(fit$coefficients[2,])))
  fit$rsid <- rownames(G_norm)[i]
  colnames(fit) <- c("est","std","t","p", "rsid")
  fit <- fit[,c("rsid", "est","std","t","p")]
  fit$logp <- round(-log10(fit$p),3)
  fit$AC <- sum(G[i,])
  fit$significant <- fit$p < (0.05/n)
  fit$sig <- ifelse(fit$significant,'*','')
  return(fit)
}))

gwas_add <- merge(gwas_add, df_effects, all.x = TRUE, by = 'rsid')
gwas_add$effects <- annotate_effects(gwas_add$beta, gwas_add$theta)

# simple linear regression (gwas)
gwas_rec <- do.call(rbind, lapply(1:nrow(G_norm_alt), function(i){
  fit <- summary(lm(y ~ G_norm_alt[i,]))
  fit <- as.data.frame(t(as.matrix(fit$coefficients[2,])))
  fit$rsid <- rownames(G_norm_alt)[i]
  colnames(fit) <- c("est","std","t","p", "rsid")
  fit <- fit[,c("rsid", "est","std","t","p")]
  fit$logp <- round(-log10(fit$p),3)
  fit$AC <- sum(G_alt[i,])
  fit$significant <- fit$p < (0.05/n)
  fit$sig <- ifelse(fit$significant,'*','')
  return(fit)
}))

gwas_rec <- merge(gwas_rec, df_effects, all.x = TRUE, by = 'rsid')
gwas_rec$effects <- annotate_effects(gwas_rec$beta, gwas_rec$theta)




betas <- as.vector(na.omit(c(gwas_rec$beta, gwas_add$beta)))
thetas <- as.vector(na.omit(c(gwas_rec$theta, gwas_add$theta)))
xlims <- max(betas, thetas) * 1.01

ps <- c(gwas_add$logp, gwas_rec$logp)
ylims <- max(ps) * 1.01


params <- paste0("Beta-h2/P: ", h2_beta, "/", pi_beta, ". Theta-h2/P: ", h2_theta, "/", pi_theta)
p1 <- ggplot(gwas_add, aes(x=beta, y=logp, size=AC, color = effects)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-xlims, xlims) +
  ylim(0, ylims) +
  xlab("Beta") +
  ylab("-Log10(P-value)") +
  labs(color = "Bonf. Sig. Marker") +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
  ) + ggtitle("Beta (Additive)", params)

p2 <- ggplot(gwas_add, aes(x=theta, y=logp, size=AC, color = effects)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-xlims, xlims) +
  ylim(0, ylims) +
  xlab("Theta") +
  ylab("-Log10(P-value)") +
  labs(color = "Bonf. Sig. Marker") +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
  ) + ggtitle("Theta (Additive)", params)

p3 <- ggplot(gwas_rec, aes(x=beta, y=logp, size=AC, color = effects)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-xlims, xlims) +
  ylim(0, ylims) +
  xlab("Beta") +
  ylab("-Log10(P-value)") +
  labs(color = "Bonf. Sig. Marker") +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
  ) + ggtitle("Beta (Recessive)", params)

p4 <- ggplot(gwas_rec, aes(x=theta, y=logp, size=AC, color = effects)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-xlims, xlims) +
  ylim(0, ylims) +
  xlab("Theta") +
  ylab("-Log10(P-value)") +
  labs(color = "Bonf. Sig. Marker") +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
  ) + ggtitle("Theta (Recessive)", params)


graphics.off()
plot_grid(p1, p2, p3, p4, labels = "AUTO", ncol = 2, nrow = 2)




