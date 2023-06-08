# a simple simulation workflow

require(ggplot2)
require(ggrepel)
require(cowplot)

make_beta <- function(M, h2, b, pi = 1){
  stopifnot((pi > 0 & pi <= 1))
  stopifnot((h2 >= 0 & h2 <= 1))
  the_mean <- b/(M*pi)
  the_var <-  h2/(M*pi)
  causal <- rbinom(M, size = 1, prob = pi)
  effect <- rnorm(M, mean = the_mean, sd = sqrt(the_var))
  return(causal * effect)
}

h2 <- 0.2
b <- 10
n_genes <- 1000
pi <- 0.5

# is mean correct?
X <- data.frame(x=replicate(5000, mean(make_beta(M = n_genes, h2 = h2, b = b, pi = pi))))
ggplot(X, aes(x=x)) + 
    geom_histogram(fill='grey60') + 
    geom_vline(xintercept = (b/n_genes), color='red')

# is variance correct
Y <- data.frame(x=replicate(5000, var(make_beta(M = n_genes, h2 = h2, b = b, pi = pi))))
ggplot(Y, aes(x=x)) + 
    geom_histogram(fill='grey60') + 
    geom_vline(xintercept = h2/(n_genes), color='green')



graphics.off()

# set up parameters
n_samples <- 10000 # samples
m_genes <- 1000 # genes
p_dosage <- 0.10

b <- 10
pi <- 0.5
h2 <- 0.15

# generate matrix
H1 <- matrix(rbinom(n_samples*m_genes, 1, p_dosage), nrow = m_genes, ncol = n_samples) #nolint
H2 <- matrix(rbinom(n_samples*m_genes, 1, p_dosage), nrow = m_genes, ncol = n_samples) #nolint
G <- H1+H2
rownames(G) <- paste0("chr0:",1:nrow(G),"A:A")

# normalize genotypes
G[G==1] <- 0
G[G==2] <- 1
stopifnot(sum(G)>2)

# only keep entries with SD > 0
keep <- apply(G, 1, sd) > 0
G <- G[keep, ]
new_m_genes <- nrow(G)

# standardized columns to have variance of 1 and mean of 0
G_norm <- (G - apply(G, 2, mean )) / apply(G, 2, sd)

# make effects
beta <- make_beta(new_m_genes, h2=h2, b=b, pi = pi)
df_beta <- data.frame(rsid=rownames(G), beta = beta)

# genetic effect contribution
y_no_noise <- as.vector(beta %*% G_norm)


print(mean(y_no_noise)) # should be approximately b
print(var(y_no_noise)) # should approximate h2 on avg
print(h2)

noise <-rnorm(n_samples, sd=sqrt(1-h2))
y <- y_no_noise + noise
print(var(y)) # should be 1
print(mean(y))

gwas_rec <- do.call(rbind, lapply(1:nrow(G_norm), function(i){
  fit <- summary(lm(y ~ G_norm[i,]))
  fit <- as.data.frame(t(as.matrix(fit$coefficients[2,])))
  fit$rsid <- rownames(G_norm)[i]
  colnames(fit) <- c("est","std","t","p", "rsid")
  fit <- fit[,c("rsid", "est","std","t","p")]
  fit$logp <- round(-log10(fit$p),3)
  #fit$AC <- sum(G_alt[i,])
  #fit$significant <- fit$p < (0.05/n)
  #fit$sig <- ifelse(fit$significant,'*','')
  return(fit)
}))

gwas_rec

gwas_rec <- merge(gwas_rec, df_beta, all.x = TRUE, by = 'rsid')

head(gwas_rec)

ggplot(gwas_rec, aes(x=beta, y=-log10(p))) + 
  geom_point() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  xlab("Beta") + 
  ylab("-log10(P)")


