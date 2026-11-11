options(stringsAsFactors = FALSE)
.libPaths(c("C:/Users/HUAWEI/AppData/Local/R/win-library/4.6", .libPaths()))
suppressMessages(library(GenomicSEM))
setwd("C:/Users/HUAWEI/WorkBuddy/2026-11-05-18-21-28")

cov3 <- readRDS("gsem/covstruc3.rds")
S <- cov3$S; V <- cov3$V
lab <- c("WHR", "COG", "EDU"); h2 <- diag(S)
R <- S / sqrt(outer(h2, h2)); dimnames(R) <- list(lab, lab)

ev <- eigen((V + t(V)) / 2, symmetric = TRUE)
cat("V eigenvalues:", sprintf("%.2e", ev$values), "\n")
Vpsd <- ev$vectors %*% diag(pmax(ev$values, 1e-12)) %*% t(ev$vectors)
L <- chol(Vpsd)

r_wc <- R["WHR","COG"]; r_we <- R["WHR","EDU"]; r_ce <- R["COG","EDU"]
part <- function(a, b, c) (a - b * c) / sqrt(1 - c^2)
r_wc_e <- part(r_wc, r_we, r_ce); r_we_c <- part(r_we, r_wc, r_ce)
att <- 1 - r_wc_e / r_wc

cat(sprintf("\nrg(WHR,COG)=%.4f  rg(WHR,EDU)=%.4f  rg(COG,EDU)=%.4f\n", r_wc, r_we, r_ce))
cat(sprintf("rg(WHR,COG|EDU)=%.4f  attenuation=%.1f%%\n", r_wc_e, 100 * att))
cat(sprintf("rg(WHR,EDU|COG)=%.4f\n", r_we_c))

theta <- c(S[1,1], S[2,1], S[3,1], S[2,2], S[3,2], S[3,3])
set.seed(20261110)
B <- 20000
Z <- matrix(rnorm(B * 6), nrow = B, ncol = 6)
draws <- sweep(Z %*% L, 2, theta, "+")

res <- matrix(NA_real_, nrow = B, ncol = 5)
bad <- 0
for (b in seq_len(B)) {
  v <- draws[b, ]
  Sb <- matrix(0, 3, 3)
  Sb[1,1] <- v[1]; Sb[2,2] <- v[4]; Sb[3,3] <- v[6]
  Sb[2,1] <- Sb[1,2] <- v[2]; Sb[3,1] <- Sb[1,3] <- v[3]; Sb[3,2] <- Sb[2,3] <- v[5]
  if (!all(is.finite(Sb)) || any(diag(Sb) <= 0)) { bad <- bad + 1; next }
  Rb <- Sb / sqrt(outer(diag(Sb), diag(Sb)))
  if (any(abs(Rb[lower.tri(Rb)]) >= 1)) { bad <- bad + 1; next }
  res[b, ] <- c(Rb[1,2], Rb[1,3], Rb[2,3],
                part(Rb[1,2], Rb[1,3], Rb[2,3]),
                part(Rb[1,3], Rb[1,2], Rb[2,3]))
}
cat("discarded draws:", bad, "\n")
res <- res[stats::complete.cases(res), ]
se <- apply(res, 2, stats::sd)
ci <- apply(res, 2, stats::quantile, probs = c(0.025, 0.975))
pv <- apply(res, 2, function(x) 2 * min(mean(x > 0), mean(x < 0)))

tab <- data.frame(quantity = c("rg_WHR_COG","rg_WHR_EDU","rg_COG_EDU",
                               "rg_WHR_COG_given_EDU","rg_WHR_EDU_given_COG"),
                  est = c(r_wc, r_we, r_ce, r_wc_e, r_we_c),
                  boot_se = se, lo = ci[1,], hi = ci[2,], p = pv)
cat("\n=== bootstrap B =", nrow(res), "===\n")
print(tab, row.names = FALSE, digits = 4)
cat(sprintf("\nattenuation = %.1f%% (95%% CI %.1f%% to %.1f%%)\n",
            100 * att, 100 * (1 - ci[2,4]/r_wc), 100 * (1 - ci[1,4]/r_wc)))
utils::write.csv(tab, "gsem/decomposition_results.csv", row.names = FALSE)

# lavaan cross-check on the genetic CORRELATION matrix (parameters only)
suppressMessages(library(lavaan))
fit <- lavaan("COG ~ a*EDU\nWHR ~ c*EDU + d*COG", sample.cov = R,
              sample.nobs = 500000, std.lv = FALSE)
pe <- parameterEstimates(fit, standardized = TRUE)
cat("\n=== lavaan cross-check (estimates only) ===\n")
print(pe[pe$op == "~", c("lhs","rhs","est","std.all")], row.names = FALSE)
d_std <- pe$std.all[pe$lhs == "WHR" & pe$rhs == "COG"]
a_std <- pe$std.all[pe$lhs == "COG" & pe$rhs == "EDU"]
cat(sprintf("residual correlation from lavaan = %.4f\n", d_std * sqrt(1 - a_std^2)))
cat("\nDONE\n")
