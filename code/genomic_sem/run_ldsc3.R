options(stringsAsFactors = FALSE)
.libPaths(c("C:/Users/HUAWEI/AppData/Local/R/win-library/4.6", .libPaths()))
suppressMessages(library(GenomicSEM))
setwd("C:/Users/HUAWEI/WorkBuddy/2026-11-05-18-21-28")

traits <- c("ldsc/ldsc-out/whradjbmi_pulit.sumstats.gz",
            "gsem/cognition.clean.sumstats.gz",
            "gsem/education.clean.sumstats.gz")
ld <- "ldsc/ldsc-ref/eur_w_ld_chr/"

cov3 <- ldsc(traits = traits,
             sample.prev = c(NA, NA, NA),
             population.prev = c(NA, NA, NA),
             ld = ld, wld = ld,
             trait.names = c("WHR", "COG", "EDU"),
             ldsc.log = "gsem/gsem3")

cat("\n=== object names ===\n")
print(names(cov3))
cat("\n=== S (genetic covariance) ===\n"); print(cov3$S)
cat("\n=== V (sampling covariance) ===\n"); print(cov3$V)
cat("\n=== N ===\n"); print(cov3$N)
cat("\n=== m ===\n"); print(cov3$m)
cat("\n=== h2 ===\n"); print(cov3$h2)
saveRDS(cov3, "gsem/covstruc3.rds")
cat("\nDONE_LDSC\n")
