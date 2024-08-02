

ldsc <- fread("data/prs/validation/ldsc_summary.txt.gz")
pred <- fread("data/prs/validation/230331_test_spiro_pgs_auc_summary.txt.gz")

setkeyv(ldsc, "phenotype")
setkeyv(pred, "phenotype")
ldsc_intercept <- ldsc[ldsc$coef == "intercept",]
ldsc_slope <- ldsc[ldsc$coef == "h2",]

mrg_slope <- merge(ldsc_slope, pred, all.x = TRUE)
mrg_intercept <- merge(ldsc_intercept, pred, all.x=TRUE)
mrg_intercept$auc_mean <- NA
mrg_intercept$auc_2_5_pct <- NA
mrg_intercept$auc_97_5_pct <- NA
mrg_intercept$auc_sd <- NA
mrg_intercept$pred_cases <- NA
mrg_intercept$pred_controls <- NA
mrg_intercept$pred_n <- NA
final <- rbind(mrg_slope, mrg_intercept)
final <- final[order(final$phenotype, final$coef),]
fwrite(final, "data/prs/validation/230418_ldsc_auc_table.txt", sep="\t")
