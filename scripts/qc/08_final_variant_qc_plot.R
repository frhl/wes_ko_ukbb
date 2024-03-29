library(ggplot2)
library(ggsci)
library(data.table)

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")

suppressPackageStartupMessages(library("argparse"))

CHR <- 1
save_figures <- TRUE

# Input files
VARIANT_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/variants/08_final_qc.variants_chr', CHR, '.tsv.bgz')

# Output files
COMBINED_VARIANT_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/variants/08_final_qc.variants.tsv')

dt <- fread(cmd = paste('zcat', VARIANT_QC_FILE), header=TRUE, sep='\t')

dt_list <- list()
dt_list[[1]] <- dt

for (CHR in c(seq(2,22), "X")) {
    # Input files
    cat(paste0("chromosome ", CHR, "\n"))
    VARIANT_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/variants/08_final_qc.variants_chr', CHR, '.tsv.bgz')
    dt_tmp <- fread(cmd = paste('zcat', VARIANT_QC_FILE), header=TRUE, sep='\t')
    setkeyv(dt_tmp, c('locus', 'alleles'))
    dt_list[[CHR]] <- dt_tmp
}

dt <- rbindlist(dt_list)

fwrite(dt, file=COMBINED_VARIANT_QC_FILE , sep='\t', quote=FALSE)

dt <- fread(COMBINED_VARIANT_QC_FILE)
dt <- dt %>% filter((qc.AF > 0) & (qc.AF < 1))
# call rate across all variants
create_pretty_hist(dt, aes(x=qc.call_rate), threshold=T_variant_call_rate, x_label='Call Rate', xlim=c(0.9,1), save_figure=save_figures,
	file=paste0(PLOTS, TRANCHE, '_08_callRate_hist'))
# cumulative call rate
create_pretty_cumulative(dt, aes(x=qc.call_rate), x_label="Call Rate", threshold=T_variant_call_rate, 
    key_label='', xlim=c(0.9,1), save_figure=save_figures, file=paste0(PLOTS, TRANCHE, '_08_callRate_cdf'))

# pHWE across all variants
p <- create_pretty_hist(dt, aes(x=qc.p_value_hwe), threshold=T_pHWE,  x_label='p(HWE)',
    xlim=c(NA,1), save_figure=FALSE) + scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )
ggsave(paste0(PLOTS, TRANCHE, '_08_pHWE_hist.pdf'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, TRANCHE, '_08_pHWE_hist.jpg'), p, width=160, height=90, units='mm', dpi=500)

# cumulative pHWE
p <- create_pretty_cumulative(dt, aes(x=qc.p_value_hwe), x_label='p(HWE)',
    threshold=T_pHWE, xlim=c(NA,1), save_figure=FALSE) + scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )
ggsave(paste0(PLOTS, TRANCHE, '_08_pHWE_cdf.pdf'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, TRANCHE, '_08_pHWE_cdf.jpg'), p, width=160, height=90, units='mm', dpi=500)



