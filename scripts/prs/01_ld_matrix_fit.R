

library(bigsnpr)
library(data.table)
library(argparse)
#options(bigstatsr.check.parallel.blas = FALSE)
#options(default.nproc.blas = NULL)

# load helpers (some functions required an active internet connection, 
# these have been re-written to rely on locally available data)
source('utils/modules/R/bigsnpr_helpers.R')


# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_bed", default=NULL, help = "Path for plink file (bed)")
parser$add_argument("--path_snp_match", default=NULL, help = "Path to a summary statistics file with matching SNPs")
parser$add_argument("--out_prefix", default=NULL, help = "Where should the results be written?")
args <- parser$parse_args()

# check input
stopifnot(file.exists(args$path_bed))
stopifnot(file.exists(args$path_snp_match))

## copied from https://choishingwan.github.io/PRS-Tutorial/ldpred/
# Get maximum amount of cores
NCORES <- nb_cores()
tmp <- tempfile(tmpdir = "data/tmp/tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL

# We want to know the ordering of samples in the bed file 
fam.order <- NULL

# preprocess the bed file (only need to do once for each data set) and attach 
#snp_readBed(path_bed)
basename <- tools::file_path_sans_ext(path_bed)
rds <- paste0(basename,'.rds')
obj.bigSNP <- snp_attach(rds)

# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")


# read summary statistics file for SNP matching
sumstats <- bigreadr::fread2(path_sumstats)

# make artificial rsid
sumstats$rsid <- 
    apply(
        data.frame(
            locus = sumstats$locus, 
            a0 = sumstats$a0,
            a1 = sumstats$a1
            ),
        1,
        paste, 
        collapse = '_')

# extract position and chromosome
locus <- do.call(rbind, strsplit(sumstats$locus, split = "\\:"))
colnames(locus) <- c('chr','pos')
sumstats <- cbind(sumstats, locus)

# calculate MAF
#sumstats$AF <- sumstats$AC / sumstats$AN
sumstats$MAF <- unlist(lapply(sumstats$AF, function(af) min(af, 1-af)))

sumstats <- 
    data.frame(
        chr = sumstats$chr,
        pos = as.integer(sumstats$pos),
        rsid = sumstats$rsid,
        a1 = sumstats$a1,
        a0 = sumstats$a0,
        n_eff = sumstats$n,
        beta_se = sumstats$standard_error,
        p = sumstats$p_value,
        beta = sumstats$beta,
        INFO = NA,
        MAF = sumstats$MAF
    )



# perform SNP matching
info_snp <- bigsnpr::snp_match(sumstats, map, strand_flip = FALSE)
#info_snp <- snp_match(sumstats, map)

# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes

# Rename the data structures
CHR <- as.numeric(gsub('chr','',map$chr))
POS <- map$pos

# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
#snp_asGeneticPos(CHR, POS, dir = ".")
POS2 <- snp_asGeneticPosLocal(CHR, POS, mapdir = "data/prs/1000-genomes-genetic-maps",genetic_map = 'hapmap')

# calculate LD
chrs <- paste0("chr",1:22)
for (chr in chrs) {
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    stopifnot(length(ind.chr) > 0)
    # Calculate the LD
    corr0 <- snp_cor(
            genotype,
            ind.col = ind.chr2,
            ncores = NCORES,
            infos.pos = POS2[ind.chr2],
            size = 3 / 1000
        )
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}
# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
        c("family.ID", "sample.ID"),
        c("FID", "IID"))





