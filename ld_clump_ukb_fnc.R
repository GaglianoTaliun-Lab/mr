#' Clump GWAS SNPs using UKB (bgen files) as LD reference.
#' Modified function ld_clump from MRCIEU/ieugwasr (https://rdrr.io/github/MRCIEU/ieugwasr/src/R/ld_clump.R) using local reference files.
#'
#' @param dat (dataframe). Dataframe containing rsIDs (rsid) and p-values (pval) to clump.
#' @param plink_bin (file path). Path to where PLINK2 is located locally.
#' @param bgen (file path). Path to where UKB bgen files are locally.
#' @param sample (file path). Path to where UKB sample file is locally.
#' @param clump_kb (integer). Physical distance threshold in kilobases. Default = 10000
#' @param clump_r2 (decimal between 0-1). r2 LD threshold for clumping. Default = 0.001
#' @param clump_p float (decimal between 0-1). Significance threshold for SNPs. Default = 0.99
#' @export 
#' 
#' 

ld_clump_ukb <- function(dat, plink_bin, bgen, sample, clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99)
{

    # Make textfile
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	fn <- tempfile()
	write.table(data.frame(SNP=dat[["rsid"]], P=dat[["pval"]]), file=fn, row.names=F, col.names=T, quote=F)

	fun2 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bgen ", shQuote(bgen, type=shell), 
 	        " 'ref-first'",
        	" --sample ", shQuote(sample, type=shell),
		" --clump ", shQuote(fn, type=shell),
		" --clump-unphased ", 
		" --clump-p1 ", clump_p, 
		" --clump-r2 ", clump_r2, 
		" --clump-kb ", clump_kb, 
		" --out ", shQuote(fn, type=shell)
	)
	system(fun2)
	res <- read.table(paste(fn, ".clumps", sep=""), header=F)
	colnames(res) <- c("CHROM", "POS", "SNP", "pval", "TOTAL","S1", "S2", "S3", "S4", "TOTAL", "SP2")
	unlink(paste(fn, "*", sep=""))
	y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
	if(nrow(y) > 0)
	{
		message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel")
	}
	return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
}