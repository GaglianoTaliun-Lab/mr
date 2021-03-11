##DESCRIPTION
#Example script to plot a Locus Zoom plot for a region on chr3:186059360-187058330 (e.g. --chr 3 --start 186059360 --end 187058330)
#Documentation: https://genome.sph.umich.edu/wiki/LocusZoom_Standalone#Download

#This example script  uses 1000 Genomes European ancestry data to calculate the correlation (LD) between nearby SNPs in the region
#This example is for human genome build 37 (aka --build hg19; if you are using build 38, use --build 38)
#The file MySummaryStats-4LZ.txt has a whitespace delimiter, and column called MARKER is chr:pos (OR it could be rsID) and the p-value columns is called PVALUE
#E.g. --markercol MARKER --delim whitespace --pvalcol PVALUE

locuszoom --build hg19 --metal MySummaryStats-4Manhattan.txt --pop EUR --source 1000G_Nov2014 --chr 3 --start 186059360 --end 187058330 --plotonly --markercol MARKER --delim whitespace --pvalcol PVALUE --source 1000G_Nov2014 --build hg19 --pop EUR
