# Calculate F-statistics for each instrumental variable and keep exposures with F-statistics > 10
# Reference: https://github.com/globalbiobankmeta/multi-ancestry-pwmr/blob/main/R%20code/instrument.R

# 'exposure_data' is a formatted exposure dataframe, where beta is located in column #2 and SE is located in column #3
exposure_data <- cbind(exposure_data,fstatistics=1)

for (s in 1:nrow(exposure_data)){
    
    z <- exposure_data[s,2]/exposure_data[s,3] # calculate Z-score with beta and SE
    pve <- z^2/(sample_size_exp+z^2) # calculate the proportion of variance explained
    # fstatistic: = (N-K-1)*PVE/(1-PVE); where N = sample size; K = N SNPs (K = 1 for each instrument):
    exposure_data[s,"fstatistics"] <- (sample_size_exp-2)*pve/(1-pve)

}

# keep only instruments with F > 10: 
exposure_data <- exposure_date[exposure_data$fstatistics>10,]
