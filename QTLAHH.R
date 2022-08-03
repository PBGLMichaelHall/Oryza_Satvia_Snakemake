#Set samples
HighBulk <- "ET_385" # SRA-run: SRR834931
LowBulk <- "ES_430" # SRA-run: SRR834927 

#set file name of the VCF file to load
file <- "freebayes~bwa~GCF_001433935.1_IRGSP-1.0~all_samples~filtered-strict~snpEff.Biallelic.vcf.gz"

#Specify which chromosomes should be included in the analysis (i.e., exclude smaller contigs)
Chroms <- c("NC_029256.1",
            "NC_029257.1",
            "NC_029258.1",
            "NC_029259.1",
            "NC_029260.1",
            "NC_029261.1",
            "NC_029262.1",
            "NC_029263.1",
            "NC_029264.1",
            "NC_029265.1",
            "NC_029266.1",
            "NC_029267.1")
df <- 
  importFromVCF(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  )


df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20,
    minTotalDepth = 70,
    maxTotalDepth = 200,
    minSampleDepth = 30,
    depthDifference = 100,
    #minGQ = 150,
    verbose = TRUE
  )

df_filt <- 
  runGprimeAnalysis(
    SNPset = df_filt,
    windowSize = 1e6,
    outlierFilter = "deltaSNP",
    filterThreshold = 0.1
  )

df_filt <- 
  runQTLseqAnalysis(
    SNPset = df_filt,
    windowSize = 1e6,
    popStruc = "F2",
    bulkSize = c(385, 430), 
    replications = 10000,
    intervals = c(95, 99)
  )

plotQTLStats(SNPset = df_filt, var = "nSNPs", plotIntervals = TRUE)
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.02)


df_filt$p1 <- df_filt$AD_ALT.LOW/df_filt$DP.LOW
df_filt$p2 <- df_filt$AD_ALT.HIGH/df_filt$DP.HIGH
df_filt3 <- df_filt %>% dplyr::filter(CHROM=="NC_029263.1")




write.csv(df_filt3,file = "QTL.csv",sep=",")

