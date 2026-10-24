```
1. Heritability Estimation using GCTABefore performing association mapping, the narrow-sense heritability ($h^2$) of each metabolite phenotype was estimated using GCTA based on the genomic relationship matrix (GRM).
#!/bin/bash
# Batch calculation of metabolite heritability (h²) using GCTA

# =============== Configuration ===============
PLINK_BED="geno_snps.qc"         # SNP data prefix (.bed/.bim/.fam required)
GRM_PREFIX="my_grm"             # Pre-calculated GRM prefix (e.g., my_grm.grm.bin)
PHENO_DIR="pheno1"              # Phenotype files directory
OUT_DIR="pheno1_snp_heritability_results" # Output directory
THREADS=4                       # Number of threads
# =============================================

# Create output directory
mkdir -p "$OUT_DIR"

# Check if GRM exists
if [ ! -f "${GRM_PREFIX}.grm.bin" ]; then
    echo "Error: GRM file ${GRM_PREFIX}.grm.bin does not exist!"
    echo "Please run first: gcta64 --bfile $PLINK_BED --make-grm --out $GRM_PREFIX"
    exit 1
fi

# Loop through all phenotype files
for pheno_file in "$PHENO_DIR"/*.txt; do
    [ -e "$pheno_file" ] || continue
    
    # Extract metabolite name (remove path and .txt extension)
    metab=$(basename "$pheno_file" .txt)
    
    echo "Processing: $metab"
    
    # Run GCTA REML
    gcta64 \
        --reml \
        --grm "$GRM_PREFIX" \
        --pheno "$pheno_file" \
        --out "$OUT_DIR/$metab" \
        --thread-num "$THREADS"
    
    # Check execution status
    if [ $? -ne 0 ]; then
        echo "Warning: Calculation failed for $metab!"
    fi
done

echo "All heritability estimates completed! Results saved in $OUT_DIR/"
```

```
2. GWAS Analysis using BOLT-LMM

Association analysis performed using linear mixed models (LMM) via BOLT-LMM:
bolt --bfile=geno_snps.qc \
     --phenoFile=pheno1/M163T376.txt \
     --modelSnps=geno_snps_pruned.prune.in \
     --maxModelSnps=6000000 \
     --phenoCol=PHENO \
     --lmmInfOnly \
     --numThreads=10 \
     --statsFile=results1_snp/M163T376

```

eqtl analysis

```
PEER Analysis and MatrixeQTL Pipelines

This section contains the R scripts used for PEER factor optimization (evaluating residual variance across different latent factor numbers) and cis/trans eQTL mapping using MatrixeQTL.
1. PEER Factor Optimization Scripts

These scripts loop through different numbers of hidden factors (from 1 to 50) using the PEER framework to evaluate residual variance patterns for peel and flesh samples respectively.
Script A: Peel PEER Optimization (peel_plot_r.r)
#!/mnt/data/qinggang.liao/02_software/01_miniconda/envs/r351/bin/Rscript
library(peer)
load("paper-section-6-4.RData") #Expression matrix after filtering according to the paper description

p_results <- data.frame(
  factor_number = seq(1:50),
  total_variance = numeric(length(seq(1:50))),
  mean_gene_variance = numeric(length(seq(1:50)))
)

for (k in seq(1:50)){
  model = PEER()
  PEER_setPhenoMean(model, expr_log_norm_p)
  PEER_setNk(model, k)
  PEER_setNmax_iterations(model, 10)
  PEER_update(model)
  residuals <- PEER_getResiduals(model)

  p_results$total_variance[k] <- var(as.vector(residuals), na.rm = TRUE)
  p_results$mean_gene_variance[k] <- mean(apply(residuals, 2, var, na.rm = TRUE))
}

save(p_results, file = "peel_peer_plot.RData")
```

```
2. MatrixeQTL Analysis Pipeline (matrix_eqtl_snp_r.r)

This script performs genome-wide cis- and trans-eQTL mapping utilizing the MatrixEQTL R package with custom thresholds and positional annotation data.
#!/mnt/data/qinggang.liao/02_software/01_R-sts/bin/Rscript
library(MatrixEQTL)

# ==================== 1. Configuration ====================
base.dir <- "/pathway/to/06_matrixeqtl"
useModel <- modelLINEAR  # Model type: modelLINEAR, modelANOVA, modelLINEAR_CROSS
pvOutputThreshold_cis <- 1e-2    # p-value threshold for cis-eQTLs
pvOutputThreshold_trans <- 1e-5  # p-value threshold for trans-eQTLs
cisDist <- 33000                 # cis-window distance threshold (bp)
output_file_name_cis <- tempfile()  
output_file_name_trans <- tempfile() 
errorCovariance <- numeric()     # Error covariance matrix (default identity matrix)

# Input file paths
SNP_file_name <- file.path(base.dir, "snp_r.txt")
expression_file_name <- file.path(base.dir, "r_expr_snp.txt")
covariates_file_name <- file.path(base.dir, "cov_r_t.txt")
snpspos_file_name <- file.path(base.dir, "snp.loc.txt")
genepos_file_name <- file.path(base.dir, "r_gene_loc.txt")

# ==================== 2. Load Data ====================
# Load Genotype Data
snps <- SlicedData$new()
snps$fileDelimiter <- "\t"
snps$fileOmitCharacters <- "NA"
snps$fileSkipRows <- 1
snps$fileSkipColumns <- 1
snps$fileSliceSize <- 2000
snps$LoadFile(SNP_file_name)

# Load Expression Data
gene <- SlicedData$new()
gene$fileDelimiter <- "\t"
gene$fileOmitCharacters <- "NA"
gene$fileSkipRows <- 1
gene$fileSkipColumns <- 1
gene$fileSliceSize <- 2000
gene$LoadFile(expression_file_name)

# Load Covariates (if available)
cvrt <- SlicedData$new()
cvrt$fileDelimiter <- "\t"
cvrt$fileOmitCharacters <- "NA"
cvrt$fileSkipRows <- 1
cvrt$fileSkipColumns <- 1
if (file.exists(covariates_file_name)) {
  cvrt$LoadFile(covariates_file_name)
}

# Load SNP and Gene Positions
snpspos <- read.table(snpspos_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos <- read.table(genepos_file_name, header = TRUE, stringsAsFactors = FALSE)

# ==================== 3. Run Matrix eQTL ====================
me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_trans,
  pvOutputThreshold = pvOutputThreshold_trans,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

# Clean up temporary files
unlink(output_file_name_cis)
unlink(output_file_name_trans)

# ==================== 4. Outputs & Results ====================
cat(sprintf("Analysis completed in: %.2f seconds\n", me$time.in.sec))
cat("Detected cis-eQTLs:\n")
print(head(me$cis$eqtl, 10))
cat("Detected trans-eQTLs:\n")
print(head(me$trans$eqtl, 10))

# Save results to output files
write.table(me$cis$eqtl, sep = '\t', quote = FALSE, file.path(base.dir, "cis_eQTL_snp_r_results.csv"), row.names = FALSE)
write.table(me$trans$eqtl, sep = '\t', quote = FALSE, file.path(base.dir, "trans_eQTL_snp_r_results.csv"), row.names = FALSE)

# ==================== 5. Visualization & Statistics ====================
png(file.path(base.dir, "p_snp_r_value_qqplot.png"))
plot(me)
dev.off()

cat(sprintf("Number of significant cis-eQTLs (p < %.2e): %d\n", pvOutputThreshold_cis, nrow(me$cis$eqtl)))
cat(sprintf("Number of significant trans-eQTLs (p < %.2e): %d\n", pvOutputThreshold_trans, nrow(me$trans$eqtl)))
```
