# st-metab-multi-omics-anlysis
Pipelines and major analyses for multi-omics investigation of phenolic acid accumulation in potato tubers.


This repository provides comprehensive scripts and pipelines for multi-omics data processing, systems genetics, and trait evaluation in potato (Solanum tuberosum):
	1	Heritability Estimation (h^2):
	◦	Batch calculation of narrow-sense heritability for metabolite phenotypes using GCTA (REML) based on a pre-calculated Genomic Relationship Matrix (GRM).
	2	Genome-Wide Association Studies (GWAS):
	◦	High-performance association mapping using linear mixed models (LMM) via BOLT-LMM.
	3	PEER Factor Optimization:
	◦	Latent hidden factor analysis and residual variance optimization across varying factor numbers (1 to 50) for peel and flesh samples.
	4	MatrixeQTL Mapping:
	◦	Genome-wide cis- and trans-eQTL discovery and statistical visualization using the MatrixEQTL R package.
	5	Window-Based SNP Distribution Analysis:
	◦	Custom Python utilities to segment chromosomes into fixed-size windows (e.g., 20kb) and quantify cis/trans eQTL hotspots based on correlation thresholds.
Requirements & Dependencies
	•	Python (>= 3.8) with standard scientific libraries (sys, os, collections).
	•	R (>= 3.5 / 4.x) with packages:
	◦	MatrixEQTL
	◦	snpStats
	◦	peer
	•	Bioinformatics Tools:
	◦	GCTA (v1.9x or later)
	◦	BOLT-LMM
	◦	PLINK (v1.9 / v2.0)

This repository covers the major analysis scripts for the project. Please feel free to reach out or communicate here if you need any further information.
