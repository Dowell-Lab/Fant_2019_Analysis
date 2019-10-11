# Presentation Outline

## Introduction
- Nascent sequencing
- PolII Pausing
  - Show some images from Charli's data
- General overview of TAF1

## Method Optimization

### Differential Expression
- Show problematic early plots
- Discuss the issue with the DESeq 2 model
- Show how the weight factor adjustment solves the issue

### Metagene Analysis
- Show plot from prebuilt R package
  - Discuss default of RPM normalization and strand-unawareness
- Show iterations on custom metagene analysis, mention speedup

### PCA
- Talk about batch correction and the options available
  - Show the difference between corrected and uncorrected
- Batch correction works either by adding error to a linear model,
    or by using Empirical Bayes (prior estimation)
	- Here we use the Limma package (Speed, 2003), which uses loess normalization

### Pause Index
- List of existing methods
- Show some comparative images from Chen 2014 Paper
