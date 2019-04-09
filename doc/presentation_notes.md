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

### Pause Index
- List of existing methods
- Show some comparative images from Chen 2014 Paper

### PCA
- Talk about batch correction and the options available
  - Show the difference between corrected and uncorrected
  - Mention the difference between the implementations of algorithms used
