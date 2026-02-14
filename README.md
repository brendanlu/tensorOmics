# tensorOmics
[![License: AGPL-3](https://img.shields.io/badge/license-AGPL--3-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![Project Status: Active](https://img.shields.io/badge/project%20status-active-brightgreen.svg)](https://example.com)

This repository contains the `tensorOmics` R package for the analysis of temporal multi-omics data using (order 3) tensor generalizations of matrix methods.

## Installation

Install the GitHub version with:

```r
install.packages("devtools")
devtools::install_github("https://github.com/brendanlu/tensorOmics")
```

## Package Website

You can view the package website [here](https://brendanlu.github.io/tensorOmics/)

## Abstract

Multi-omics studies capture comprehensive molecular profiles across biological layers to understand complex biological processes. A central challenge is integrating information across heterogeneous data types to identify coordinated molecular responses, particularly when measurements are collected longitudinally. Traditional integration methods can be broadly classified as unsupervised (exploring patterns without phenotypic information) or supervised (discriminating between groups or predicting outcomes). These approaches rely predominantly on matrix-based techniques that concatenate or project data into lower-dimensional spaces. However, matrix methods struggle with longitudinal data, as flattening multi-dimensional structures obscures temporal trajectories and violates independence assumptions. Tensor-based methods preserve the natural multi-way structure of longitudinal data but existing approaches are predominantly unsupervised, cannot incorporate phenotypic responses for discriminant analysis, and lack frameworks for integrating multiple omics layers.

We present tensorOmics, a comprehensive framework for longitudinal omics analysis using tensor factorisation. The framework encompasses supervised and unsupervised methods for both single-omic (tensor PCA, tensor PLS discriminant analysis) and multi-omic settings (tensor PLS, block tensor PLS, block tensor PLS discriminant analysis). This unified approach captures coordinated responses across biological layers while preserving temporal structure. We validated tensorOmics through three case studies: antibiotic perturbation experiments, anaerobic digestion systems, and fecal microbiota transplantation. These applications demonstrate tensorOmics differentiates treatment groups, captures time-dependent molecular signatures, and reveals multi-layer coordinated responses that cross-sectional methods miss.

## Case studies

You can view case studies [here](https://github.com/SarithaKodikara/tensorOmics_manuscript?tab=readme-ov-file)

## Link to full manuscript

You can view the preprint [here](https://www.biorxiv.org/content/10.64898/2026.02.10.705179v1)
