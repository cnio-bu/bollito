# Snakemake workflow: sc-test-data

[![Snakemake](https://img.shields.io/badge/snakemake-≥4.2-brightgreen.svg)](https://snakemake.bitbucket.io)
[![build status](https://gitlab.com/bu_cnio/sc-test-data/badges/master/pipeline.svg)](https://gitlab.com/bu_cnio/sc-test-data/commits/master)


This workflow creates small test datasets for single cell data analyses. The generated data is available in the folders `ref` and `reads`, such that the repository can be directly used as a git submodule for continuous integration tests.

Based on https://github.com/snakemake-workflows/ngs-test-data
