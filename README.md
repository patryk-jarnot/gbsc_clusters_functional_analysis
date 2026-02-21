# Clusters Functional Analysis

A tools for functional analysis of protein clusters.

## Description

This project contains Python scripts for:
- Downloading Gene Ontology (GO) annotations for proteins from GBSC clusters
- Performing statistical functional analysis of clusters
- Generating analysis reports

## Requirements

- Python 3.7+
- Libraries: `requests`, `pandas`, `numpy`, `scipy`, `sqlite3`, `statsmodels`


## Installation

Create and activate a virtual environment:
```bash
python3 -m venv venv
source venv/bin/activate
```

Install the package:
```bash
pip3 install git+https://github.com/patryk-jarnot/gbsc_clusters_functional_analysis.git
```

## Usage

Run the analysis:

```bash
python3 -m cfanalysis.cfpipeline -g [path_to_clusters]
```

Example:

```bash
python3 -m cfanalysis.cfpipeline -g input/clusters_acc
```

