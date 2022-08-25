# RNAseq

![Status](https://img.shields.io/badge/status-alpha-red)

Latest Release:
* Github:![Github Release](https://img.shields.io/badge/release-v1-blue)

Code for analyzing publicly available RNAseq data using Salmon quasi mapper algorithm.

## SALMON-v.2.1.0
### 2.1.0

	* New features:

	1) Speed/memory optimization.  Faster.


## Install

The common way to install Salmon is through
[PYPI](https://pypi.org/project/salmon/)
* x86_64

## Usage

Example for regular peak calling on TF ChIP-seq:

`macs3 callpeak -t ChIP.bam -c Control.bam -f BAM -g hs -n test -B -q 0.01`


Subcommand | Description
-----------|----------
[`callpeak`](./docs/callpeak.md) | Main MACS3 Function to call peaks from alignment results.
[`bdgpeakcall`](./docs/bdgpeakcall.md) | Call peaks from bedGraph output.
