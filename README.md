# q2-bias-sim

A [QIIME 2](https://qiime2.org) plugin for applying realistic measurement biases to microbial community proportions. Useful for benchmarking source tracking and other microbiome pipelines.

## Actions

| Action | Type | Description |
|---|---|---|
| `multiplicative-bias` | Method | Apply per-taxon log-normal multiplicative bias (simulates PCR/extraction efficiency differences) |
| `detection-threshold` | Method | Zero out taxa below a detection limit and renormalize (simulates sequencing dropout) |
| `contamination` | Method | Mix a contaminant profile into samples at a specified fraction |
| `resample-counts` | Method | Draw integer counts from proportions via multinomial sampling |
| `apply-bias-pipeline` | Pipeline | Chain any combination of the above steps in one command |
| `summarize-bias` | Visualizer | Side-by-side heatmaps and scatter plot showing before/after bias effect |

## Installation

```bash
pip install -e . --no-deps
```

Install into your QIIME 2 conda environment. Verify with `qiime info`.

## Example

```bash
qiime bias-sim apply-bias-pipeline \
  --i-proportions clean_proportions.qza \
  --p-steps multiplicative threshold \
  --p-multiplicative-log-scale 1.0 \
  --p-threshold-min-abundance 0.001 \
  --p-seed 42 \
  --o-biased-proportions biased.qza
```
