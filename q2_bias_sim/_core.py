"""
Bias/perturbation functions — copied from mock_communities.py (Part 2).
Edit mock_communities.py (for notebooks) and this file independently.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def apply_multiplicative_bias(
    proportions: pd.DataFrame,
    log_scale: float = 0.5,
    bias_factors: np.ndarray | None = None,
    seed: int | None = None,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Apply per-taxon log-normal multiplicative bias, then renormalize."""
    rng = np.random.default_rng(seed)
    n_taxa = proportions.shape[1]

    if bias_factors is None:
        bias_factors = rng.lognormal(mean=0.0, sigma=log_scale, size=n_taxa)

    biased = proportions.values * bias_factors[np.newaxis, :]
    row_sums = biased.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, 1.0, row_sums)
    biased = biased / row_sums

    return (
        pd.DataFrame(biased, index=proportions.index, columns=proportions.columns),
        bias_factors,
    )


def apply_detection_threshold(
    proportions: pd.DataFrame,
    min_abundance: float = 1e-4,
    stochastic: bool = False,
    steepness: float = 2.0,
    seed: int | None = None,
) -> pd.DataFrame:
    """Zero out taxa below a detection threshold, then renormalize."""
    rng = np.random.default_rng(seed)
    vals = proportions.values.copy()

    if stochastic:
        log_p = np.log(vals + 1e-30)
        log_thresh = np.log(min_abundance)
        detect_prob = 1.0 / (1.0 + np.exp(-steepness * (log_p - log_thresh)))
        mask = rng.random(vals.shape) < detect_prob
        vals = vals * mask
    else:
        vals[vals < min_abundance] = 0.0

    row_sums = vals.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, 1.0, row_sums)
    vals = vals / row_sums

    return pd.DataFrame(vals, index=proportions.index, columns=proportions.columns)


def apply_contamination(
    proportions: pd.DataFrame,
    contaminant: np.ndarray | pd.Series | None = None,
    fraction: float = 0.01,
    per_sample: bool = False,
    seed: int | None = None,
) -> pd.DataFrame:
    """Mix in a contaminant community at a given fraction."""
    rng = np.random.default_rng(seed)
    n_samples, n_taxa = proportions.shape

    if contaminant is None:
        contaminant = np.ones(n_taxa) / n_taxa
    else:
        contaminant = np.asarray(contaminant, dtype=float)
        if contaminant.shape[0] != n_taxa:
            raise ValueError(
                f"contaminant length ({contaminant.shape[0]}) != n_taxa ({n_taxa})"
            )
        s = contaminant.sum()
        if s > 0:
            contaminant = contaminant / s

    if per_sample:
        a = 2.0
        b = a * (1.0 / max(fraction, 1e-8) - 1.0)
        fractions = rng.beta(a, b, size=n_samples)[:, np.newaxis]
    else:
        fractions = fraction

    mixed = (1.0 - fractions) * proportions.values + fractions * contaminant[np.newaxis, :]
    return pd.DataFrame(mixed, index=proportions.index, columns=proportions.columns)


def plot_bias_effect(
    before: pd.DataFrame,
    after: pd.DataFrame,
    title: str | None = None,
    figsize: tuple = (12, 5),
    log_scale: bool = True,
):
    """Side-by-side heatmaps showing bias effect on community proportions."""
    fig, axes = plt.subplots(1, 3, figsize=figsize,
                              gridspec_kw={"width_ratios": [1, 1, 0.8]})
    if title:
        fig.suptitle(title)

    for ax, data, panel_title in zip(
        axes[:2], [before, after], ["Before bias", "After bias"]
    ):
        vals = data.values.copy()
        if log_scale:
            vals = np.log10(vals + 1e-6)
        im = ax.imshow(vals, aspect="auto", cmap="viridis", interpolation="none")
        ax.set_title(panel_title)
        ax.set_xlabel("Taxa")
        ax.set_ylabel("Samples")
        cb = fig.colorbar(im, ax=ax, fraction=0.04)
        cb.set_label("log₁₀(proportion + 1e-6)" if log_scale else "proportion")

    ax_sc = axes[2]
    mean_before = np.maximum(before.mean(axis=0).values, 1e-6)
    mean_after = np.maximum(after.mean(axis=0).values, 1e-6)
    ax_sc.scatter(mean_before, mean_after, s=15, alpha=0.6, edgecolors="k", linewidths=0.3)
    lims = [1e-6, max(mean_before.max(), mean_after.max()) * 1.2]
    ax_sc.plot(lims, lims, "k--", alpha=0.3, linewidth=0.8)
    ax_sc.set_xscale("log")
    ax_sc.set_yscale("log")
    ax_sc.set_xlabel("Mean proportion (before)")
    ax_sc.set_ylabel("Mean proportion (after)")
    ax_sc.set_title("Per-taxon shift")

    rect = [0, 0, 1, 0.88] if title else [0, 0, 1, 1]
    fig.tight_layout(rect=rect)
    return fig, axes
