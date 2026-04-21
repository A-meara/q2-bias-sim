import biom
import numpy as np
import pandas as pd

from ._core import (
    apply_multiplicative_bias,
    apply_detection_threshold,
    apply_contamination,
)


def _biom_to_df(table: biom.Table) -> pd.DataFrame:
    return pd.DataFrame(
        table.matrix_data.toarray().T,
        index=table.ids('sample'),
        columns=table.ids('observation'),
    )


def _df_to_biom(df: pd.DataFrame) -> biom.Table:
    return biom.Table(
        df.T.values.astype(float),
        observation_ids=list(df.columns),
        sample_ids=list(df.index),
    )


def multiplicative_bias(
    proportions: biom.Table,
    log_scale: float = 0.5,
    seed: int = None,
) -> biom.Table:
    """Apply per-taxon log-normal multiplicative bias."""
    df = _biom_to_df(proportions)
    biased_df, _ = apply_multiplicative_bias(df, log_scale=log_scale, seed=seed)
    return _df_to_biom(biased_df)


def detection_threshold(
    proportions: biom.Table,
    min_abundance: float = 1e-4,
    stochastic: bool = False,
    steepness: float = 2.0,
    seed: int = None,
) -> biom.Table:
    """Zero out taxa below a detection threshold, then renormalize."""
    df = _biom_to_df(proportions)
    out = apply_detection_threshold(
        df,
        min_abundance=min_abundance,
        stochastic=stochastic,
        steepness=steepness,
        seed=seed,
    )
    return _df_to_biom(out)


def contamination(
    proportions: biom.Table,
    contaminant_profile: biom.Table = None,
    fraction: float = 0.01,
    per_sample: bool = False,
    seed: int = None,
) -> biom.Table:
    """Mix in a contaminant community at a given fraction.

    If contaminant_profile is not provided, a uniform profile over all taxa
    is used.
    """
    df = _biom_to_df(proportions)

    cont_arr = None
    if contaminant_profile is not None:
        cont_df = _biom_to_df(contaminant_profile)
        # Reindex to match taxa order; missing taxa get 0
        cont_df = cont_df.reindex(columns=df.columns, fill_value=0.0)
        # Use the mean across all samples in the contaminant table
        cont_arr = cont_df.values.mean(axis=0)

    out = apply_contamination(
        df,
        contaminant=cont_arr,
        fraction=fraction,
        per_sample=per_sample,
        seed=seed,
    )
    return _df_to_biom(out)


def resample_counts(
    proportions: biom.Table,
    library_size: int = 10_000,
    seed: int = None,
) -> biom.Table:
    """Draw integer counts from proportions via multinomial sampling."""
    rng = np.random.default_rng(seed)
    df = _biom_to_df(proportions)
    counts = np.array([rng.multinomial(library_size, row) for row in df.values])
    counts_df = pd.DataFrame(counts, index=df.index, columns=df.columns)
    return biom.Table(
        counts_df.T.values.astype(int),
        observation_ids=list(counts_df.columns),
        sample_ids=list(counts_df.index),
    )
