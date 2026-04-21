"""Pipeline action: chains individual bias methods via QIIME2 ctx."""

_VALID_STEPS = {'multiplicative', 'threshold', 'contamination'}


def apply_bias_pipeline(
    ctx,
    proportions,
    steps,
    multiplicative_log_scale=0.5,
    multiplicative_seed=None,
    threshold_min_abundance=1e-4,
    threshold_stochastic=False,
    threshold_steepness=2.0,
    threshold_seed=None,
    contamination_fraction=0.01,
    contamination_per_sample=False,
    contamination_seed=None,
    seed=None,
):
    """Apply an ordered sequence of bias steps.

    Steps are applied in the order given by the ``steps`` parameter.
    Each step's provenance is captured independently by QIIME2.
    """
    invalid = set(steps) - _VALID_STEPS
    if invalid:
        raise ValueError(
            f"Unknown step(s): {invalid}. Valid steps: {_VALID_STEPS}"
        )

    mult = ctx.get_action('bias_sim', 'multiplicative_bias')
    thresh = ctx.get_action('bias_sim', 'detection_threshold')
    cont = ctx.get_action('bias_sim', 'contamination')

    current = proportions
    for i, step in enumerate(steps):
        # base seed overrides per-step seeds with deterministic offsets
        step_seed = (seed + i + 1) if seed is not None else None

        if step == 'multiplicative':
            s = step_seed if seed is not None else multiplicative_seed
            result = mult(
                proportions=current,
                log_scale=multiplicative_log_scale,
                seed=s,
            )
            current = result.biased_proportions

        elif step == 'threshold':
            s = step_seed if seed is not None else threshold_seed
            result = thresh(
                proportions=current,
                min_abundance=threshold_min_abundance,
                stochastic=threshold_stochastic,
                steepness=threshold_steepness,
                seed=s,
            )
            current = result.thresholded_proportions

        elif step == 'contamination':
            s = step_seed if seed is not None else contamination_seed
            result = cont(
                proportions=current,
                fraction=contamination_fraction,
                per_sample=contamination_per_sample,
                seed=s,
            )
            current = result.contaminated_proportions

    return (current,)
