from qiime2.plugin import (
    Plugin, Int, Float, Bool, Str, List, Range,
)
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency

from . import __version__
from ._methods import multiplicative_bias, detection_threshold, contamination, resample_counts
from ._pipelines import apply_bias_pipeline
from ._visualizers import summarize_bias

try:
    from qiime2.plugin import Optional as Q2Optional
    _OPTIONAL_TABLE = Q2Optional[FeatureTable[RelativeFrequency]]
except ImportError:
    # Older QIIME2 versions: fall back to required input
    _OPTIONAL_TABLE = FeatureTable[RelativeFrequency]

plugin = Plugin(
    name='bias-sim',
    version=__version__,
    website='https://github.com/example/q2-bias-sim',
    package='q2_bias_sim',
    description=(
        'Apply realistic measurement biases (multiplicative, detection threshold, '
        'contamination) to microbial community proportions.'
    ),
    short_description='Microbial community bias simulator.',
)

plugin.methods.register_function(
    function=multiplicative_bias,
    inputs={
        'proportions': FeatureTable[RelativeFrequency],
    },
    parameters={
        'log_scale': Float % Range(0, None),
        'seed': Int,
    },
    outputs=[
        ('biased_proportions', FeatureTable[RelativeFrequency]),
    ],
    input_descriptions={
        'proportions': 'Input relative abundance table.',
    },
    parameter_descriptions={
        'log_scale': (
            'Std dev of the log-normal bias factor distribution. '
            '0 = no bias, larger = more distortion.'
        ),
        'seed': 'Random seed.',
    },
    output_descriptions={
        'biased_proportions': 'Proportions after multiplicative bias and renormalization.',
    },
    name='Apply multiplicative bias',
    description=(
        'Apply per-taxon log-normal multiplicative bias to simulate PCR/extraction '
        'efficiency differences, then renormalize rows to the simplex.'
    ),
)

plugin.methods.register_function(
    function=detection_threshold,
    inputs={
        'proportions': FeatureTable[RelativeFrequency],
    },
    parameters={
        'min_abundance': Float % Range(0, 1, inclusive_end=True),
        'stochastic': Bool,
        'steepness': Float % Range(0, None),
        'seed': Int,
    },
    outputs=[
        ('thresholded_proportions', FeatureTable[RelativeFrequency]),
    ],
    input_descriptions={
        'proportions': 'Input relative abundance table.',
    },
    parameter_descriptions={
        'min_abundance': 'Taxa with abundance below this value are zeroed out.',
        'stochastic': (
            'If true, use a sigmoid detection probability instead of a hard cutoff.'
        ),
        'steepness': 'Sigmoid steepness when stochastic=True. Higher = sharper cutoff.',
        'seed': 'Random seed.',
    },
    output_descriptions={
        'thresholded_proportions': 'Proportions with sub-threshold taxa removed and renormalized.',
    },
    name='Apply detection threshold',
    description=(
        'Zero out taxa below a detection threshold, then renormalize. '
        'Simulates sequencing dropout of rare taxa.'
    ),
)

plugin.methods.register_function(
    function=contamination,
    inputs={
        'proportions': FeatureTable[RelativeFrequency],
        'contaminant_profile': _OPTIONAL_TABLE,
    },
    parameters={
        'fraction': Float % Range(0, 1, inclusive_end=True),
        'per_sample': Bool,
        'seed': Int,
    },
    outputs=[
        ('contaminated_proportions', FeatureTable[RelativeFrequency]),
    ],
    input_descriptions={
        'proportions': 'Input relative abundance table.',
        'contaminant_profile': (
            'A FeatureTable[RelativeFrequency] whose mean row is used as the '
            'contaminant profile. If not provided, a uniform profile over all '
            'taxa is used.'
        ),
    },
    parameter_descriptions={
        'fraction': (
            'Contamination fraction. '
            'Result = (1 - fraction) * original + fraction * contaminant.'
        ),
        'per_sample': (
            'If true, fraction varies per sample from a Beta distribution '
            'with mean equal to fraction.'
        ),
        'seed': 'Random seed.',
    },
    output_descriptions={
        'contaminated_proportions': 'Proportions after contamination mixing.',
    },
    name='Apply contamination',
    description=(
        'Mix a contaminant community profile into sample proportions at a '
        'specified fraction. Simulates reagent or environmental contamination.'
    ),
)

plugin.methods.register_function(
    function=resample_counts,
    inputs={
        'proportions': FeatureTable[RelativeFrequency],
    },
    parameters={
        'library_size': Int % Range(1, None),
        'seed': Int,
    },
    outputs=[
        ('counts', FeatureTable[Frequency]),
    ],
    input_descriptions={
        'proportions': 'Relative abundance table to resample from.',
    },
    parameter_descriptions={
        'library_size': 'Sequencing depth per sample (multinomial draw depth).',
        'seed': 'Random seed.',
    },
    output_descriptions={
        'counts': 'Integer count table from multinomial resampling.',
    },
    name='Resample counts from proportions',
    description=(
        'Draw integer sequencing counts from a proportions table via '
        'multinomial sampling. Use this to convert biased proportions back '
        'to a count table for downstream QIIME 2 analyses.'
    ),
)

plugin.pipelines.register_function(
    function=apply_bias_pipeline,
    inputs={
        'proportions': FeatureTable[RelativeFrequency],
    },
    parameters={
        'steps': List[Str],
        'multiplicative_log_scale': Float % Range(0, None),
        'multiplicative_seed': Int,
        'threshold_min_abundance': Float % Range(0, 1, inclusive_end=True),
        'threshold_stochastic': Bool,
        'threshold_steepness': Float % Range(0, None),
        'threshold_seed': Int,
        'contamination_fraction': Float % Range(0, 1, inclusive_end=True),
        'contamination_per_sample': Bool,
        'contamination_seed': Int,
        'seed': Int,
    },
    outputs=[
        ('biased_proportions', FeatureTable[RelativeFrequency]),
    ],
    input_descriptions={
        'proportions': 'Clean input proportions.',
    },
    parameter_descriptions={
        'steps': (
            "Ordered list of bias steps to apply. "
            "Valid values: 'multiplicative', 'threshold', 'contamination'. "
            "Example: --p-steps multiplicative threshold"
        ),
        'multiplicative_log_scale': 'log_scale for the multiplicative step.',
        'multiplicative_seed': 'Seed for the multiplicative step (overridden by seed).',
        'threshold_min_abundance': 'min_abundance for the threshold step.',
        'threshold_stochastic': 'stochastic flag for the threshold step.',
        'threshold_steepness': 'steepness for the threshold step.',
        'threshold_seed': 'Seed for the threshold step (overridden by seed).',
        'contamination_fraction': 'fraction for the contamination step.',
        'contamination_per_sample': 'per_sample flag for the contamination step.',
        'contamination_seed': 'Seed for the contamination step (overridden by seed).',
        'seed': (
            'Base random seed. When set, overrides all per-step seeds with '
            'deterministic offsets (seed+1, seed+2, …).'
        ),
    },
    output_descriptions={
        'biased_proportions': 'Final proportions after all bias steps.',
    },
    name='Apply bias pipeline',
    description=(
        'Apply a composable, ordered sequence of bias steps '
        '(multiplicative, threshold, contamination) to a proportions table. '
        'Each step is recorded individually in QIIME 2 provenance.'
    ),
)

plugin.visualizers.register_function(
    function=summarize_bias,
    inputs={
        'before': FeatureTable[RelativeFrequency],
        'after': FeatureTable[RelativeFrequency],
    },
    parameters={
        'log_scale': Bool,
        'title': Str,
    },
    input_descriptions={
        'before': 'Clean proportions before bias.',
        'after': 'Biased proportions after pipeline.',
    },
    parameter_descriptions={
        'log_scale': 'If true, display log₁₀-scaled values in heatmaps.',
        'title': 'Optional figure title.',
    },
    name='Summarize bias effect',
    description=(
        'Generate side-by-side heatmaps (before vs. after) and a per-taxon '
        'mean-shift scatter plot showing the effect of bias on proportions.'
    ),
)
