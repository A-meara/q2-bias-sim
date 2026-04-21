import os

import biom
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ._core import plot_bias_effect


def _biom_to_df(table: biom.Table) -> pd.DataFrame:
    return pd.DataFrame(
        table.matrix_data.toarray().T,
        index=table.ids('sample'),
        columns=table.ids('observation'),
    )


def summarize_bias(output_dir: str, before: biom.Table, after: biom.Table,
                   log_scale: bool = True, title: str = None) -> None:
    """Side-by-side heatmaps and per-taxon scatter showing bias effect."""
    before_df = _biom_to_df(before)
    after_df = _biom_to_df(after)

    fig, _ = plot_bias_effect(
        before_df, after_df, title=title, log_scale=log_scale
    )
    fig.savefig(
        os.path.join(output_dir, 'bias_effect.png'), dpi=150, bbox_inches='tight'
    )
    plt.close(fig)

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write(
            '<!DOCTYPE html><html><head><title>Bias Effect Summary</title>'
            '<style>body{font-family:sans-serif;padding:20px}'
            'img{max-width:100%;height:auto}</style></head>'
            '<body><h2>Bias Effect Summary</h2>'
            '<img src="bias_effect.png" alt="Bias effect plot"/>'
            '</body></html>'
        )
