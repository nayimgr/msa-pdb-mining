"""Shared colour choices for the visual renderers.

Centralised so the static image and the interactive HTML stay in lockstep:
the depth heatmap uses seaborn's **mako_r** colormap (light = shallow coverage,
dark = deep), and the secondary-structure cartoon uses one colour per element
(helix / strand).
"""

from __future__ import annotations

import functools

from matplotlib.colors import Colormap

# Secondary-structure cartoon colours (shared by image.py and html.py).
SS_HELIX_COLOR = "#cc3333"   # α-helix — red cylinder
SS_STRAND_COLOR = "#f5b700"  # β-strand — gold arrow


@functools.lru_cache(maxsize=1)
def depth_cmap() -> Colormap:
    """The mako colormap (from seaborn) used for structural-coverage depth.

    seaborn is imported lazily and the result memoised, so the (heavier) import
    is paid once, only when a visual renderer actually runs.
    """
    import seaborn as sns

    return sns.color_palette("mako_r", as_cmap=True)
