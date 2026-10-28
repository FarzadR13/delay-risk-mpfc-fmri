"""Plot helpers shared by the non-brain figure panels."""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from common import C

BLUE, GREEN = "#4C72B0", "#55A868"


def ols_band(ax, x, y, color="black", band="#D9D9D9", n_grid=100, line_width=2):
    """Ordinary least-squares line with the 95% confidence band of the fitted mean."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    b, a = np.polyfit(x, y, 1)
    g = np.linspace(x.min(), x.max(), n_grid)
    s = np.sqrt(np.sum((y - (a + b * x)) ** 2) / (n - 2))
    se = s * np.sqrt(1 / n + (g - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2))
    tc = stats.t.ppf(0.975, n - 2)
    ax.fill_between(g, a + b * g - tc * se, a + b * g + tc * se, color=band, alpha=0.9, lw=0, zorder=1)
    ax.plot(g, a + b * g, color=color, lw=line_width, zorder=2)


def p_text(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


def clean_axes(ax):
    ax.spines[["top", "right"]].set_visible(False)


def save_figure(fig, name):
    C.FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(C.FIGURES / f"{name}.{ext}", dpi=300, bbox_inches="tight", pad_inches=0.1, facecolor=fig.get_facecolor())
    plt.close(fig)
