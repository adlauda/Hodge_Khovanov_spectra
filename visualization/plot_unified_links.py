"""Generate outputs/plots/unified_plot_links.png."""

from pathlib import Path
import re
import sqlite3

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "plots" / "unified_plot_links.png"


def crossings(name):
    match = re.search(r"\d+", str(name))
    return int(match.group()) if match else None


def load_data():
    with sqlite3.connect(ROOT / "databases" / "link_research.db") as conn:
        minima = pd.read_sql_query(
            """
            SELECT links.name, MIN(bidegrees.smallest_nonzero) AS min_gap
            FROM links
            JOIN bidegrees ON links.id = bidegrees.link_id
            WHERE bidegrees.smallest_nonzero IS NOT NULL
              AND bidegrees.smallest_nonzero > 0
            GROUP BY links.id
            """,
            conn,
        )
        all_gaps = pd.read_sql_query(
            """
            SELECT links.name, bidegrees.smallest_nonzero AS gap
            FROM links
            JOIN bidegrees ON links.id = bidegrees.link_id
            WHERE bidegrees.smallest_nonzero IS NOT NULL
              AND bidegrees.smallest_nonzero > 0
            """,
            conn,
        )
    for frame, gap_col in [(minima, "min_gap"), (all_gaps, "gap")]:
        frame["crossings"] = frame["name"].map(crossings)
        frame.dropna(subset=["crossings", gap_col], inplace=True)
        frame["crossings"] = frame["crossings"].astype(int)
    return minima, all_gaps


def lorentzian(x, amp, ctr, wid, y0):
    return y0 + amp * wid**2 / ((x - ctr) ** 2 + wid**2)


def log_lorentzian(x, amp, ctr, wid, y0):
    return np.log(np.maximum(lorentzian(x, amp, ctr, wid, y0), 1e-15))


def r2_log(y_true, y_pred):
    mask = (y_true > 0) & (y_pred > 0)
    if mask.sum() < 2:
        return np.nan
    log_true = np.log(y_true[mask])
    log_pred = np.log(y_pred[mask])
    ss_res = np.sum((log_true - log_pred) ** 2)
    ss_tot = np.sum((log_true - np.mean(log_true)) ** 2)
    return np.nan if np.isclose(ss_tot, 0) else 1 - ss_res / ss_tot


def main():
    minima, all_gaps = load_data()
    minima = minima[minima["crossings"] <= 11]
    all_gaps = all_gaps[all_gaps["crossings"] <= 11]
    frontier = minima.loc[minima.groupby("crossings")["min_gap"].idxmin()].sort_values("crossings")
    avg_min = minima.groupby("crossings", as_index=False)["min_gap"].mean()
    avg_all = all_gaps.groupby("crossings", as_index=False)["gap"].mean()
    x_plot = np.linspace(minima["crossings"].min(), minima["crossings"].max(), 400)

    x_blue = avg_min["crossings"].to_numpy(float)
    y_blue = avg_min["min_gap"].to_numpy(float)
    p0 = (max(y_blue) - min(y_blue), np.median(x_blue), max(1.0, np.std(x_blue) / 2), max(min(y_blue), 1e-10))
    bounds = ([0, x_blue.min() - 10, 1e-8, 1e-12], [np.inf, x_blue.max() + 10, np.inf, np.inf])
    p_lor, _ = curve_fit(log_lorentzian, x_blue, np.log(y_blue), p0=p0, bounds=bounds, maxfev=50000)
    y_lor = lorentzian(x_plot, *p_lor)
    blue_r2 = r2_log(y_blue, lorentzian(x_blue, *p_lor))

    x_purple = avg_all["crossings"].to_numpy(float)
    y_purple = avg_all["gap"].to_numpy(float)
    slope_p, intercept_p = np.polyfit(x_purple, np.log(y_purple), 1)
    y_purple_fit = np.exp(slope_p * x_plot + intercept_p)
    purple_r2 = r2_log(y_purple, np.exp(slope_p * x_purple + intercept_p))

    x_green = frontier["crossings"].to_numpy(float)
    y_green = frontier["min_gap"].to_numpy(float)
    slope_g, intercept_g = np.polyfit(x_green, np.log(y_green), 1)
    y_green_fit = np.exp(slope_g * x_plot + intercept_g)
    green_r2 = r2_log(y_green, np.exp(slope_g * x_green + intercept_g))

    fig, ax = plt.subplots(figsize=(14, 9))
    ax.scatter(all_gaps["crossings"], all_gaps["gap"], alpha=0.015, color="gray", s=10, label="All encountered gaps", rasterized=True)
    ax.scatter(frontier["crossings"], frontier["min_gap"], color="green", marker="^", s=90, label="Minimal spectral gap")
    ax.scatter(avg_min["crossings"], avg_min["min_gap"], color="blue", marker="s", s=70, label="Average minimal spectral gap")
    ax.scatter(avg_all["crossings"], avg_all["gap"], color="purple", marker="o", s=70, label="Average spectral gap")
    ax.plot(x_plot, y_lor, color="blue", linestyle="--", linewidth=2.2, label=f"Lorentzian fit, R^2={blue_r2:.3f}")
    ax.plot(x_plot, y_purple_fit, color="purple", linewidth=2.2, label=f"Exponential fit, R^2={purple_r2:.3f}")
    ax.plot(x_plot, y_green_fit, color="green", linewidth=2.2, label=f"Exponential fit, R^2={green_r2:.3f}")
    ax.set_yscale("log")
    ax.set_xlabel("Crossing number")
    ax.set_ylabel("Spectral gap")
    ax.set_title("Link Khovanov Hodge Laplacian Spectral Gaps")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(loc="lower left", framealpha=0.9)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
