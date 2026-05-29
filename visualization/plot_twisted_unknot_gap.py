"""Generate outputs/plots/twisted_unknot_gap.png."""

from pathlib import Path
import sqlite3

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "plots" / "twisted_unknot_gap.png"


def main():
    with sqlite3.connect(ROOT / "databases" / "twisted_unknot_research.db") as conn:
        data = pd.read_sql_query(
            """
            SELECT twisted_unknots.crossings,
                   MIN(bidegrees.smallest_nonzero) AS min_gap
            FROM twisted_unknots
            JOIN bidegrees ON twisted_unknots.id = bidegrees.twisted_unknot_id
            WHERE bidegrees.smallest_nonzero IS NOT NULL
              AND bidegrees.smallest_nonzero > 0
            GROUP BY twisted_unknots.id, twisted_unknots.crossings
            ORDER BY twisted_unknots.crossings
            """,
            conn,
        )
    if data.empty:
        raise SystemExit("no twisted-unknot spectral gaps found")

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(data["crossings"], data["min_gap"], marker="o", linewidth=2.2)
    ax.set_yscale("log")
    ax.set_xlabel("n")
    ax.set_ylabel("Minimal spectral gap")
    ax.set_title("Twisted Unknot Spectral Gap")
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
