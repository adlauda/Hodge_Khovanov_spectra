"""Reproduce conservative numerical figures without rerunning an eigensolver.

Run from any directory. Input SQLite databases are opened read-only. The knot
database is a checked local derivative made by clean_knot_database.py; the other
databases are unchanged snapshots from SOURCE_COMMIT. The manifest records
source and working-copy SHA-256 hashes. New PNGs and derived tables are written under
this submission directory. 'Complete' means the source database contains a
non-NULL khovanov_polynomial, not an independent proof of numerical correctness.
"""
from pathlib import Path
import datetime
import hashlib
import json
import os
import re
import sqlite3

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parents[2]
REPO_DATABASES = REPO_ROOT / "databases"
SNAPSHOTS = ROOT / "snapshots"
if not SNAPSHOTS.exists():
    SNAPSHOTS = ROOT.parents[1] / "Backup/2026-09-08-R4/reproducibility/snapshots"
CLEANED_KNOTS = ROOT / "data/knot_research.db"
if not CLEANED_KNOTS.is_file():
    CLEANED_KNOTS = REPO_DATABASES / "knot_research.db"
FIGURES = ROOT.parent / "figures"
SOURCE_REPO = "https://github.com/adlauda/Hodge_Khovanov_spectra"
SOURCE_COMMIT = "ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea"
SOURCE_SHA256 = {
    "knot_research.db": "918832fc674a8b92aeef88f464bad18c4a21478ef90fefa2e57f52aed1d1e21a",
    "knot_research_14.db": "2c3e9d55e764e6bc10a01abacdc00b850ec281e68536ba83f6c5302656c5a3c2",
    "link_research.db": "6c0f6ee56d9fdc327c7076cc1fd637e6f5758c0961e4a17db776a16a0fb1919c",
    "twisted_unknot_research.db": "72656193c687f42556a48c6466432c0891fe5d6f353f294be548a23febef07e3",
}
BLUE, ORANGE, INK = "#28608a", "#b86520", "#252a30"
plt.rcParams.update({"font.size": 10, "axes.labelsize": 11,
                     "axes.spines.top": False, "axes.spines.right": False})


def crossing_number(name):
    match = re.search(r"\d+", str(name))
    return int(match.group()) if match else None


def is_alternating(name):
    n = crossing_number(name)
    if n >= 11:
        return bool(re.match(r"^\d+a", str(name)))
    index = int(str(name).split("_")[1])
    return index <= {8: 18, 9: 41, 10: 123}.get(n, 10**9)


def load_database(filename, table, foreign_key):
    derivation = None
    if filename == "knot_research.db":
        cleanup = json.loads((ROOT / "DATABASE-CLEANUP.json").read_text(encoding="utf-8"))
        if (cleanup["source_commit"] != SOURCE_COMMIT
                or cleanup["source_sha256"] != SOURCE_SHA256[filename]
                or cleanup["remaining_duplicate_groups"] != 0):
            raise ValueError("Unexpected cleanup provenance")
        path = CLEANED_KNOTS
        expected_digest = cleanup["output_sha256"]
        derivation = {"local_file": Path(os.path.relpath(path, ROOT)).as_posix(),
                      "record": "DATABASE-CLEANUP.json",
                      "source_snapshot_sha256": SOURCE_SHA256[filename],
                      "selection_rule": cleanup["selection_rule"]}
    else:
        path = SNAPSHOTS / filename
        if not path.is_file():
            path = REPO_DATABASES / filename
        expected_digest = SOURCE_SHA256[filename]
    if not path.is_file():
        raise FileNotFoundError(f"Required numerical input not found: {path}")
    with path.open("rb") as source_file:
        digest = hashlib.file_digest(source_file, "sha256").hexdigest()
    if digest != expected_digest:
        raise ValueError(f"Database hash mismatch for {filename}: {digest}")
    with sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True) as conn:
        frame = pd.read_sql_query(f"""
            SELECT t.id, t.name,
                   t.khovanov_polynomial IS NOT NULL AS marked_complete,
                   COUNT(b.id) AS stored_blocks,
                   SUM(CASE WHEN b.smallest_nonzero > 0 THEN 1 ELSE 0 END) AS positive_blocks,
                   MIN(CASE WHEN b.smallest_nonzero > 0 THEN b.smallest_nonzero END) AS observed_min
            FROM {table} t LEFT JOIN bidegrees b ON b.{foreign_key}=t.id
            GROUP BY t.id ORDER BY t.id
        """, conn)
        types = dict(conn.execute("SELECT typeof(betti), COUNT(*) FROM bidegrees GROUP BY typeof(betti)"))
        block_counts = dict(zip(("total", "null_gap", "nonpositive_gap", "positive_below_1e7"),
            conn.execute("""SELECT COUNT(*), SUM(smallest_nonzero IS NULL),
                 SUM(smallest_nonzero <= 0), SUM(smallest_nonzero > 0 AND smallest_nonzero < 1e-7)
                 FROM bidegrees""").fetchone()))
        duplicate_groups = conn.execute(f"""SELECT COUNT(*) FROM (
            SELECT {foreign_key},h,q FROM bidegrees GROUP BY {foreign_key},h,q HAVING COUNT(*)>1)
        """).fetchone()[0]
        if table == "twisted_unknots":
            crossings = dict(conn.execute("SELECT name,crossings FROM twisted_unknots"))
            frame["crossings"] = frame.name.map(crossings)
        else:
            frame["crossings"] = frame.name.map(crossing_number)
        frontier_blocks = pd.read_sql_query(f"""SELECT t.name,b.h,b.q,b.dimension,b.smallest_nonzero,
            t.khovanov_polynomial IS NOT NULL AS marked_complete
            FROM {table} t JOIN bidegrees b ON b.{foreign_key}=t.id
            WHERE b.smallest_nonzero > 0 ORDER BY b.smallest_nonzero LIMIT 12""", conn)
    frame["database"] = filename
    info = {
        "file": filename, "source_commit": SOURCE_COMMIT,
        "source_url": f"https://raw.githubusercontent.com/adlauda/Hodge_Khovanov_spectra/{SOURCE_COMMIT}/databases/{filename}",
        "source_sha256": SOURCE_SHA256[filename],
        "working_file": Path(os.path.relpath(path, ROOT)).as_posix(),
        "working_sha256": digest,
        "sha256": digest,
        "bytes": path.stat().st_size, "registered_diagrams": int(len(frame)),
        "marked_complete": int(frame.marked_complete.sum()),
        "with_any_blocks": int((frame.stored_blocks > 0).sum()),
        "partial_with_blocks": int(((frame.marked_complete == 0) & (frame.stored_blocks > 0)).sum()),
        "with_no_blocks": frame.loc[frame.stored_blocks == 0, "name"].tolist(),
        "betti_storage_types": types, "block_counts": block_counts,
        "duplicate_bidegree_groups": duplicate_groups,
        "lowest_stored_blocks": frontier_blocks.to_dict("records"),
    }
    if derivation:
        info["local_derivative"] = derivation
    return frame, info


def summarize(frame):
    records = []
    for n, all_n in frame.groupby("crossings"):
        complete = all_n[(all_n.marked_complete == 1) & all_n.observed_min.notna()]
        values = complete.observed_min.to_numpy(float)
        q10, q25, q50, q75, q90 = np.quantile(values, [.1, .25, .5, .75, .9])
        minimum = float(values.min())
        ties = complete.loc[np.isclose(complete.observed_min, minimum, rtol=0, atol=1e-9), "name"].tolist()
        records.append(dict(crossings=int(n), registered=int(len(all_n)),
            complete_used=int(len(complete)), partial_with_blocks=int(((all_n.marked_complete == 0) & (all_n.stored_blocks > 0)).sum()),
            no_blocks=int((all_n.stored_blocks == 0).sum()),
            q10=float(q10), q25=float(q25), median=float(q50), q75=float(q75), q90=float(q90),
            minimum=minimum, all_stored_minimum=float(all_n.observed_min.min()), minimizers="; ".join(ties)))
    return pd.DataFrame(records)


def style_axis(ax):
    ax.set_yscale("log")
    ax.set_xlabel("Diagram crossing number")
    ax.set_ylabel("Smallest reported positive gap\nper diagram")
    ax.grid(True, axis="y", which="both", alpha=.2, linewidth=.6)


def distribution_plot(summary, kind, filename):
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    x = summary.crossings.to_numpy()
    ax.fill_between(x, summary.q10, summary.q90, color=BLUE, alpha=.12, label="10th–90th percentiles")
    ax.fill_between(x, summary.q25, summary.q75, color=BLUE, alpha=.25, label="25th–75th percentiles")
    ax.plot(x, summary["median"], "o-", color=BLUE, ms=4.5, lw=1.5, label="Median")
    ax.plot(x, summary.minimum, "v--", color=INK, ms=5, lw=1, label="Smallest observed value")
    ax.set_xticks(x)
    style_axis(ax)
    if kind == "knots":
        through13 = summary[summary.crossings <= 13]
        at14 = summary[summary.crossings == 14]
        ax.set_title(f"Completed diagrams: {through13.complete_used.sum():,} through 13 crossings; "
                     f"{at14.complete_used.sum():,} at 14", fontsize=9, pad=11)
        ax.axvline(13.5, color=".55", lw=.8, ls=":")
        ax.text(14, .985, "sample", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8)
    else:
        ax.set_title("Completed diagrams: 1,424 links", fontsize=9, pad=11)
    ax.legend(loc="lower left", fontsize=9, frameon=True, framealpha=.95)
    fig.tight_layout()
    fig.savefig(FIGURES / filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def twisted_plot(frame):
    d = frame[(frame.marked_complete == 1) & frame.observed_min.notna()].sort_values("crossings")
    n = d.crossings.to_numpy(float)
    grid = np.linspace(n.min(), n.max(), 250)
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    ax.plot(grid, 4 * np.sin(np.pi / (2 * (grid + 2)))**2, color=BLUE, lw=1.8,
            label="Exact gap in bidegree $(0,3-n)$")
    ax.plot(n, d.observed_min, "o", mfc="white", mec=INK, mew=1.3, ms=5,
            label="Smallest reported gap over stored bidegrees")
    ax.set_yscale("log")
    ax.set_xlabel("Number of twists $n$")
    ax.set_ylabel("Spectral gap")
    ax.set_xticks(n)
    ax.grid(True, which="both", alpha=.2, linewidth=.6)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES / "twisted_unknot_bidegree_gap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def alternating_plot(frame):
    d = frame[frame.name.map(is_alternating)]
    summary = summarize(d)
    fig, ax = plt.subplots(figsize=(7.5, 5.4))
    complete = d[(d.marked_complete == 1) & d.observed_min.notna()]
    jitter = np.random.default_rng(20260907).uniform(-.09, .09, len(complete))
    ax.scatter(complete.crossings + jitter, complete.observed_min, s=6, alpha=.18, color=".45", lw=0,
               label="Stored per-diagram minima; complete records only")
    ax.plot(summary.crossings, summary.minimum, "*", ms=10, color=INK, mfc="white", mew=1.2,
            label="Smallest observed value at each crossing number")
    odd = np.arange(3, 14, 2)
    even = np.arange(4, 15, 2)
    ax.plot(odd, 4*np.sin(np.pi/(2*odd))**2, "-o", color=BLUE, lw=1.2, ms=3,
            label=r"$T(2,n)$: exact gap in $(n-1,3n-2)$")
    ax.plot(even, 4*np.sin(np.pi/(2*(even-1)))**2, "-s", color=ORANGE, lw=1.2, ms=3,
            label=r"$\mathrm{Tw}_n$: exact gap in $(n-2,3n-9)$")
    exact = 4*np.sin(np.pi/26)**2
    ax.plot([14], [exact], "D", mfc="white", mec=ORANGE, mew=1.5, ms=7)
    ax.annotate("Formula for the unsampled twist diagram", xy=(14, exact), xytext=(8.2, .041),
                fontsize=8, arrowprops={"arrowstyle": "-", "color": ".5"})
    ax.axvline(13.5, color=".55", lw=.8, ls=":")
    ax.text(14, .985, "sample", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8)
    style_axis(ax)
    ax.set_ylim(bottom=.035)
    ax.set_xticks(range(3, 15))
    ax.legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(.5, -.15), ncol=1)
    fig.tight_layout()
    fig.savefig(FIGURES / "alternating_distinguished_bidegrees.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return summary


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    configs = [("knot_research.db", "knots", "knot_id"),
               ("knot_research_14.db", "knots", "knot_id"),
               ("link_research.db", "links", "link_id"),
               ("twisted_unknot_research.db", "twisted_unknots", "twisted_unknot_id")]
    frames, manifest = {}, {"repository": SOURCE_REPO, "commit": SOURCE_COMMIT,
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "selection_policy": "Only diagrams marked complete with a positive stored gap enter distribution statistics. Partial diagrams are counted but excluded. No eigensolver was run.",
        "sources": []}
    for config in configs:
        frame, info = load_database(*config)
        frames[config[0]] = frame
        manifest["sources"].append(info)
        frame.to_csv(ROOT / (config[0].replace(".db", "_diagram_summary.csv")), index=False)
    knots = pd.concat([frames["knot_research.db"], frames["knot_research_14.db"]], ignore_index=True)
    summaries = {"knots": summarize(knots), "links": summarize(frames["link_research.db"])}
    for name, summary in summaries.items():
        summary.to_csv(ROOT / f"{name}_gap_summary.csv", index=False)
        distribution_plot(summary, name, f"observed_{name}_gap_distribution.png")
    twisted_plot(frames["twisted_unknot_research.db"])
    alt = alternating_plot(knots)
    alt.to_csv(ROOT / "alternating_gap_summary.csv", index=False)
    manifest["figures"] = [{"file": p.name, "sha256": hashlib.file_digest(p.open("rb"), "sha256").hexdigest()}
                           for p in sorted(FIGURES.glob("*.png")) if p.name in {
                               "observed_knots_gap_distribution.png", "observed_links_gap_distribution.png",
                               "twisted_unknot_bidegree_gap.png", "alternating_distinguished_bidegrees.png"}]
    manifest["auxiliary_sources"] = []
    for name, expected, description in [
        ("snapshots/knotinfo.csv", "f5d0a015177423c67e5fa2c428e2b6524ad04a4f9302c2c3f8ad2737d221340b",
         "Diagram input used only for the separate orientation check; same pinned repository commit."),
        ("frontier_diagnostics.json", None,
         "Conservatively relabeled derivative of the previous local recomputation record; source SHA256 73736d545eb63e15be35f9ed82cd246040260990f88570b19155f4c6be901747; not regenerated by this plotting script.")]:
        p = SNAPSHOTS / "knotinfo.csv" if name == "snapshots/knotinfo.csv" else ROOT / name
        if name == "snapshots/knotinfo.csv" and not p.is_file():
            p = REPO_ROOT / "data/knotinfo.csv"
        if p.exists():
            with p.open("rb") as f:
                digest = hashlib.file_digest(f, "sha256").hexdigest()
            normalized_line_endings = False
            if expected and digest != expected:
                # Git may check out CSV text with CRLF on Windows.
                normalized = hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
                if normalized != expected:
                    raise ValueError(f"Auxiliary snapshot hash mismatch: {name}")
                normalized_line_endings = True
            auxiliary = {"file": Path(os.path.relpath(p, ROOT)).as_posix(),
                         "sha256": digest, "description": description}
            if expected:
                auxiliary["source_sha256"] = expected
                auxiliary["line_endings_normalized_for_source_check"] = normalized_line_endings
            manifest["auxiliary_sources"].append(auxiliary)
    (ROOT / "NUMERICS-MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"coverage": [{k: d[k] for k in ("file", "registered_diagrams", "marked_complete", "with_any_blocks", "partial_with_blocks", "duplicate_bidegree_groups")} for d in manifest["sources"]],
        "figures": manifest["figures"]}, indent=2))


if __name__ == "__main__":
    main()
