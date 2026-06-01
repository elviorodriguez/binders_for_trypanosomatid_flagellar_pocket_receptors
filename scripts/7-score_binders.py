"""
score_designs.py — Miniprotein Binder Design Scoring Tool
==========================================================

Scores miniprotein binder designs based on two structural quality metrics:
  - min_interaction_pae  : minimum Predicted Aligned Error at the binder/target
                           interface (lower = better predicted contact)
  - rmsd_to_rfdiff       : RMSD of the AF2 model vs. the original RFdiffusion
                           backbone (lower = more faithful to the designed fold)

For each design the script runs a rank-weighted scoring scheme:
  1. For every AF2 model rank (1-5) it computes the Euclidean distance from the
     ideal point (0, 0) in (min_interaction_pae, rmsd_to_rfdiff) space.
  2. Each rank contributes its distance multiplied by a predefined weight
     (rank 1 carries the most weight because it is the highest-confidence model).
  3. The raw score is the reciprocal of the total weighted distance, so a smaller
     distance from (0,0) → a higher score → a better design.
  4. All per-rank rows that share a design_id receive the same raw_score; the
     output file collapses to one row per design_id, sorted best-first.

Usage
-----
    python score_designs.py <input_csv> [options]

Positional arguments
--------------------
    input_csv
        Path to the comma-separated input file produced by the design analysis
        pipeline.  Required columns:
            design_id           – unique string identifier for each design
            pae_rank            – integer 1-5 (AF2 model rank)
            min_interaction_pae – float, lower is better
            rmsd_to_rfdiff      – float, lower is better

        Optional columns (not used by the scoring algorithm but accepted):
            target, rank, mean_binder_plddt, mean_binder_intrachain_pae,
            binder_length, target_length, model_file, ptm, iptm

Optional arguments
------------------
    -o / --output PATH
        Explicit path for the output CSV.
        Default: same directory as <input_csv>, named "designs_ranking.csv".

    -s / --sep CHAR
        Column separator used in the input file.  Default: "," (comma).
        Use -s $'\\t' for TSV files.

    -v / --verbose
        Print a short summary table (top-10 designs) after scoring.

    -h / --help
        Show this help message and exit.

Examples
--------
    # Minimal — output goes next to the input file
    python score_designs.py ./results/design_analysis_results.csv

    # Explicit output path
    python score_designs.py data/design_analysis_results.csv -o data/ranked.csv

    # Tab-separated input + verbose summary
    python score_designs.py runs/designs.tsv -s $'\\t' -v

Output
------
    A CSV with two columns (design_id, raw_score), one row per unique design,
    sorted by raw_score descending (best design first).
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from scipy.spatial.distance import euclidean


# ---------------------------------------------------------------------------
# Core scoring logic
# ---------------------------------------------------------------------------

def score_designs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Score miniprotein binder designs based on min_interaction_pae and
    rmsd_to_rfdiff.

    The function calculates the Euclidean distance from the ideal point (0, 0)
    for each AF2 model rank, applies rank-importance weights, and returns the
    input DataFrame augmented with a 'raw_score' column.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain at minimum the columns:
            design_id           – unique design identifier
            pae_rank            – integer 1-5 (AF2 model rank within a design)
            min_interaction_pae – float (lower = better)
            rmsd_to_rfdiff      – float (lower = better)

    Returns
    -------
    pandas.DataFrame
        A copy of *df* with an additional 'raw_score' column.
        All rows belonging to the same design_id share the same raw_score.
        Higher raw_score means a better design.

    Notes
    -----
    Rank weights (contribution of each AF2 model rank to the final score):
        rank 1 → 0.10  (highest-confidence model, largest weight)
        rank 2 → 0.08
        rank 3 → 0.06
        rank 4 → 0.04
        rank 5 → 0.02

    raw_score = 1 / sum_over_ranks(weight_rank * euclidean_distance_rank)

    A design with distance 0 for all ranks would score +inf; in practice this
    never occurs because both metrics are strictly positive for real structures.
    Designs with missing rank data are still scored using the ranks present.
    """

    # Rank → weight mapping
    weights = {1: 0.10, 2: 0.08, 3: 0.06, 4: 0.04, 5: 0.02}

    result_df = df.copy()
    design_scores: dict = {}

    for design_id in df["design_id"].unique():
        design_data = df[df["design_id"] == design_id]

        weighted_sum = 0.0
        total_weight = 0.0

        for rank in range(1, 6):
            rank_data = design_data[design_data["pae_rank"] == rank]

            if not rank_data.empty:
                row = rank_data.iloc[0]
                distance = euclidean(
                    [0, 0],
                    [row["min_interaction_pae"], row["rmsd_to_rfdiff"]],
                )
                weighted_sum += distance * weights[rank]
                total_weight += weights[rank]

        if total_weight > 0 and weighted_sum > 0:
            raw_score = 1.0 / weighted_sum
        else:
            raw_score = 0.0

        design_scores[design_id] = raw_score

    result_df["raw_score"] = result_df["design_id"].map(design_scores)
    return result_df


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="score_designs.py",
        description=(
            "Score miniprotein binder designs using a rank-weighted Euclidean "
            "distance metric and write a ranked CSV."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,  # re-use the module docstring as the extended help
    )

    parser.add_argument(
        "input_csv",
        metavar="input_csv",
        type=Path,
        help="Path to the design analysis CSV file.",
    )
    parser.add_argument(
        "-o", "--output",
        metavar="PATH",
        type=Path,
        default=None,
        help=(
            "Output CSV path.  Defaults to <input_csv_dir>/designs_ranking.csv"
        ),
    )
    parser.add_argument(
        "-s", "--sep",
        metavar="CHAR",
        default=",",
        help="Column separator in the input file (default: ',').",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print a top-10 summary table to stdout after scoring.",
    )
    return parser


def resolve_output_path(input_path: Path, output_arg) -> Path:
    """Return the output path, defaulting to the input file's directory."""
    if output_arg is not None:
        return output_arg
    return input_path.parent / "designs_ranking.csv"


def validate_columns(df: pd.DataFrame, required: list[str]) -> None:
    """Raise SystemExit with a clear message if required columns are missing."""
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(
            f"[ERROR] The input file is missing required column(s): "
            f"{', '.join(missing)}\n"
            f"Found columns: {', '.join(df.columns.tolist())}",
            file=sys.stderr,
        )
        sys.exit(1)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    # --- Validate input path ------------------------------------------------
    input_path: Path = args.input_csv.resolve()
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # --- Load data ----------------------------------------------------------
    print(f"[INFO] Loading: {input_path}")
    try:
        df = pd.read_csv(input_path, sep=args.sep)
    except Exception as exc:
        print(f"[ERROR] Could not read CSV: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Loaded {len(df):,} rows, {df.shape[1]} columns.")

    # --- Validate required columns ------------------------------------------
    required_cols = ["design_id", "pae_rank", "min_interaction_pae", "rmsd_to_rfdiff"]
    validate_columns(df, required_cols)

    # --- Score --------------------------------------------------------------
    print("[INFO] Scoring designs …")
    scored_df = score_designs(df)

    # --- Collapse to one row per design, sorted best-first ------------------
    unique_ranked = (
        scored_df[["design_id", "raw_score"]]
        .drop_duplicates()
        .sort_values("raw_score", ascending=False)
        .reset_index(drop=True)
    )
    print(f"[INFO] Scored {len(unique_ranked):,} unique designs.")

    # --- Write output -------------------------------------------------------
    output_path = resolve_output_path(input_path, args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    unique_ranked.to_csv(output_path, index=False)
    print(f"[INFO] Rankings written to: {output_path}")

    # --- Optional summary ---------------------------------------------------
    if args.verbose:
        top_n = min(10, len(unique_ranked))
        print(f"\nTop {top_n} designs:")
        print(unique_ranked.head(top_n).to_string(index=False))


if __name__ == "__main__":
    main()
