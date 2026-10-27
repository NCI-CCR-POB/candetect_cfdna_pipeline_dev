#!/usr/bin/env python3
"""Generate an HTML report for cfDNA WGS ichorCNA and voting outputs."""

import argparse
import base64
import csv
import datetime as dt
import html
import os
import re
import shutil
import subprocess
from pathlib import Path


VOTING_COLUMNS = [
    "library",
    "cohort",
    "Tumor_Fraction",
    "ichor_cutoff",
    "ichor_pred_class",
    "nmf_prob",
    "med_score_cutoff_nmf",
    "med_score_pred_class_nmf",
    "best_n_components_NMF",
    "binwise_prob",
    "med_score_cutoff_binwise",
    "med_score_pred_class_binwise",
    "num_positive_calls",
]

COLUMN_LABELS = {
    "library": "Library",
    "cohort": "Cohort",
    "Tumor_Fraction": "ichorCNA Tumor_Fraction",
    "ichor_cutoff": "ichorCNA Cutoff",
    "ichor_pred_class": "ichorCNA Call",
    "nmf_prob": "NMF Score",
    "med_score_cutoff_nmf": "NMF Cutoff",
    "med_score_pred_class_nmf": "NMF Call",
    "best_n_components_NMF": "NMF Components",
    "binwise_prob": "Binwise Score",
    "med_score_cutoff_binwise": "Binwise Cutoff",
    "med_score_pred_class_binwise": "Binwise Call",
    "num_positive_calls": "Positive Calls",
}

NUMERIC_COLUMNS = {
    "Tumor_Fraction",
    "ichor_cutoff",
    "nmf_prob",
    "med_score_cutoff_nmf",
    "binwise_prob",
    "med_score_cutoff_binwise",
}

METHOD_SUMMARY_COLUMNS = [
    {
        "method": "ichorCNA",
        "source_column": "ichor_pred_class",
        "display_label": "ichor_pred_class",
    },
    {
        "method": "NMF",
        "source_column": "med_score_pred_class_nmf",
        "display_label": "nmf_med_score_pred_class",
    },
    {
        "method": "Binwise",
        "source_column": "med_score_pred_class_binwise",
        "display_label": "binwise_med_score_pred_class",
    },
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build an HTML report from cfDNA WGS ichorCNA and voting outputs."
    )
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Pipeline result directory, for example nf1_training_set_new.",
    )
    parser.add_argument(
        "--output",
        help=(
            "HTML output path. Defaults to "
            "<data-dir>/analysis/cfdna_wgs/report/cfdna_wgs_report.html."
        ),
    )
    parser.add_argument(
        "--tf",
        help="Tumor fraction CSV. Defaults to <data-dir>/analysis/cfdna_wgs/ichorCNA/tf.csv.",
    )
    parser.add_argument(
        "--voting-dir",
        help=(
            "Voting output directory. Defaults to "
            "<data-dir>/analysis/cfdna_wgs/plots/voting."
        ),
    )
    parser.add_argument(
        "--plots-dir",
        help=(
            "Plot output directory. Defaults to "
            "<data-dir>/analysis/cfdna_wgs/plots."
        ),
    )
    parser.add_argument(
        "--ichor-dir",
        help="ichorCNA output directory. Defaults to <data-dir>/analysis/cfdna_wgs/ichorCNA.",
    )
    parser.add_argument(
        "--title",
        default="cfDNA WGS Report",
        help="Report title shown at the top of the HTML page.",
    )
    parser.add_argument(
        "--no-thumbnails",
        action="store_true",
        help="Skip PDF-to-PNG thumbnail generation for ichorCNA PDFs.",
    )
    parser.add_argument(
        "--thumbnail-size",
        type=int,
        default=900,
        help="Maximum pixel size for generated ichorCNA thumbnails.",
    )
    parser.add_argument(
        "--max-table-rows",
        type=int,
        default=0,
        help="Limit rows displayed per voting table. Use 0 to show all rows.",
    )
    return parser.parse_args()


def esc(value):
    return html.escape(str(value), quote=True)


def normalize_library(value):
    return str(value or "").strip().upper()


def safe_id(value):
    safe = re.sub(r"[^A-Za-z0-9_-]+", "-", str(value)).strip("-")
    return safe or "section"


def read_csv_rows(path):
    if not path or not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def format_fraction(value):
    if value is None or value == "":
        return "NA"
    try:
        return f"{float(value):.4f}".rstrip("0").rstrip(".")
    except ValueError:
        return str(value)


def format_table_value(column, value):
    if value is None or value == "":
        return "NA"
    if column in NUMERIC_COLUMNS:
        return format_fraction(value)
    return str(value)


def read_tumor_fraction_map(tf_path):
    rows = read_csv_rows(tf_path)
    tumor_fraction = {}
    for row in rows:
        library = (
            row.get("Sample_ID")
            or row.get("sample_id")
            or row.get("library")
            or row.get("Library")
            or row.get("sample")
        )
        fraction = row.get("Tumor_Fraction") or row.get("tumor_fraction")
        if library:
            tumor_fraction[normalize_library(library)] = fraction
    return rows, tumor_fraction


def rel_url(target, base_dir):
    rel = os.path.relpath(target, base_dir)
    return Path(rel).as_posix()


def image_data_uri(path):
    if not path or not path.exists():
        return ""
    suffix = path.suffix.lower()
    mime = "image/png" if suffix == ".png" else "image/jpeg"
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{encoded}"


def voting_key(path):
    name = path.name
    name = re.sub(r"^roc_calls_heatmap_", "", name)
    name = re.sub(r"_majority_vote\.csv$", "", name)
    name = re.sub(r"\.(png|pdf)$", "", name)
    return name


def auc_key(path):
    name = path.name
    name = re.sub(r"^auc_comparison_", "", name)
    name = re.sub(r"\.(png|csv)$", "", name)
    return name


def roc_key(path):
    name = path.name
    name = re.sub(r"^roc_comparison_", "", name)
    name = re.sub(r"\.png$", "", name)
    return name


def find_voting_outputs(voting_dir):
    pngs = {}
    csvs = {}
    if voting_dir.exists():
        for path in sorted(voting_dir.glob("roc_calls_heatmap_*_v_*.png")):
            if not path.name.startswith("._"):
                pngs[voting_key(path)] = path
        for path in sorted(voting_dir.glob("roc_calls_heatmap_*_v_*_majority_vote.csv")):
            if not path.name.startswith("._"):
                csvs[voting_key(path)] = path
    keys = sorted(set(pngs) | set(csvs))
    return [{"key": key, "png": pngs.get(key), "csv": csvs.get(key)} for key in keys]


def find_auc_outputs(plots_dir):
    pngs = {}
    csvs = {}
    if plots_dir.exists():
        for path in sorted(plots_dir.glob("auc_comparison_*_v_*.png")):
            if not path.name.startswith("._"):
                pngs[auc_key(path)] = path
        for path in sorted(plots_dir.glob("auc_comparison_*_v_*.csv")):
            if not path.name.startswith("._"):
                csvs[auc_key(path)] = path
    keys = sorted(set(pngs) | set(csvs))
    return [{"key": key, "png": pngs.get(key), "csv": csvs.get(key)} for key in keys]


def find_roc_outputs(plots_dir):
    pngs = {}
    if plots_dir.exists():
        for path in sorted(plots_dir.glob("roc_comparison_*_v_*.png")):
            if not path.name.startswith("._"):
                pngs[roc_key(path)] = path
    return [{"key": key, "png": pngs.get(key)} for key in sorted(pngs)]


def find_ichor_pdfs(ichor_dir):
    if not ichor_dir.exists():
        return []
    return sorted(
        path
        for path in ichor_dir.rglob("*_genomeWide.pdf")
        if path.is_file() and not path.name.startswith("._")
    )


def find_pdf_converter():
    for converter in ("pdftoppm", "magick", "convert"):
        if shutil.which(converter):
            return converter
    return None


def make_pdf_thumbnail(pdf_path, thumb_dir, converter, thumbnail_size):
    if not converter:
        return None

    thumb_dir.mkdir(parents=True, exist_ok=True)
    thumb_path = thumb_dir / f"{pdf_path.stem}.png"
    if thumb_path.exists() and thumb_path.stat().st_mtime >= pdf_path.stat().st_mtime:
        return thumb_path

    if converter == "pdftoppm":
        prefix = thumb_dir / pdf_path.stem
        command = [
            "pdftoppm",
            "-f",
            "1",
            "-singlefile",
            "-png",
            "-scale-to",
            str(thumbnail_size),
            str(pdf_path),
            str(prefix),
        ]
    elif converter == "magick":
        command = [
            "magick",
            "-density",
            "120",
            f"{pdf_path}[0]",
            "-thumbnail",
            f"{thumbnail_size}x{thumbnail_size}",
            str(thumb_path),
        ]
    else:
        command = [
            "convert",
            "-density",
            "120",
            f"{pdf_path}[0]",
            "-thumbnail",
            f"{thumbnail_size}x{thumbnail_size}",
            str(thumb_path),
        ]

    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (OSError, subprocess.CalledProcessError):
        return None

    return thumb_path if thumb_path.exists() else None


def render_metric(value, label):
    return f"""
      <div class="metric">
        <strong>{esc(value)}</strong>
        <span>{esc(label)}</span>
      </div>
    """


def render_warning_list(warnings):
    if not warnings:
        return ""
    items = "\n".join(f"<li>{esc(warning)}</li>" for warning in warnings)
    return f"""
    <section class="panel warning-panel">
      <h2>Report Warnings</h2>
      <ul>{items}</ul>
    </section>
    """


def comparison_cohorts(comparison_key, rows):
    if "_v_" in comparison_key:
        cohort1, cohort2 = comparison_key.split("_v_", 1)
        return cohort1, cohort2

    cohorts = []
    for row in rows:
        cohort = str(row.get("cohort", "")).strip()
        if cohort and cohort not in cohorts:
            cohorts.append(cohort)
    if len(cohorts) >= 2:
        return cohorts[0], cohorts[1]
    return "cohort1", "cohort2"


def parse_binary_call(value):
    if value is None:
        return None
    text = str(value).strip().strip('"')
    if text in {"0", "0.0"}:
        return 0
    if text in {"1", "1.0"}:
        return 1
    try:
        numeric = float(text)
    except ValueError:
        return None
    if numeric == 0:
        return 0
    if numeric == 1:
        return 1
    return None


def prediction_from_call(call, cohort1, cohort2):
    if call == 0:
        return cohort1
    if call == 1:
        return cohort2
    return None


def prediction_from_positive_calls(value, cohort1, cohort2):
    if value is None or value == "":
        return None
    try:
        positive_calls = float(value)
    except ValueError:
        return None
    return cohort2 if positive_calls >= 2 else cohort1


def accuracy_text(correct, total):
    if total <= 0:
        return "NA"
    return f"{correct / total * 100:.1f}%"


def count_text(correct, total):
    return f"{correct}/{total}"


def accuracy_cell(correct, total):
    if total <= 0:
        return "NA"
    return f"{count_text(correct, total)} ({accuracy_text(correct, total)})"


def summarize_predictions(rows, cohort1, cohort2, predictor):
    summary = {
        cohort1: {"correct": 0, "total": 0},
        cohort2: {"correct": 0, "total": 0},
    }

    for row in rows:
        actual = str(row.get("cohort", "")).strip()
        if actual not in summary:
            continue
        summary[actual]["total"] += 1
        predicted = predictor(row)
        if predicted == actual:
            summary[actual]["correct"] += 1

    correct = summary[cohort1]["correct"] + summary[cohort2]["correct"]
    total = summary[cohort1]["total"] + summary[cohort2]["total"]
    return {
        "cohort1": summary[cohort1],
        "cohort2": summary[cohort2],
        "correct": correct,
        "total": total,
    }


def cohort_count_cell(cohort, stats):
    return f"{cohort}: {count_text(stats['correct'], stats['total'])}"


def cohort_accuracy_cell(cohort, stats):
    return f"{cohort}: {accuracy_cell(stats['correct'], stats['total'])}"


def summarize_prediction_rows(comparison_key, rows):
    cohort1, cohort2 = comparison_cohorts(comparison_key, rows)
    comparison_label = comparison_key.replace("_v_", " v ")
    matrix_row = {"Comparison": comparison_label}
    method_rows = []

    for method in METHOD_SUMMARY_COLUMNS:
        source_column = method["source_column"]
        stats = summarize_predictions(
            rows=rows,
            cohort1=cohort1,
            cohort2=cohort2,
            predictor=lambda row, column=source_column: prediction_from_call(
                parse_binary_call(row.get(column)), cohort1, cohort2
            ),
        )
        matrix_row[method["method"]] = accuracy_cell(stats["correct"], stats["total"])
        method_rows.append(
            {
                "Comparison": comparison_label,
                "Method": method["method"],
                "Cohort 1 Correct/Total": cohort_accuracy_cell(cohort1, stats["cohort1"]),
                "Cohort 2 Correct/Total": cohort_accuracy_cell(cohort2, stats["cohort2"]),
                "Total Correct/Total": accuracy_cell(stats["correct"], stats["total"]),
            }
        )

    majority_stats = summarize_predictions(
        rows=rows,
        cohort1=cohort1,
        cohort2=cohort2,
        predictor=lambda row: prediction_from_positive_calls(
            row.get("num_positive_calls"), cohort1, cohort2
        ),
    )
    matrix_row["Majority Vote"] = accuracy_cell(
        majority_stats["correct"], majority_stats["total"]
    )
    majority_row = {
        "Comparison": comparison_label,
        "Cohort 1 Correct/Total": cohort_accuracy_cell(
            cohort1, majority_stats["cohort1"]
        ),
        "Cohort 2 Correct/Total": cohort_accuracy_cell(
            cohort2, majority_stats["cohort2"]
        ),
        "Total Correct/Total": accuracy_cell(
            majority_stats["correct"], majority_stats["total"]
        ),
    }

    return {
        "matrix": matrix_row,
        "methods": method_rows,
        "majority": majority_row,
    }


def render_static_table(rows, columns):
    if not rows:
        return '<p class="muted">No rows available.</p>'

    header = "".join(f"<th>{esc(column)}</th>" for column in columns)
    body_rows = []
    for row in rows:
        cells = "".join(f"<td>{esc(row.get(column, 'NA'))}</td>" for column in columns)
        body_rows.append(f"<tr>{cells}</tr>")

    return f"""
      <div class="table-wrap summary-table-wrap">
        <table>
          <thead><tr>{header}</tr></thead>
          <tbody>{''.join(body_rows)}</tbody>
        </table>
      </div>
    """


def render_prediction_summary(voting_outputs):
    matrix_rows = []
    method_rows = []
    majority_rows = []

    for item in voting_outputs:
        if not item["csv"] or not item["csv"].exists():
            continue
        rows = read_csv_rows(item["csv"])
        summary = summarize_prediction_rows(item["key"], rows)
        matrix_rows.append(summary["matrix"])
        method_rows.extend(summary["methods"])
        majority_rows.append(summary["majority"])

    if not matrix_rows:
        return """
        <section class="panel prediction-summary">
          <h2>Prediction Summary</h2>
          <p class="muted">No majority-vote CSV files were found for prediction summaries.</p>
        </section>
        """

    matrix_columns = ["Comparison", "ichorCNA", "NMF", "Binwise", "Majority Vote"]
    method_columns = [
        "Comparison",
        "Method",
        "Cohort 1 Correct/Total",
        "Cohort 2 Correct/Total",
        "Total Correct/Total",
    ]
    majority_columns = [
        "Comparison",
        "Cohort 1 Correct/Total",
        "Cohort 2 Correct/Total",
        "Total Correct/Total",
    ]

    return f"""
    <section class="panel prediction-summary">
      <h2>Prediction Summary</h2>
      <p class="muted">Calls use 0 for cohort 1 and 1 for cohort 2. Majority-vote criterion: num_positive_calls &gt;= 2 predicts cohort 2; num_positive_calls &lt; 2 predicts cohort 1.</p>

      <h3>At-A-Glance Accuracy Matrix</h3>
      {render_static_table(matrix_rows, matrix_columns)}

      <h3>Per-Method Prediction Summary</h3>
      {render_static_table(method_rows, method_columns)}

      <h3>Overall Majority-Vote Prediction Summary</h3>
      {render_static_table(majority_rows, majority_columns)}
    </section>
    """


def render_voting_table(rows, table_id, max_rows):
    if not rows:
        return '<p class="muted">No majority-vote CSV rows found for this comparison.</p>'

    columns = [column for column in VOTING_COLUMNS if column in rows[0]]
    if not columns:
        columns = list(rows[0].keys())

    display_rows = rows if max_rows <= 0 else rows[:max_rows]
    hidden_note = ""
    if max_rows > 0 and len(rows) > max_rows:
        hidden_note = (
            f'<p class="muted">Showing {max_rows} of {len(rows)} rows. '
            "Increase --max-table-rows to show more.</p>"
        )

    header = "".join(f"<th>{esc(COLUMN_LABELS.get(col, col))}</th>" for col in columns)
    body_rows = []
    for row in display_rows:
        search_text = " ".join(str(row.get(col, "")) for col in columns)
        cells = "".join(
            f"<td>{esc(format_table_value(col, row.get(col, '')))}</td>" for col in columns
        )
        body_rows.append(f'<tr data-search="{esc(search_text.lower())}">{cells}</tr>')

    return f"""
      <div class="table-toolbar">
        <label>
          Filter table
          <input type="search" placeholder="LIB252, MPNST, 3" oninput="filterTable('{esc(table_id)}', this.value)">
        </label>
      </div>
      {hidden_note}
      <div class="table-wrap">
        <table id="{esc(table_id)}">
          <thead><tr>{header}</tr></thead>
          <tbody>{''.join(body_rows)}</tbody>
        </table>
      </div>
    """


def render_voting_sections(voting_outputs, output_dir, max_table_rows):
    if not voting_outputs:
        return '<p class="muted">No voting heatmaps or majority-vote CSV files were found.</p>'

    sections = []
    for item in voting_outputs:
        key = item["key"]
        label = key.replace("_v_", " v ")
        table_id = f"table-{safe_id(key)}"
        rows = read_csv_rows(item["csv"]) if item["csv"] else []

        if item["png"] and item["png"].exists():
            heatmap = (
                f'<img class="heatmap" src="{esc(image_data_uri(item["png"]))}" '
                f'alt="Voting heatmap for {esc(label)}">'
            )
        else:
            heatmap = '<div class="image-placeholder">No voting PNG found</div>'

        csv_link = ""
        if item["csv"] and item["csv"].exists():
            csv_link = (
                f'<a class="small-link" href="{esc(rel_url(item["csv"], output_dir))}">'
                "Open CSV</a>"
            )

        sections.append(
            f"""
            <section class="panel voting-panel">
              <div class="section-heading">
                <div>
                  <h2>Voting: {esc(label)}</h2>
                  <p class="muted">{len(rows)} samples in majority-vote table</p>
                </div>
                {csv_link}
              </div>
              {heatmap}
              <h3>Majority-Vote Calls</h3>
              {render_voting_table(rows, table_id, max_table_rows)}
            </section>
            """
        )
    return "\n".join(sections)


def render_auc_sections(auc_outputs, output_dir):
    if not auc_outputs:
        return """
        <section class="panel">
          <h2>AUC Comparison Plots</h2>
          <p class="muted">No AUC comparison PNG files were found.</p>
        </section>
        """

    cards = []
    for item in auc_outputs:
        key = item["key"]
        label = key.replace("_v_", " v ")
        if item["png"] and item["png"].exists():
            plot = (
                f'<img class="heatmap auc-plot" src="{esc(image_data_uri(item["png"]))}" '
                f'alt="AUC comparison plot for {esc(label)}">'
            )
        else:
            plot = '<div class="image-placeholder">No AUC PNG found</div>'

        csv_link = ""
        if item["csv"] and item["csv"].exists():
            csv_link = (
                f'<a class="small-link" href="{esc(rel_url(item["csv"], output_dir))}">'
                "Open AUC CSV</a>"
            )

        cards.append(
            f"""
            <article class="auc-card">
              <div class="section-heading">
                <div>
                  <h3>{esc(label)}</h3>
                </div>
                {csv_link}
              </div>
              {plot}
            </article>
            """
        )

    return f"""
    <section class="panel auc-panel">
      <h2>AUC Comparison Plots</h2>
      <p class="muted">AUC comparison PNGs are embedded in this HTML. CSV links point back to the original result files.</p>
      <div class="auc-grid">
        {''.join(cards)}
      </div>
    </section>
    """


def render_roc_sections(roc_outputs, output_dir):
    if not roc_outputs:
        return """
        <section class="panel">
          <h2>ROC Curve Plots</h2>
          <p class="muted">No ROC comparison PNG files were found.</p>
        </section>
        """

    cards = []
    for item in roc_outputs:
        key = item["key"]
        label = key.replace("_v_", " v ")
        if item["png"] and item["png"].exists():
            plot = (
                f'<img class="heatmap auc-plot" src="{esc(image_data_uri(item["png"]))}" '
                f'alt="ROC curve plot for {esc(label)}">'
            )
        else:
            plot = '<div class="image-placeholder">No ROC PNG found</div>'

        png_link = ""
        if item["png"] and item["png"].exists():
            png_link = (
                f'<a class="small-link" href="{esc(rel_url(item["png"], output_dir))}">'
                "Open ROC PNG</a>"
            )

        cards.append(
            f"""
            <article class="auc-card">
              <div class="section-heading">
                <div>
                  <h3>{esc(label)}</h3>
                </div>
                {png_link}
              </div>
              {plot}
            </article>
            """
        )

    return f"""
    <section class="panel auc-panel">
      <h2>ROC Curve Plots</h2>
      <p class="muted">ROC curve PNGs are embedded in this HTML.</p>
      <div class="auc-grid">
        {''.join(cards)}
      </div>
    </section>
    """


def render_ichor_cards(ichor_pdfs, tumor_fraction, output_dir, thumb_dir, converter, args):
    if not ichor_pdfs:
        return '<p class="muted">No ichorCNA genome-wide PDFs were found.</p>'

    def tumor_fraction_sort_key(pdf_path):
        library = pdf_path.stem.replace("_genomeWide", "")
        library_key = normalize_library(library)
        raw_fraction = tumor_fraction.get(library_key)
        try:
            fraction = float(raw_fraction)
        except (TypeError, ValueError):
            return (1, 0, library_key)
        return (0, -fraction, library_key)

    cards = []
    for pdf_path in sorted(ichor_pdfs, key=tumor_fraction_sort_key):
        library = pdf_path.stem.replace("_genomeWide", "")
        library_key = normalize_library(library)
        fraction = format_fraction(tumor_fraction.get(library_key))
        thumb_path = None
        if not args.no_thumbnails:
            thumb_path = make_pdf_thumbnail(
                pdf_path=pdf_path,
                thumb_dir=thumb_dir,
                converter=converter,
                thumbnail_size=args.thumbnail_size,
            )

        if thumb_path:
            thumb_html = (
                f'<img class="thumbnail" src="{esc(image_data_uri(thumb_path))}" '
                f'alt="ichorCNA thumbnail for {esc(library_key)}">'
            )
        else:
            thumb_html = '<div class="thumbnail placeholder">PDF</div>'

        cards.append(
            f"""
            <article class="ichor-card" data-search="{esc(library_key.lower())}">
              <a href="{esc(rel_url(pdf_path, output_dir))}">{thumb_html}</a>
              <div class="ichor-meta">
                <strong>{esc(library_key)}</strong>
                <span>Tumor fraction: {esc(fraction)}</span>
                <a href="{esc(rel_url(pdf_path, output_dir))}">Open PDF</a>
              </div>
            </article>
            """
        )
    return "\n".join(cards)


def build_html(
    args,
    paths,
    counts,
    warnings,
    prediction_summary_html,
    auc_html,
    roc_html,
    voting_html,
    ichor_html,
    converter,
):
    generated_at = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    converter_text = "disabled" if args.no_thumbnails else (converter or "not available")
    metrics = "\n".join(
        [
            render_metric(counts["samples"], "Samples in tf.csv"),
            render_metric(counts["voting"], "Voting comparisons"),
            render_metric(counts["ichor"], "ichorCNA PDFs"),
            render_metric(converter_text, "PDF thumbnail converter"),
        ]
    )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{esc(args.title)}</title>
  <style>
    :root {{
      --bg: #f6f7f9;
      --panel: #ffffff;
      --text: #1f2933;
      --muted: #5c6672;
      --border: #d9dee5;
      --accent: #0f766e;
      --accent-2: #2563eb;
      --warn: #9a3412;
    }}
    * {{
      box-sizing: border-box;
    }}
    body {{
      margin: 0;
      font-family: Arial, Helvetica, sans-serif;
      color: var(--text);
      background: var(--bg);
      line-height: 1.45;
    }}
    header {{
      background: #ffffff;
      border-bottom: 1px solid var(--border);
      padding: 24px 32px;
    }}
    main {{
      max-width: 1320px;
      margin: 0 auto;
      padding: 28px 24px 48px;
    }}
    h1, h2, h3 {{
      margin: 0;
      letter-spacing: 0;
    }}
    h1 {{
      font-size: 30px;
    }}
    h2 {{
      font-size: 20px;
    }}
    h3 {{
      font-size: 16px;
      margin: 18px 0 10px;
    }}
    a {{
      color: var(--accent-2);
      text-decoration: none;
      font-weight: 600;
    }}
    a:hover {{
      text-decoration: underline;
    }}
    .subtitle, .muted {{
      color: var(--muted);
    }}
    .subtitle {{
      margin-top: 6px;
    }}
    .summary {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 12px;
      margin-bottom: 24px;
    }}
    .metric, .panel {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 6px;
    }}
    .metric {{
      padding: 14px 16px;
    }}
    .metric strong {{
      display: block;
      font-size: 25px;
      color: var(--accent);
      overflow-wrap: anywhere;
    }}
    .metric span {{
      color: var(--muted);
      font-size: 13px;
    }}
    .panel {{
      padding: 18px;
      margin-bottom: 22px;
    }}
    .section-heading {{
      display: flex;
      align-items: flex-start;
      justify-content: space-between;
      gap: 16px;
      margin-bottom: 14px;
    }}
    .small-link {{
      white-space: nowrap;
      font-size: 13px;
    }}
    .paths code {{
      display: block;
      margin-top: 6px;
      padding: 8px 10px;
      background: #eef2f6;
      border-radius: 4px;
      overflow-wrap: anywhere;
      font-size: 12px;
    }}
    .warning-panel {{
      border-color: #fed7aa;
      color: var(--warn);
    }}
    .heatmap {{
      display: block;
      max-width: 100%;
      height: auto;
      border: 1px solid var(--border);
      border-radius: 4px;
      background: #ffffff;
      margin: 10px 0 16px;
    }}
    .auc-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(min(360px, 100%), 1fr));
      gap: 18px;
      margin-top: 12px;
    }}
    .auc-card {{
      min-width: 0;
      border-top: 1px solid var(--border);
      padding-top: 12px;
    }}
    .auc-plot {{
      margin-bottom: 0;
    }}
    .image-placeholder {{
      height: 220px;
      border: 1px dashed var(--border);
      border-radius: 4px;
      display: flex;
      align-items: center;
      justify-content: center;
      color: var(--muted);
      background: #f9fafb;
      margin: 10px 0 16px;
    }}
    .table-toolbar {{
      display: flex;
      justify-content: flex-end;
      margin-bottom: 8px;
    }}
    label {{
      color: var(--muted);
      font-size: 13px;
    }}
    input[type="search"] {{
      margin-left: 8px;
      padding: 7px 9px;
      min-width: 220px;
      border: 1px solid var(--border);
      border-radius: 4px;
      font-size: 13px;
    }}
    .table-wrap {{
      overflow: auto;
      border: 1px solid var(--border);
      border-radius: 4px;
      max-height: 420px;
    }}
    .summary-table-wrap {{
      max-height: none;
      margin-bottom: 18px;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 12px;
      background: #ffffff;
    }}
    th, td {{
      padding: 6px 8px;
      border-bottom: 1px solid var(--border);
      text-align: left;
      white-space: nowrap;
    }}
    th {{
      position: sticky;
      top: 0;
      background: #eef2f6;
      z-index: 1;
    }}
    .ichor-toolbar {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      gap: 12px;
      margin: 14px 0;
    }}
    .ichor-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(190px, 1fr));
      gap: 14px;
    }}
    .ichor-card {{
      border: 1px solid var(--border);
      border-radius: 6px;
      overflow: hidden;
      background: #ffffff;
    }}
    .thumbnail {{
      width: 100%;
      aspect-ratio: 4 / 3;
      object-fit: contain;
      background: #f2f4f7;
      border-bottom: 1px solid var(--border);
      display: flex;
      align-items: center;
      justify-content: center;
      color: var(--muted);
      font-weight: 700;
    }}
    .ichor-meta {{
      padding: 10px;
      display: grid;
      gap: 4px;
      font-size: 13px;
    }}
    footer {{
      color: var(--muted);
      font-size: 12px;
      margin-top: 28px;
    }}
    @media (max-width: 820px) {{
      header {{
        padding: 20px;
      }}
      main {{
        padding: 20px 14px 36px;
      }}
      .summary {{
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }}
      .section-heading, .ichor-toolbar {{
        display: block;
      }}
      input[type="search"] {{
        width: 100%;
        margin: 8px 0 0;
      }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>{esc(args.title)}</h1>
    <div class="subtitle">Generated {esc(generated_at)} from {esc(paths["data_dir"])}</div>
  </header>
  <main>
    <section class="summary">
      {metrics}
    </section>

    {render_warning_list(warnings)}

    <section class="panel paths">
      <h2>Input Locations</h2>
      <code>{esc(paths["tf"])}</code>
      <code>{esc(paths["voting_dir"])}</code>
      <code>{esc(paths["plots_dir"])}</code>
      <code>{esc(paths["ichor_dir"])}</code>
    </section>

    {prediction_summary_html}

    {auc_html}

    {roc_html}

    <section class="panel">
      <h2>Voting Results</h2>
      <p class="muted">Heatmaps are embedded in this HTML. CSV links point back to the original result files.</p>
    </section>
    {voting_html}

    <section class="panel">
      <div class="section-heading">
        <div>
          <h2>ichorCNA Genome-Wide Figures</h2>
          <p class="muted">Thumbnail overviews are embedded in this HTML and sorted by descending tumor fraction. Open PDF links point to the original genome-wide PDF files.</p>
        </div>
      </div>
      <div class="ichor-toolbar">
        <label>
          Filter libraries
          <input type="search" placeholder="LIB252" oninput="filterCards(this.value)">
        </label>
      </div>
      <div class="ichor-grid" id="ichor-grid">
        {ichor_html}
      </div>
    </section>

    <footer>
      Report generated by candetect_cfdna_pipeline/scripts/create_html_report.py.
    </footer>
  </main>
  <script>
    function filterTable(tableId, query) {{
      const table = document.getElementById(tableId);
      if (!table) return;
      const needle = query.trim().toLowerCase();
      table.querySelectorAll("tbody tr").forEach(row => {{
        row.hidden = needle && !row.dataset.search.includes(needle);
      }});
    }}

    function filterCards(query) {{
      const needle = query.trim().toLowerCase();
      document.querySelectorAll(".ichor-card").forEach(card => {{
        card.hidden = needle && !card.dataset.search.includes(needle);
      }});
    }}
  </script>
</body>
</html>
"""


def main():
    args = parse_args()
    data_dir = Path(args.data_dir).expanduser().resolve()
    output = (
        Path(args.output).expanduser().resolve()
        if args.output
        else data_dir / "analysis" / "cfdna_wgs" / "report" / "cfdna_wgs_report.html"
    )
    tf_path = (
        Path(args.tf).expanduser().resolve()
        if args.tf
        else data_dir / "analysis" / "cfdna_wgs" / "ichorCNA" / "tf.csv"
    )
    voting_dir = (
        Path(args.voting_dir).expanduser().resolve()
        if args.voting_dir
        else data_dir / "analysis" / "cfdna_wgs" / "plots" / "voting"
    )
    plots_dir = (
        Path(args.plots_dir).expanduser().resolve()
        if args.plots_dir
        else data_dir / "analysis" / "cfdna_wgs" / "plots"
    )
    ichor_dir = (
        Path(args.ichor_dir).expanduser().resolve()
        if args.ichor_dir
        else data_dir / "analysis" / "cfdna_wgs" / "ichorCNA"
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    thumb_dir = output.parent / "assets" / "ichor_thumbnails"
    output_dir = output.parent

    warnings = []
    if not data_dir.exists():
        raise SystemExit(f"ERROR: data directory does not exist: {data_dir}")
    if not tf_path.exists():
        warnings.append(f"Tumor fraction file was not found: {tf_path}")
    if not voting_dir.exists():
        warnings.append(f"Voting directory was not found: {voting_dir}")
    if not plots_dir.exists():
        warnings.append(f"Plot directory was not found: {plots_dir}")
    if not ichor_dir.exists():
        warnings.append(f"ichorCNA directory was not found: {ichor_dir}")

    tf_rows, tumor_fraction = read_tumor_fraction_map(tf_path)
    voting_outputs = find_voting_outputs(voting_dir)
    auc_outputs = find_auc_outputs(plots_dir)
    roc_outputs = find_roc_outputs(plots_dir)
    ichor_pdfs = find_ichor_pdfs(ichor_dir)

    converter = None if args.no_thumbnails else find_pdf_converter()
    if not args.no_thumbnails and not converter:
        warnings.append(
            "No PDF converter found. Install pdftoppm, ImageMagick magick, or convert "
            "to generate ichorCNA thumbnails. PDF links are still included."
        )

    paths = {
        "data_dir": str(data_dir),
        "tf": str(tf_path),
        "voting_dir": str(voting_dir),
        "plots_dir": str(plots_dir),
        "ichor_dir": str(ichor_dir),
    }
    counts = {
        "samples": len(tf_rows),
        "voting": len(voting_outputs),
        "ichor": len(ichor_pdfs),
    }

    voting_html = render_voting_sections(
        voting_outputs=voting_outputs,
        output_dir=output_dir,
        max_table_rows=args.max_table_rows,
    )
    prediction_summary_html = render_prediction_summary(voting_outputs)
    auc_html = render_auc_sections(
        auc_outputs=auc_outputs,
        output_dir=output_dir,
    )
    roc_html = render_roc_sections(
        roc_outputs=roc_outputs,
        output_dir=output_dir,
    )
    ichor_html = render_ichor_cards(
        ichor_pdfs=ichor_pdfs,
        tumor_fraction=tumor_fraction,
        output_dir=output_dir,
        thumb_dir=thumb_dir,
        converter=converter,
        args=args,
    )

    report_html = build_html(
        args=args,
        paths=paths,
        counts=counts,
        warnings=warnings,
        prediction_summary_html=prediction_summary_html,
        auc_html=auc_html,
        roc_html=roc_html,
        voting_html=voting_html,
        ichor_html=ichor_html,
        converter=converter,
    )
    output.write_text(report_html, encoding="utf-8")

    print(f"Report written to: {output}")
    print(f"Samples in tf.csv: {counts['samples']}")
    print(f"Voting comparisons: {counts['voting']}")
    print(f"ichorCNA PDFs: {counts['ichor']}")
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"  - {warning}")


if __name__ == "__main__":
    main()
