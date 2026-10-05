#!/usr/bin/env python3
"""Update README task totals and the generated domain coverage chart."""

import re
from html import escape
from pathlib import Path

# Directories to exclude from task counts (template/demo tasks)
EXCLUDED_DIRS = {"hello-world"}

DOMAINS = [
    ("AGING_BIOLOGY_COUNT", "Aging Biology", "tasks/aging-biology", "#3D8A2D"),
    ("CLINICAL_AGING_COUNT", "Clinical Aging", "tasks/clinical-aging", "#C23637"),
    ("GEROSCIENCE_PHARMACOLOGY_COUNT", "Geroscience Pharmacology", "tasks/geroscience-pharmacology", "#0168A9"),
    ("AGING_QUANTITATIVE_COUNT", "Quantitative Aging", "tasks/aging-quantitative", "#C59527"),
]

DOMAIN_FIELDS = {
    "AGING_BIOLOGY_COUNT": [
        ("CELLULAR_SENESCENCE_COUNT", "Cellular Senescence"),
        ("EPIGENETIC_CLOCKS_COUNT", "Epigenetic Clocks"),
        ("PROTEOSTASIS_COUNT", "Proteostasis"),
        ("MITOCHONDRIAL_FUNCTION_COUNT", "Mitochondrial Function"),
        ("INFLAMMAGING_COUNT", "Inflammaging"),
    ],
    "CLINICAL_AGING_COUNT": [
        ("GERIATRIC_ASSESSMENT_COUNT", "Geriatric Assessment"),
        ("FRAILTY_COUNT", "Frailty"),
        ("MULTIMORBIDITY_COUNT", "Multimorbidity"),
        ("CLINICAL_TRIALS_COUNT", "Clinical Trials"),
    ],
    "GEROSCIENCE_PHARMACOLOGY_COUNT": [
        ("SENOLYTICS_COUNT", "Senolytics"),
        ("GEROPROTECTORS_COUNT", "Geroprotectors"),
        ("DOSE_RESPONSE_COUNT", "Dose & Response"),
    ],
    "AGING_QUANTITATIVE_COUNT": [
        ("SURVIVAL_ANALYSIS_COUNT", "Survival Analysis"),
        ("CAUSAL_INFERENCE_COUNT", "Causal Inference"),
        ("OMICS_INTEGRATION_COUNT", "Omics Integration"),
        ("BIOMARKER_VALIDATION_COUNT", "Biomarker Validation"),
    ],
}


def count_tasks(field_path):
    """Count task directories in a field (dirs containing task.toml)"""
    path = Path(field_path)
    if not path.exists():
        return 0
    return sum(
        1 for d in path.iterdir()
        if (
            d.is_dir()
            and not d.name.startswith('.')
            and d.name not in EXCLUDED_DIRS
            and (d / 'task.toml').exists()
        )
    )


def count_domain_tasks(domain_path):
    """Count task.toml files recursively within a domain."""
    path = Path(domain_path)
    if not path.exists():
        return 0
    return sum(
        1
        for task_file in path.glob("*/*/task.toml")
        if task_file.parent.name not in EXCLUDED_DIRS
    )


def render_coverage_svg(domain_counts, field_counts, total):
    """Render grouped field bars using a shared task-count scale."""
    elements = []
    chart_width = 1200
    summary_height = 44
    block_height = 168
    header_height = 36
    footer_height = header_height
    row_height = 33
    bar_x = 430
    bar_width = 700
    max_count = max(
        10,
        max(
            field_counts[marker]
            for fields in DOMAIN_FIELDS.values()
            for marker, _ in fields
        ),
    )
    elements.extend(
        [
            f'  <rect class="summary-hdr" x="0" y="0" width="{chart_width}" height="{summary_height}"/>',
            f'  <text class="summary-title" x="{chart_width / 2:g}" y="29" text-anchor="middle" font-size="18">LONGEVITY-BENCH 0.1</text>',
        ]
    )

    for domain_index, (marker, label, _, color) in enumerate(DOMAINS):
        y = summary_height + domain_index * block_height
        count = domain_counts[marker]
        elements.extend(
            [
                f'  <rect class="group-hdr" x="0" y="{y}" width="{chart_width}" height="{header_height}" style="fill:{color};fill-opacity:0.10"/>',
                f'  <text class="group-caption" x="12" y="{y + 24}" font-size="17" style="fill:{color}">{escape(label)}</text>',
                f'  <text class="group-count" x="{chart_width - 12}" y="{y + 24}" text-anchor="end" font-size="17" style="fill:{color}">{count}</text>',
            ]
        )

        fields = sorted(
            DOMAIN_FIELDS[marker],
            key=lambda field: -field_counts[field[0]],
        )
        for field_index, (field_marker, field_label) in enumerate(fields):
            field_count = field_counts[field_marker]
            row_y = y + header_height + field_index * row_height
            field_bar_width = round(bar_width * field_count / max_count)
            elements.extend(
                [
                    f'  <rect class="field-row" x="0" y="{row_y}" width="{chart_width}" height="{row_height}"/>',
                    f'  <text class="field-name" x="12" y="{row_y + 22}" font-size="15">{escape(field_label)}</text>',
                    f'  <rect class="bar-track" x="{bar_x}" y="{row_y + 10}" width="{bar_width}" height="13"/>',
                    f'  <rect x="{bar_x}" y="{row_y + 10}" width="{field_bar_width}" height="13" fill="{color}"/>',
                    f'  <text class="field-count" x="{chart_width - 12}" y="{row_y + 22}" text-anchor="end" font-size="15">{field_count}</text>',
                ]
            )

    footer_y = summary_height + len(DOMAINS) * block_height
    chart_height = footer_y + footer_height
    footer_color = "#0F766E"
    elements.extend(
        [
            f'  <rect class="group-hdr" x="0" y="{footer_y}" width="{chart_width}" height="{footer_height}" style="fill:{footer_color};fill-opacity:0.10"/>',
            f'  <text class="group-caption" x="12" y="{footer_y + 24}" font-size="17" style="fill:{footer_color}">Total</text>',
            f'  <text class="group-count" x="{chart_width - 12}" y="{footer_y + 24}" text-anchor="end" font-size="17" style="fill:{footer_color}">{total} TASKS</text>',
        ]
    )

    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {chart_width} {chart_height}" preserveAspectRatio="xMidYMid meet" role="img" aria-labelledby="title description">',
            '  <title id="title">Longevity-Bench task coverage by domain and field</title>',
            f'  <desc id="description">{total} active tasks shown as grouped field bars across four aging and longevity domains.</desc>',
            "  <style>",
            '    text { font-family: "Google Sans Code", ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; fill: #0a0a0a; }',
            "    .summary-title { font-weight: 700; }",
            "    .group-caption { font-weight: 600; letter-spacing: 0.04em; text-transform: uppercase; }",
            "    .group-count { font-weight: 700; }",
            "    .field-name { font-weight: 500; }",
            "    .field-count { font-weight: 700; }",
            "    .group-hdr, .field-row { fill: transparent; }",
            "    .bar-track { fill: #f4f4f5; }",
            "    .summary-hdr { fill: transparent; }",
            "    @media (prefers-color-scheme: dark) {",
            "      text { fill: #f0f6fc; }",
            "      .bar-track { fill: #3f3f46; }",
            "    }",
            "  </style>",
            *elements,
            "</svg>",
            "",
        ]
    )


def main():
    # README marker -> task directory path
    fields = {
        # aging-biology
        'CELLULAR_SENESCENCE_COUNT':      'tasks/aging-biology/cellular-senescence',
        'EPIGENETIC_CLOCKS_COUNT':        'tasks/aging-biology/epigenetic-clocks',
        'PROTEOSTASIS_COUNT':             'tasks/aging-biology/proteostasis',
        'MITOCHONDRIAL_FUNCTION_COUNT':   'tasks/aging-biology/mitochondrial-function',
        'INFLAMMAGING_COUNT':             'tasks/aging-biology/inflammaging',
        # clinical-aging
        'GERIATRIC_ASSESSMENT_COUNT':     'tasks/clinical-aging/geriatric-assessment',
        'FRAILTY_COUNT':                  'tasks/clinical-aging/frailty',
        'MULTIMORBIDITY_COUNT':           'tasks/clinical-aging/multimorbidity',
        'CLINICAL_TRIALS_COUNT':          'tasks/clinical-aging/clinical-trials',
        # geroscience-pharmacology
        'SENOLYTICS_COUNT':               'tasks/geroscience-pharmacology/senolytics',
        'GEROPROTECTORS_COUNT':           'tasks/geroscience-pharmacology/geroprotectors',
        'DOSE_RESPONSE_COUNT':            'tasks/geroscience-pharmacology/dose-response',
        # aging-quantitative
        'SURVIVAL_ANALYSIS_COUNT':        'tasks/aging-quantitative/survival-analysis',
        'CAUSAL_INFERENCE_COUNT':         'tasks/aging-quantitative/causal-inference',
        'OMICS_INTEGRATION_COUNT':        'tasks/aging-quantitative/omics-integration',
        'BIOMARKER_VALIDATION_COUNT':     'tasks/aging-quantitative/biomarker-validation',
    }

    counts = {marker: count_tasks(path) for marker, path in fields.items()}
    domain_counts = {
        marker: count_domain_tasks(path)
        for marker, _, path, _ in DOMAINS
    }
    counts.update(domain_counts)
    total = sum(domain_counts.values())
    counts['TOTAL_COUNT'] = total

    readme_path = Path('README.md')
    content = readme_path.read_text(encoding="utf-8")

    for marker, count in counts.items():
        if marker == 'TOTAL_COUNT':
            pattern = rf'(<!--{marker}-->\s*\*\*)\d+(\*\*)'
            replacement = rf'\g<1>{count}\2'
        else:
            pattern = rf'(<!--{marker}-->)\s*\d+'
            replacement = rf'\g<1> {count}'
        content = re.sub(pattern, replacement, content)

    # Update the badge
    badge_pattern = r'(\[!\[Tasks\]\(https://img\.shields\.io/badge/tasks-)\d+(-white\.svg\)\])'
    content = re.sub(badge_pattern, rf'\g<1>{total}\2', content)

    # newline="\n" is load-bearing: this script rewrites README.md and
    # assets/task-coverage.svg, which are tracked files. Path.write_text defaults
    # to newline=None, which translates "\n" to os.linesep -- on Windows that
    # silently converts both files to CRLF on every run.
    readme_path.write_text(content, encoding="utf-8", newline="\n")
    coverage_path = Path('assets/task-coverage.svg')
    coverage_path.parent.mkdir(parents=True, exist_ok=True)
    coverage_path.write_text(render_coverage_svg(domain_counts, counts, total),
                             encoding="utf-8", newline="\n")

    print("Updated task counts and assets/task-coverage.svg")
    print(f"  Total tasks: {total}")
    for marker, label, _, _ in DOMAINS:
        print(f"  {label}: {domain_counts[marker]}")


if __name__ == '__main__':
    main()
