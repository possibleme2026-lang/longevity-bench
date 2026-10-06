#!/usr/bin/env python3
"""Enforce that every byte a task ships is traceable to a source.

The failure this check exists to prevent is not a missing file -- it is a
plausible one. A fixture can be typed by hand, labelled "real", and cited to a
reference that does not exist; nothing about the task looks wrong afterwards. So
provenance is declared in a machine-readable file, each shipped file is bound to
its declaration by SHA256, and anything claiming to be real must name a catalogue
entry that a network verifier already confirmed.

Origin is recorded PER FILE rather than per task, because tasks mix them: a
simulated trial table sits next to a hand-written statistical analysis plan, and
calling both "simulated" would be as wrong as calling both "real".

`authoring/provenance/sources.json` (or `environment/data/sources.json`):

    {
      "task": "senolytic-effect-estimation",
      "retrieved_utc": "2026-10-06T00:00:00Z",
      "catalogue": "docs/sources/verified.json",
      "generator": "authoring/provenance/generate_trial.py",
      "files": [
        {"path": "environment/data/trial_data.csv",
         "sha256": "e3b0c442...", "bytes": 12289,
         "origin": "simulated",
         "transform": "24 months of recruitment, seed 20259198"},

        {"path": "environment/data/protocol.json",
         "sha256": "9f2c...", "bytes": 4672,
         "origin": "authored",
         "note": "Hand-written analysis plan. Contains no measured values; every
                  number in it is a design choice made by the task author."}
      ]
    }

Origins and what each one costs:

  real       the file must name a `source_label` present in the catalogue with
             all its declared identifiers resolved, and a `transform` saying how
             the shipped bytes were derived from that source. An unresolved
             citation is not a source.
  simulated  the task must name a committed `generator`, relative to the task
             directory, that produces it. A synthetic fixture with no generator
             is not acceptable.
  authored   a hand-written specification, schema, or schedule. Requires a `note`
             stating that it carries no measured values.

Usage:
    python3 ci_checks/check-data-provenance.py <task_dir> [<task_dir> ...]
    python3 ci_checks/check-data-provenance.py --catalogue docs/sources/verified.json tasks/...
    python3 ci_checks/check-data-provenance.py --summary tasks/...   # origin census only
"""
import hashlib
import json
import os
import re
import sys

DEFAULT_CATALOGUE = "docs/sources/verified.json"
ORIGINS = ("real", "simulated", "authored")

# Fixtures that document the data rather than being data.
EXEMPT_NAMES = {"DATA_SOURCES.md", "README.md", ".gitkeep", "sources.json"}
PROVENANCE_CANDIDATES = (
    "authoring/provenance/sources.json",
    "environment/data/sources.json",
)
PLACEHOLDERS = re.compile(
    r"(TODO|FIXME|XXX|example\.invalid|lorem ipsum|CHANGEME)", re.I)
HEX64 = re.compile(r"^[0-9a-f]{64}$")
MIN_NOTE = 20


def repo_root(start):
    """Walk up from a task directory until the directory holding ci_checks/."""
    path = os.path.abspath(start)
    while True:
        if os.path.isdir(os.path.join(path, "ci_checks")):
            return path
        parent = os.path.dirname(path)
        if parent == path:
            return None
        path = parent


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def data_files(task_dir):
    root = os.path.join(task_dir, "environment", "data")
    if not os.path.isdir(root):
        return []
    out = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for name in filenames:
            if name in EXEMPT_NAMES:
                continue
            out.append(os.path.relpath(os.path.join(dirpath, name), task_dir)
                       .replace(os.sep, "/"))
    return sorted(out)


def load_provenance(task_dir):
    for rel in PROVENANCE_CANDIDATES:
        path = os.path.join(task_dir, rel)
        if os.path.isfile(path):
            with open(path, encoding="utf-8") as fh:
                return rel, json.load(fh)
    return None, None


def catalogue_labels(root, catalogue_rel):
    """Labels in the catalogue whose declared identifiers all resolved."""
    if not root or not catalogue_rel:
        return None, "no catalogue declared"
    path = os.path.join(root, catalogue_rel)
    if not os.path.isfile(path):
        return None, f"catalogue {catalogue_rel} not found"
    try:
        with open(path, encoding="utf-8") as fh:
            entries = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        return None, f"catalogue unreadable ({type(exc).__name__})"
    return {e["label"] for e in entries
            if e.get("all_declared_identifiers_resolve") and e.get("label")}, None


def check_task(task_dir, default_catalogue):
    problems, census = [], {"real": 0, "simulated": 0, "authored": 0}

    files = data_files(task_dir)
    rel, prov = load_provenance(task_dir)

    if not files:
        return problems, census, "no data files"
    if prov is None:
        problems.append(
            f"environment/data/ holds {len(files)} file(s) but no provenance is declared; "
            f"expected one of {', '.join(PROVENANCE_CANDIDATES)} accounting for every one")
        return problems, census, None

    hit = PLACEHOLDERS.search(json.dumps(prov))
    if hit:
        problems.append(f"{rel} contains a placeholder marker {hit.group(0)!r}")

    declared = prov.get("files")
    if not isinstance(declared, list) or not declared:
        problems.append(f"{rel}: files must be a non-empty list")
        return problems, census, None

    by_path = {}
    for entry in declared:
        path = (entry.get("path") or "").replace(os.sep, "/")
        if not path:
            problems.append(f"{rel}: an entry has no path")
            continue
        by_path[path] = entry

    for rel_file in files:
        if rel_file not in by_path:
            problems.append(f"{rel_file} is shipped but not declared in {rel}")
    for path in by_path:
        if path not in files:
            problems.append(f"{rel} declares {path}, which is not in environment/data/")

    needs_generator = any(e.get("origin") == "simulated" for e in by_path.values())
    generator, root = prov.get("generator"), repo_root(task_dir)
    if needs_generator:
        if not generator:
            problems.append(f"{rel}: a file is 'simulated' but no generator is named")
        elif not os.path.isfile(os.path.join(task_dir, generator)):
            # Task-relative, like files[].path -- not repo-relative.
            problems.append(f"{rel}: generator {generator!r} does not exist")

    labels = err = None
    if any(e.get("origin") == "real" for e in by_path.values()):
        labels, err = catalogue_labels(root, prov.get("catalogue") or default_catalogue)
        if err:
            problems.append(f"{rel}: cannot use catalogue: {err}")

    for rel_file in files:
        entry = by_path.get(rel_file)
        if not entry:
            continue
        origin = entry.get("origin")
        if origin not in ORIGINS:
            problems.append(f"{rel_file}: origin must be one of {ORIGINS}, got {origin!r}")
            continue
        census[origin] += 1

        digest = entry.get("sha256")
        if not digest or not HEX64.match(str(digest)):
            problems.append(f"{rel_file}: sha256 missing or malformed in {rel}")
        else:
            actual = sha256_of(os.path.join(task_dir, rel_file))
            if actual != digest:
                problems.append(
                    f"{rel_file}: sha256 does not match the shipped bytes "
                    f"(declared {digest[:12]}..., actual {actual[:12]}...)")

        if origin == "real":
            label = entry.get("source_label")
            if not label:
                problems.append(f"{rel_file}: origin is 'real' but no source_label is given")
            elif labels is not None and label not in labels:
                problems.append(
                    f"{rel_file}: source_label {label!r} is not a fully verified entry in the "
                    f"catalogue; a citation that has not resolved is not a source")
            if not entry.get("transform"):
                problems.append(
                    f"{rel_file}: origin is 'real' but no transform is recorded, so the "
                    f"shipped bytes cannot be traced to what was downloaded")
        elif origin == "authored":
            note = entry.get("note") or ""
            if len(note) < MIN_NOTE:
                problems.append(
                    f"{rel_file}: origin is 'authored' but note is missing or too short; an "
                    f"authored file must say what it is and that it holds no measured values")

    return problems, census, None


def main():
    args = list(sys.argv[1:])
    catalogue, summary = DEFAULT_CATALOGUE, False
    if "--catalogue" in args:
        i = args.index("--catalogue")
        catalogue = args[i + 1]
        del args[i:i + 2]
    if "--summary" in args:
        summary = True
        args.remove("--summary")

    if not args:
        print("usage: check-data-provenance.py [--catalogue PATH] [--summary] <task_dir> ...")
        return 2

    failed, total = 0, {"real": 0, "simulated": 0, "authored": 0}
    for task_dir in args:
        slug = os.path.basename(os.path.abspath(task_dir))
        problems, census, note = check_task(task_dir, catalogue)
        for key in total:
            total[key] += census.get(key, 0)
        if problems:
            failed += 1
            print(f"FAIL {slug}")
            for p in problems:
                print(f"  - {p}")
        elif not summary:
            counts = ", ".join(f"{v} {k}" for k, v in census.items() if v)
            print(f"ok   {slug:<44} {counts or 'no files'}"
                  + (f"  ({note})" if note else ""))

    if failed:
        print(f"\n{failed} task(s) ship data that is not traceable to a source.")
        return 1
    print(f"\nAll {len(args)} task(s) declare traceable provenance for their data.")
    print(f"Suite census: {total['real']} real, {total['simulated']} simulated, "
          f"{total['authored']} authored file(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
