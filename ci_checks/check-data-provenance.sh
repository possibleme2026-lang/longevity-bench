#!/bin/bash
# Thin wrapper so this joins the CHECKS table in static-checks.yml, which invokes every
# static check as `bash <script> <task_dir>`. The analysis needs to hash files and read
# JSON, so it lives in check-data-provenance.py next to this file.
set -uo pipefail
exec python3 "$(dirname "${BASH_SOURCE[0]}")/check-data-provenance.py" "$@"
