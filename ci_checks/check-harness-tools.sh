#!/bin/bash
# The agent harness needs tmux inside the task image.
#
# harbor's Terminus-style agents run the model in a tmux session, and harbor
# installs tmux at trial setup when `tmux -V` fails. That install is best-effort
# and capped at 240s, and on a slim base image it falls back to compiling tmux
# from source. With three trials building concurrently on this host, one in three
# died during setup with "Failed to start tmux session" -- a trial lost to the
# harness rather than to the model, which corrupts a score in a way nothing in
# the logs calls out.
#
# So tmux must be baked at image build time, where the network and the budget are
# reliable. asciinema is what harbor records the session with; it is reported as a
# warning rather than a failure because recording is optional and a base image may
# legitimately provide only tmux.
#
# A task whose base image already ships tmux may opt out with a comment anywhere
# in environment/Dockerfile:
#
#     # harness-tools: provided by base image
#
# Arguments: task directories (e.g., tasks/my-task) or no args to check all
set -uo pipefail

if [ $# -eq 0 ]; then
    FILES_TO_CHECK=$(find tasks -type f -name "Dockerfile" -path "*/environment/*" 2>/dev/null || true)
else
    FILES_TO_CHECK=""
    for task_dir in "$@"; do
        if [ -d "$task_dir" ] && [ -f "$task_dir/environment/Dockerfile" ]; then
            FILES_TO_CHECK="$FILES_TO_CHECK $task_dir/environment/Dockerfile"
        fi
    done
fi

if [ -z "$FILES_TO_CHECK" ]; then
    echo "No environment Dockerfiles to check"
    exit 0
fi

FAILED=0
for file in $FILES_TO_CHECK; do
    [ -f "$file" ] || continue

    if grep -qE '^\s*#\s*harness-tools:\s*provided' "$file"; then
        echo "ok   $file: harness tools declared as provided by the base image"
        continue
    fi

    # An install line that names tmux: apt/apk/dnf/yum/zypper/microdnf/pacman.
    if grep -qE '(apt-get|apt|apk|dnf|yum|zypper|microdnf|pacman)[^#]*\btmux\b' "$file"; then
        if ! grep -qE '(apt-get|apt|apk|dnf|yum|zypper|microdnf|pacman)[^#]*\basciinema\b' "$file"; then
            echo "Warn $file: tmux is installed but asciinema is not (session recording will be unavailable)"
        fi
        continue
    fi

    echo "FAIL $file: tmux is not installed in the task image."
    echo "      The agent harness runs the model in a tmux session and only falls back to"
    echo "      installing tmux at trial setup, best-effort and capped at 240s. On a slim"
    echo "      base that fallback compiled tmux from source and lost 1 trial in 3 at setup."
    echo "      Add, before the data COPY:"
    echo "        RUN apt-get update \\"
    echo "            && apt-get install -y --no-install-recommends tmux asciinema \\"
    echo "            && rm -rf /var/lib/apt/lists/*"
    FAILED=1
done

exit $FAILED
