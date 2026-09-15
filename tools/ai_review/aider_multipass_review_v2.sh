#!/usr/bin/env bash
# aider_multipass_review_v2.sh
#
# RUN THIS on the machine where aider/VS Code live (the laptop), talking to
# your desktop's Ollama server over Tailscale. Requires that you already
# ran setup_ollama_personas.sh/.ps1 ON THE DESKTOP to create:
#   bio-coder-7b, uiux-designer-7b, doc-explainer-7b
#
# and that .aider.model.settings.yml and .aider.model.metadata.json are
# sitting in your repo root (aider auto-discovers them there).
#
# Usage (from repo root):
#   export OLLAMA_API_BASE=http://100.108.166.70:11434
#   ./tools/ai_review/aider_multipass_review_v2.sh app/single_cell
#
# NOTE: macOS ships bash 3.2 by default, which has a known limitation:
# under `set -u`, expanding "${array[@]}" on an EMPTY array throws
# "unbound variable" even though the array is legitimately empty. This
# script avoids bash arrays entirely and uses plain string word-splitting
# instead, which works fine on bash 3.2.

set -euo pipefail

MODULE_DIR="${1:?Usage: $0 <path/to/module_dir>}"
export OLLAMA_API_BASE="${OLLAMA_API_BASE:-http://100.108.166.70:11434}"

MAP_TOKENS="${MAP_TOKENS:-2048}"

# Shared files (relative to app/, since we cd into dirname(MODULE_DIR)
# before invoking aider) that each persona should be able to SEE but not
# edit. Update this list as your repo structure evolves.
READONLY_SHARED="project_manager.py"

REPO_ROOT="$(pwd)"
REVIEW_LOG="${REPO_ROOT}/review_notes_$(basename "$MODULE_DIR")_$(date +%Y%m%d_%H%M%S).md"

PARENT_DIR="$(dirname "$MODULE_DIR")"
TARGET_NAME="$(basename "$MODULE_DIR")"

echo "=== Using Ollama at ${OLLAMA_API_BASE} ==="
echo "=== Module: ${MODULE_DIR} (cwd for aider will be: ${PARENT_DIR}) ==="

# Build read-only args as a plain string (bash-3.2 safe), only including
# files that actually exist relative to PARENT_DIR.
READONLY_ARGS=""
for f in $READONLY_SHARED; do
  if [ -f "${PARENT_DIR}/${f}" ]; then
    READONLY_ARGS="${READONLY_ARGS} --read ${f}"
  fi
done

if [ -z "$READONLY_ARGS" ]; then
  echo "(No shared read-only files found under ${PARENT_DIR}/ matching: ${READONLY_SHARED} -- continuing without them.)"
fi

# --- Pass 1: Bio-Coder -- REVIEW ONLY, NO EDITS -----------------------------
echo ""
echo "--- Pass 1/3: Bio-Coder review (non-destructive) ---"
(
  cd "$PARENT_DIR"
  # shellcheck disable=SC2086
  aider \
    --model "ollama_chat/bio-coder-7b" \
    --subtree-only \
    --map-tokens "${MAP_TOKENS}" \
    --no-auto-commits \
    ${READONLY_ARGS} \
    --message "/ask Review every file in ${TARGET_NAME} for bioinformatics and
statistical correctness: cohort/sample handling, edge cases, genomics data
structure misuse, and reproducibility risks. List concrete issues with
file names and line references. Do not propose a full rewrite -- just the
review." \
    "$TARGET_NAME"
) > "$REVIEW_LOG" 2>&1 || true

echo "Bio-Coder review saved to: ${REVIEW_LOG}"
echo "*** Read this before continuing -- it will NOT auto-apply anything. ***"

# --- Pass 2: UI/UX pass -- EDITS APPLIED, auto-committed --------------------
echo ""
echo "--- Pass 2/3: UI/UX pass (editing, auto-commit) ---"
(
  cd "$PARENT_DIR"
  # shellcheck disable=SC2086
  aider \
    --model "ollama_chat/uiux-designer-7b" \
    --subtree-only \
    --map-tokens "${MAP_TOKENS}" \
    --yes-always \
    --auto-commits \
    ${READONLY_ARGS} \
    --message "Review this module strictly for Streamlit UI/UX: control
layout, dashboard state handling, rendering efficiency, label/help-text
clarity, and consistency with the rest of the app's manager/workspace
pattern. Do not change any statistical or bioinformatics logic." \
    "$TARGET_NAME"
)

# --- Pass 3: Doc-Explainer pass -- EDITS APPLIED, auto-committed ------------
echo ""
echo "--- Pass 3/3: Doc-Explainer pass (editing, auto-commit) ---"
(
  cd "$PARENT_DIR"
  # shellcheck disable=SC2086
  aider \
    --model "ollama_chat/doc-explainer-7b" \
    --subtree-only \
    --map-tokens "${MAP_TOKENS}" \
    --yes-always \
    --auto-commits \
    ${READONLY_ARGS} \
    --message "Add or improve plain-language docstrings, inline comments, and
Streamlit help= tooltips so a non-expert user or new collaborator can
understand what each function/control does and why. Do not change any
logic at all -- documentation and comments only." \
    "$TARGET_NAME"
)

echo ""
echo "=== All passes complete for ${MODULE_DIR} ==="
echo "Review the git log/diff before pushing:"
echo "  git log --oneline -5"
echo "  git diff HEAD~2"
echo "Bio-Coder's non-applied review notes are in: ${REVIEW_LOG}"
