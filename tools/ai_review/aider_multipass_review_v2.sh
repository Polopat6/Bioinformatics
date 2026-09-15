#!/usr/bin/env bash
# aider_multipass_review_v2.sh
#
# RUN THIS on the machine where aider/VS Code live (the laptop), talking to
# your desktop's Ollama server over Tailscale. Requires that you already
# ran setup_ollama_personas.sh ON THE DESKTOP to create:
#   bio-coder-7b, uiux-designer-7b, doc-explainer-7b
#
# and that .aider.model.settings.yml and .aider.model.metadata.json are
# sitting in your repo root (same directory as this script, or your repo's
# git root -- aider auto-discovers them there).
#
# Usage:
#   export OLLAMA_API_BASE=http://100.108.166.70:11434
#   ./aider_multipass_review_v2.sh app/single_cell
#   ./aider_multipass_review_v2.sh app/bulk_rnaseq
#   ./aider_multipass_review_v2.sh app/ontology_analysis

set -euo pipefail

MODULE_DIR="${1:?Usage: $0 <path/to/module_dir>}"
export OLLAMA_API_BASE="${OLLAMA_API_BASE:-http://100.108.166.70:11434}"

# 16k context is tight for a 7B model doing whole-file rewrites -- keep the
# repo map modest so there's room left for the actual file + reply.
MAP_TOKENS="${MAP_TOKENS:-2048}"

# Shared/common files each persona should be able to SEE (read-only) but
# not edit, so cross-module calls still make sense. Adjust to your repo.
READONLY_SHARED="project_manager.py reference_manager.py"

REVIEW_LOG="review_notes_$(basename "$MODULE_DIR")_$(date +%Y%m%d_%H%M%S).md"

echo "=== Using Ollama at ${OLLAMA_API_BASE} ==="
echo "=== Module: ${MODULE_DIR} ==="

readonly_flags=()
for f in $READONLY_SHARED; do
  [ -f "$f" ] && readonly_flags+=(--read "$f")
done

# --- Pass 1: Bio-Coder -- REVIEW ONLY, NO EDITS -----------------------------
# Uses aider's /ask mode so the 7B model can flag concerns without touching
# files. Treat this as a checklist to review yourself, not an auto-fix.
echo ""
echo "--- Pass 1/3: Bio-Coder review (non-destructive) ---"
aider \
  --model "ollama_chat/bio-coder-7b" \
  --subtree-only \
  --map-tokens "${MAP_TOKENS}" \
  --no-auto-commits \
  "${readonly_flags[@]}" \
  --message "/ask Review every file in ${MODULE_DIR} for bioinformatics and
statistical correctness: cohort/sample handling, edge cases, genomics data
structure misuse, and reproducibility risks. List concrete issues with
file names and line references. Do not propose a full rewrite -- just the
review." \
  "$MODULE_DIR" > "$REVIEW_LOG" 2>&1 || true

echo "Bio-Coder review saved to: ${REVIEW_LOG}"
echo "*** Read this before continuing -- it will NOT auto-apply anything. ***"

# --- Pass 2: UI/UX pass -- EDITS APPLIED, auto-committed --------------------
echo ""
echo "--- Pass 2/3: UI/UX pass (editing, auto-commit) ---"
aider \
  --model "ollama_chat/uiux-designer-7b" \
  --subtree-only \
  --map-tokens "${MAP_TOKENS}" \
  --yes-always \
  --auto-commits \
  "${readonly_flags[@]}" \
  --message "Review this module strictly for Streamlit UI/UX: control
layout, dashboard state handling, rendering efficiency, label/help-text
clarity, and consistency with the rest of the app's manager/workspace
pattern. Do not change any statistical or bioinformatics logic." \
  "$MODULE_DIR"

# --- Pass 3: Doc-Explainer pass -- EDITS APPLIED, auto-committed ------------
echo ""
echo "--- Pass 3/3: Doc-Explainer pass (editing, auto-commit) ---"
aider \
  --model "ollama_chat/doc-explainer-7b" \
  --subtree-only \
  --map-tokens "${MAP_TOKENS}" \
  --yes-always \
  --auto-commits \
  "${readonly_flags[@]}" \
  --message "Add or improve plain-language docstrings, inline comments, and
Streamlit help= tooltips so a non-expert user or new collaborator can
understand what each function/control does and why. Do not change any
logic at all -- documentation and comments only." \
  "$MODULE_DIR"

echo ""
echo "=== All passes complete for ${MODULE_DIR} ==="
echo "Review the git log/diff before pushing:"
echo "  git log --oneline -5"
echo "  git diff HEAD~2"
echo "Bio-Coder's non-applied review notes are in: ${REVIEW_LOG}"
