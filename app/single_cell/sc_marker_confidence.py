"""
single_cell/sc_marker_confidence.py

Solves a real, confirmed usability problem in Step 8e: assigning final
cell-type labels required scrolling back and forth between Step 8b's
mean-marker-score table and Step 8e's own editor to remember which
manual marker panel scored highest for each cluster -- this real gap
directly caused at least one confirmed real mislabeling (a
Megakaryocyte cluster with a clear, strong manual score of 1.78 in
Step 8b ended up saved as "CD4 T cells" in the final labels, only
caught afterward via ground-truth validation against Kang et al.
2018's own published labels).

classify_marker_confidence() below assigns one of 4 plain-language
confidence tiers to a cluster's TOP-scoring manual marker panel, based
on TWO things -- not just "which number is biggest":
  1. Is the top score itself comfortably positive (a real signal), not
     just "least negative of a uniformly weak/negative bunch"?
  2. Is the top score clearly SEPARATED from the runner-up (a genuine
     standout), or only narrowly ahead of a close competitor (a
     genuinely ambiguous case -- e.g. this exact dataset's own
     confirmed CD8 T cell / NK cell transcriptional overlap)?

Intended to be surfaced as a new read-only column in Step 8e's own
cell-type-assignment table, positioned between the CellTypist
suggestion column and the editable final-label column -- so a user can
see the manual-scoring evidence for each cluster WITHOUT needing to
scroll back to Step 8b at all.
"""

HIGH_THRESHOLD = 0.5
POTENTIAL_THRESHOLD = 0.1
MIN_SEPARATION_FOR_HIGH = 0.15


def classify_marker_confidence(sorted_scores):
    """
    sorted_scores: list of (cell_type_name, score) tuples, SORTED
        descending by score.

    Returns (top_cell_type: str or None, confidence: str, reason: str).

    confidence is one of "High", "Potential", "Low", "No signal" -- see
    this module's own docstring for the full rationale behind each
    tier. Returns (None, "No signal", ...) if sorted_scores is empty.
    """
    if not sorted_scores:
        return None, "No signal", "No marker panels have been scored yet."

    top_name, top_score = sorted_scores[0]
    second_score = sorted_scores[1][1] if len(sorted_scores) > 1 else float("-inf")
    separation = top_score - second_score

    if top_score <= 0:
        return top_name, "No signal", (
            f"Top score ({top_name}: {top_score:.2f}) is zero or negative -- no real "
            "positive evidence for any defined panel."
        )

    if top_score >= HIGH_THRESHOLD and separation >= MIN_SEPARATION_FOR_HIGH:
        return top_name, "High", (
            f"{top_name} ({top_score:.2f}) is comfortably positive and clearly ahead of "
            f"the runner-up (by {separation:.2f}) -- a strong, trustworthy standout."
        )

    if top_score >= POTENTIAL_THRESHOLD:
        if separation < MIN_SEPARATION_FOR_HIGH and second_score > 0:
            return top_name, "Potential", (
                f"{top_name} ({top_score:.2f}) is the top score, but only narrowly ahead "
                f"of a close competitor (by {separation:.2f}) -- worth double-checking "
                "against the marker gene table before trusting this alone."
            )
        return top_name, "Potential", (
            f"{top_name} ({top_score:.2f}) is positive but modest -- a plausible lead, "
            "not yet a strong standout."
        )

    return top_name, "Low", (
        f"{top_name} ({top_score:.2f}) is only weakly positive (barely above zero) -- "
        "treat as a weak hint only."
    )


CONFIDENCE_ICONS = {
    "High": "🟢 High",
    "Potential": "🟡 Potential",
    "Low": "🟠 Low",
    "No signal": "⚪ No signal",
}
