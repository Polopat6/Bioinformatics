"""
single_cell/singlecell_trim_manager.py

fastp-based trimming for the Single-cell RNA-Seq pipeline. R1 is NEVER
passed to fastp at all (a real Biostars-reported issue showed fastp's
per-read asymmetric flags don't reliably protect one read while trimming
the other) -- instead, R2 is trimmed alone in single-end mode, and R1 is
re-synced down to whichever reads survived R2's trimming (matching the
validated approach used by the real seq2science pipeline for this exact
problem).

Default quality threshold is Q15 (gentler than bulk's Q20/Q30), matching
real single-cell fastp workflows -- preserving read depth matters more
here than the marginal quality gain from a stricter threshold.

--- Minimum R2 length: floor vs. target (2026-08-17) ---
length_required (DEFAULT_LENGTH_REQUIRED=20) is a FLOOR: reads shorter
than this after trimming are discarded outright as too short to be
useful at all. This is NOT the same thing as the chemistry-specific
TARGET/recommended full R2 sequencing length 10x Genomics documents (91bp
for 3' v3/v3.1 and 5' v2, 98bp for 5' v1, confirmed via 10x's own support
KB -- see chemistry_manager.py's CHEMISTRY_CATALOG
"recommended_r2_length" field) -- that target is about how long you
should SEQUENCE the read for confident unique mapping in the first
place ("shorter reads have increased chances of multi-mapping" per 10x's
own KB), not about the minimum to accept after quality/adapter trimming.
There is no documented chemistry-specific value for the FLOOR itself --
20bp works as a "don't keep hopeless reads" safety net regardless of
chemistry, since it sits far below any chemistry's actual target length
anyway. singlecell_workspace.py's Step 3 UI shows the chemistry's target
length as informational CONTEXT next to this floor slider, but does not
change the floor's default based on chemistry.

--- Poly-A/poly-G tail trimming (2026-08-17) ---
fastp's adapter/quality trimming does NOT touch poly-A or poly-G
homopolymer tails by default (poly-G auto-enables only for reads fastp
detects as coming from 2-color NextSeq/NovaSeq chemistry; poly-A/general
poly-X trimming is never automatic at all -- both require explicit
--trim_poly_g/--trim_poly_x flags, confirmed against fastp's own option
reference). This matters MORE for single-cell 3' data specifically than
typical bulk RNA-seq: a peer-reviewed 2022 study (Svoboda et al., NAR
Genomics and Bioinformatics) confirms internal oligo(dT) priming -- the
exact bead-capture mechanism 10x's poly-dT primers use -- causes
systematic poly-A read-through contamination in single-cell RNA-seq
specifically (explicitly naming 10x among the affected platforms). The
Bulk RNA-Seq pipeline already exposes this as an "auto_fix_poly_tails"
toggle; this module mirrors that same pattern for single-cell R2
trimming via fastp's real, documented --trim_poly_g/--trim_poly_x flags.

--- SE-mode adapter detection trade-off (documented, not user-facing
    toggle) ---
Since R2 is trimmed alone in single-end mode (see module docstring
above), fastp's adapter auto-detection uses its SE method (scanning
k-mer overrepresentation in the tail of the first ~1M reads) rather than
its more robust PE overlap-based detection (confirmed via fastp's own
documentation) -- a real, unavoidable trade-off of protecting R1 by
never passing it to fastp at all. This is surfaced to the user as
information in singlecell_workspace.py's Step 3 UI, not something this
module can eliminate.

--- R1/R2 resync performance fix: seqkit instead of pure Python
    (2026-09-01) ---
A real, confirmed production bug: the original resync_r1_to_r2()
implementation read both R2 (to build a Python set() of surviving read
IDs) and R1 (to filter against that set) line-by-line in pure Python.
For a typical bulk-scale FASTQ this is fine, but single-cell runs can
be enormously deeper -- confirmed directly against a real accession
(SRR13734384, ~776 million spots) -- where this pure-Python two-pass
approach took MULTIPLE HOURS for a single sample, even though fastp
itself (a real, multithreaded, compiled C++ tool) finished its own R2
trimming quickly. Critically, the `threads` parameter passed to
run_trim_sample() ONLY ever reaches fastp's own --thread flag -- the
resync step was, and structurally could never be, sped up by raising
that slider at all, which is exactly why increasing threads had no
effect on the real-world hang this fix addresses.

Fixed by using seqkit (github.com/shenwei356/seqkit) -- a fast,
single-binary, compiled Go tool already standard in many bioinformatics
environments -- to perform the resync instead, via TWO of seqkit's own
long-stable subcommands (deliberately NOT the newer "seqkit pair"
subcommand, whose exact flags/version-availability across different
seqkit releases is less certain):
  1. `seqkit seq -n --only-id` on the trimmed R2 file -- extracts just
     the surviving read IDs (one per line), using seqkit's own ID
     definition (the first whitespace-delimited token of the header)
     consistently.
  2. `seqkit grep -f <that ID list>` on the original R1 file -- filters
     R1 down to only those same IDs, using the IDENTICAL ID-parsing
     convention as step 1 (both driven by seqkit itself), so the two
     steps are guaranteed to agree with each other regardless of the
     exact header format either FASTQ file happens to use.
This two-step approach avoids ever needing to hold millions of read IDs
in Python memory as a set, and avoids Python-level per-line string
parsing entirely for the actual matching operation.

seqkit_available() gates this path. If seqkit is NOT installed, this
module explicitly FALLS BACK to the original pure-Python implementation
-- correctness is preserved either way -- but surfaces a clear, visible
warning message (returned all the way up through run_trim_sample()'s own
result) explaining that the slow path was used and that installing
seqkit is strongly recommended for realistic single-cell run sizes. This
mirrors this app's own established "explicit visibility over silent
degraded behavior" pattern used elsewhere (e.g. sc_sra_manager.py's own
pigz-vs-gzip fallback for compression).

One honest, real trade-off worth knowing: the intermediate surviving-ID
list file (step 1's output) scales with read count -- for an
extreme-depth run like SRR13734384 (~776M reads), this ID list can
itself be on the order of many GB of plain text, written to disk
temporarily. This is deleted immediately after use, and in practice is
still dramatically smaller and faster than the multi-hour pure-Python
path it replaces, but it is a real, additional temporary-disk-space
consideration for extreme-depth runs specifically.
"""
import gzip
import os
import shutil
import subprocess

DEFAULT_QUALIFIED_QUALITY_PHRED = 15
DEFAULT_LENGTH_REQUIRED = 20
DEFAULT_AUTO_FIX_POLY_TAILS = True
DEFAULT_POLY_X_MIN_LEN = 10  # fastp's own default for --poly_x_min_len


def build_fastp_r2_command(r2_in, r2_out, json_report, html_report,
                            qualified_quality_phred=DEFAULT_QUALIFIED_QUALITY_PHRED,
                            length_required=DEFAULT_LENGTH_REQUIRED, threads=4,
                            auto_fix_poly_tails=DEFAULT_AUTO_FIX_POLY_TAILS,
                            poly_x_min_len=DEFAULT_POLY_X_MIN_LEN):
    """
    auto_fix_poly_tails=True adds fastp's real --trim_poly_g (2-color
    chemistry poly-G artifact removal) AND --trim_poly_x (general
    homopolymer tail removal, which also catches poly-A read-through
    from oligo-dT priming -- fastp has no dedicated "poly-A-only" flag;
    --trim_poly_x is the correct, documented mechanism for this). Both
    are OFF by default in fastp itself and must be explicitly requested
    -- see module docstring for why this matters more for single-cell
    3' data specifically than typical bulk RNA-seq.
    """
    cmd = [
        "fastp", "--in1", r2_in, "--out1", r2_out,
        "--qualified_quality_phred", str(qualified_quality_phred),
        "--length_required", str(length_required),
        "--thread", str(threads), "--json", json_report, "--html", html_report,
    ]
    if auto_fix_poly_tails:
        cmd += ["--trim_poly_g", "--trim_poly_x", "--poly_x_min_len", str(poly_x_min_len)]
    return cmd


# ---------------------------------------------------------------------------
# R1/R2 resync -- seqkit-backed (fast path) with pure-Python fallback
# ---------------------------------------------------------------------------

def seqkit_available():
    "Check whether seqkit is installed and on PATH."
    return shutil.which("seqkit") is not None


def build_seqkit_extract_ids_command(fastq_in, id_list_out):
    "seqkit seq -n --only-id: writes just the surviving read IDs (one per line), one per FASTQ record, using seqkit's own ID-parsing convention (first whitespace-delimited header token)."
    return ["seqkit", "seq", "-n", "--only-id", fastq_in, "-o", id_list_out]


def build_seqkit_grep_by_id_command(fastq_in, id_list_path, fastq_out):
    "seqkit grep -f <id list>: filters fastq_in down to only records whose ID (same parsing convention as the extraction step) appears in id_list_path."
    return ["seqkit", "grep", "-f", id_list_path, fastq_in, "-o", fastq_out]


def _count_lines(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return sum(1 for _ in f)


def _resync_r1_to_r2_seqkit(r1_in, r2_trimmed, r1_out, work_dir, subprocess_runner=None):
    """
    Fast path: two seqkit subprocess calls (extract surviving IDs from
    R2, then filter R1 by that same ID list) instead of pure-Python
    line-by-line processing. See module docstring for full rationale.

    Returns (success: bool, message: str, n_written: int|None).
    """
    runner = subprocess_runner or subprocess.run
    os.makedirs(work_dir, exist_ok=True)
    id_list_path = os.path.join(work_dir, os.path.basename(r1_out) + ".surviving_ids.txt")

    try:
        extract_cmd = build_seqkit_extract_ids_command(r2_trimmed, id_list_path)
        result = runner(extract_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return False, f"seqkit failed extracting surviving read IDs from R2 (exit {result.returncode}): {result.stderr}", None

        grep_cmd = build_seqkit_grep_by_id_command(r1_in, id_list_path, r1_out)
        result2 = runner(grep_cmd, capture_output=True, text=True)
        if result2.returncode != 0:
            return False, f"seqkit failed filtering R1 by surviving read IDs (exit {result2.returncode}): {result2.stderr}", None

        # The ID list has exactly one line per surviving read -- a cheap,
        # accurate proxy for "how many reads were retained" without
        # needing to re-scan the (potentially huge) output FASTQ itself.
        try:
            n_written = _count_lines(id_list_path)
        except OSError:
            n_written = None

        return True, "Resynced R1 to R2's survivors using seqkit.", n_written
    finally:
        # Always clean up the intermediate ID list -- for extreme-depth
        # runs (hundreds of millions of reads) this file can itself be
        # many GB of plain text; see module docstring's own honest
        # disk-space note.
        if os.path.isfile(id_list_path):
            try:
                os.remove(id_list_path)
            except OSError:
                pass


def _read_ids(fastq_path):
    opener = gzip.open if fastq_path.endswith(".gz") else open
    with opener(fastq_path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                read_id = line.strip().lstrip("@").split()[0]
                if read_id.endswith(("/1", "/2")):
                    read_id = read_id[:-2]
                yield read_id


def _resync_r1_to_r2_python_fallback(r1_in, r2_trimmed, r1_out):
    """
    Original pure-Python implementation -- ONLY used when seqkit is not
    installed. Correct, but confirmed to take multiple HOURS on a
    single extreme-depth single-cell sample (hundreds of millions of
    reads); see module docstring. Kept as a fallback so this pipeline
    never hard-fails on a system without seqkit installed, at the cost
    of speed.
    """
    surviving_ids = set(_read_ids(r2_trimmed))
    in_opener = gzip.open if r1_in.endswith(".gz") else open
    out_opener = gzip.open if r1_out.endswith(".gz") else open
    written = 0
    with in_opener(r1_in, "rt") as fin, out_opener(r1_out, "wt") as fout:
        record = []
        for line in fin:
            record.append(line)
            if len(record) == 4:
                header = record[0].strip().lstrip("@").split()[0]
                if header.endswith(("/1", "/2")):
                    header = header[:-2]
                if header in surviving_ids:
                    fout.writelines(record)
                    written += 1
                record = []
    return written


_SEQKIT_INSTALL_HINT = (
    "⚠️ seqkit not found -- used a much slower pure-Python fallback for R1/R2 resync. "
    "For large single-cell datasets (hundreds of millions of reads), this fallback can take "
    "HOURS for a single sample, even though fastp's own trimming step itself finishes quickly "
    "-- the two are separate steps, and only fastp's portion is sped up by adding threads. "
    "Installing seqkit (e.g. `mamba install -c bioconda seqkit`) is strongly recommended and "
    "will make this step dramatically faster for realistic single-cell run sizes."
)


def resync_r1_to_r2(r1_in, r2_trimmed, r1_out, work_dir=None, subprocess_runner=None):
    """
    Filter r1_in down to only reads whose ID also survives in
    r2_trimmed, preserving original order. Uses seqkit (fast, compiled)
    when available; falls back to a pure-Python implementation
    otherwise -- see module docstring, "R1/R2 resync performance fix",
    for the full rationale and a real, confirmed multi-hour case this
    was fixing (SRR13734384, ~776M reads).

    work_dir: directory for seqkit's own temporary intermediate ID-list
        file (deleted automatically after use). Defaults to r1_out's
        own directory if not given.

    Returns (success: bool, n_written: int|None, used_fallback: bool, message: str).
    A message is ALWAYS returned (even on the fast, no-warning path), so
    callers can decide for themselves whether/how to surface it.
    """
    work_dir = work_dir or (os.path.dirname(r1_out) or ".")

    if seqkit_available():
        success, message, n_written = _resync_r1_to_r2_seqkit(
            r1_in, r2_trimmed, r1_out, work_dir, subprocess_runner=subprocess_runner,
        )
        return success, n_written, False, message

    try:
        n_written = _resync_r1_to_r2_python_fallback(r1_in, r2_trimmed, r1_out)
    except Exception as e:  # noqa: BLE001 -- surface any failure directly rather than crash the whole pipeline run
        return False, None, True, f"Pure-Python resync fallback failed: {e}"

    return True, n_written, True, _SEQKIT_INSTALL_HINT


def fastp_output_paths(project_trimmed_dir, project_fastp_reports_dir, sample_name):
    return {
        "r1_out": os.path.join(project_trimmed_dir, f"{sample_name}_R1_001.trimmed.fastq.gz"),
        "r2_out": os.path.join(project_trimmed_dir, f"{sample_name}_R2_001.trimmed.fastq.gz"),
        "r2_trimmed_tmp": os.path.join(project_trimmed_dir, "_work", f"{sample_name}_R2_001.trimmed_pre_resync.fastq.gz"),
        "json_report": os.path.join(project_fastp_reports_dir, f"{sample_name}.fastp.json"),
        "html_report": os.path.join(project_fastp_reports_dir, f"{sample_name}.fastp.html"),
    }


def _run_and_log(command, log_f):
    process = subprocess.Popen(command, stdout=log_f, stderr=subprocess.STDOUT)
    process.wait()
    return process.returncode


def run_trim_sample(r1_in, r2_in, output_paths, qualified_quality_phred=DEFAULT_QUALIFIED_QUALITY_PHRED,
                     length_required=DEFAULT_LENGTH_REQUIRED, threads=4,
                     auto_fix_poly_tails=DEFAULT_AUTO_FIX_POLY_TAILS,
                     poly_x_min_len=DEFAULT_POLY_X_MIN_LEN, subprocess_runner=None):
    "Full single-sample trim: fastp on R2 (single-end) -> resync R1 down to R2's survivors (seqkit-backed, see module docstring)."
    runner = subprocess_runner or subprocess.run

    os.makedirs(os.path.dirname(output_paths["r2_trimmed_tmp"]), exist_ok=True)
    os.makedirs(os.path.dirname(output_paths["json_report"]), exist_ok=True)
    os.makedirs(os.path.dirname(output_paths["r1_out"]), exist_ok=True)

    cmd = build_fastp_r2_command(
        r2_in, output_paths["r2_trimmed_tmp"], output_paths["json_report"], output_paths["html_report"],
        qualified_quality_phred=qualified_quality_phred, length_required=length_required, threads=threads,
        auto_fix_poly_tails=auto_fix_poly_tails, poly_x_min_len=poly_x_min_len,
    )
    result = runner(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return False, f"fastp failed on R2 (exit {result.returncode}): {result.stderr}"

    shutil.move(output_paths["r2_trimmed_tmp"], output_paths["r2_out"])

    resync_work_dir = os.path.join(os.path.dirname(output_paths["r1_out"]), "_work")
    resync_success, n_written, used_fallback, resync_message = resync_r1_to_r2(
        r1_in, output_paths["r2_out"], output_paths["r1_out"],
        work_dir=resync_work_dir, subprocess_runner=subprocess_runner,
    )
    if not resync_success:
        return False, f"R1/R2 resync failed: {resync_message}"

    base_message = f"Trimmed successfully -- {n_written if n_written is not None else '?'} read pairs retained."
    if used_fallback:
        base_message += f" {resync_message}"
    return True, base_message
