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


def _read_ids(fastq_path):
    opener = gzip.open if fastq_path.endswith(".gz") else open
    with opener(fastq_path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                read_id = line.strip().lstrip("@").split()[0]
                if read_id.endswith(("/1", "/2")):
                    read_id = read_id[:-2]
                yield read_id


def resync_r1_to_r2(r1_in, r2_trimmed, r1_out):
    "Filter r1_in down to only reads whose ID also survives in r2_trimmed, preserving original order. Returns the number of reads written."
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
    "Full single-sample trim: fastp on R2 (single-end) -> resync R1 down to R2's survivors."
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
    n_written = resync_r1_to_r2(r1_in, output_paths["r2_out"], output_paths["r1_out"])
    return True, f"Trimmed successfully -- {n_written} read pairs retained."
