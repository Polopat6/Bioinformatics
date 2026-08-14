"""
sra_manager.py

Look up and download publicly available sequencing runs directly from
NCBI's Sequence Read Archive (SRA), given any valid NCBI accession —
a run (SRR/ERR/DRR...), an experiment, a study (SRP...), or a whole
BioProject (PRJNA...).

Motivation: users testing this tool with public datasets can easily
grab data that *looks* like RNA-seq but is actually a different assay
type (e.g. whole-genome sequencing, ChIP-seq) — this is indistinguishable
from real RNA-seq FASTQ files just by looking at them, and only shows up
later as an unexpectedly low Salmon/STAR mapping rate. This module
surfaces each run's official "Library Strategy" metadata field up front,
so a mismatch (e.g. "WGS" when RNA-Seq was expected) is flagged clearly
before any time is spent downloading or quantifying the wrong data.

This module also extracts each sample's SAMPLE_ATTRIBUTES (the
depositor-provided characteristics fields — e.g. "treatment", "genotype",
"media", "strain" — the same fields shown on a GEO record's sample
page), so a usable metadata table can be built automatically for samples
fetched from SRA, rather than requiring the user to manually type it in.

--- Known limitation ---
Not every SRA record has SAMPLE_ATTRIBUTES populated, even when the
underlying study was originally deposited via GEO (which usually has
rich per-sample "characteristics" fields). GEO's characteristics don't
always propagate into SRA's own XML schema — this is a known gap in how
NCBI links the two databases, not a bug in this code. When this happens,
we fall back to:
  1. The SAMPLE-level TITLE, if present
  2. The EXPERIMENT-level TITLE, if the SAMPLE title is missing (GEO-
     derived records sometimes place the original GEO sample title here
     instead)
  3. The GEO accession (GSM number), if present, so the user can look up
     the original GEO record's characteristics table directly — this is
     the ground-truth source when SRA's own docsum is genuinely bare.

--- API note ---
This uses NCBI's E-utilities (esearch + efetch against the "sra"
database), which is the current, actively maintained way to query SRA
metadata. An older endpoint
(trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?...rettype=runinfo) that many
older tools/tutorials reference has been deprecated by NCBI and now
returns "400 Bad Request" for any query — E-utilities is the correct
modern replacement and is what this module uses.

--- Multiple accessions at once ---
esearch's term syntax does not reliably support treating a raw,
comma-separated list of accession numbers as an OR-query the way one
might expect (e.g. "SRR1,SRR2" is not guaranteed to behave as "SRR1 OR
SRR2") -- so looking up a BATCH of individual accessions (as opposed to
one study/BioProject accession that itself expands to many runs) is
handled by lookup_multiple_accessions() below: each individual term is
esearch'd separately (cheap, and NCBI's default rate limit of ~3
requests/second comfortably supports a batch of even a few dozen
accessions), but every term's resulting UIDs are combined into a single
efetch call, so a batch of many accessions still only costs ONE
metadata-fetch request rather than one per accession.

Kept as its own module for the same reason as reference_manager.py and
quantification_manager.py: this involves network requests, XML parsing,
and subprocess execution that is logically distinct from the Streamlit
UI/workflow code.
"""

import gzip
import os
import re
import shutil
import ssl
import subprocess
import time
import urllib.parse
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed

# Same SSL context workaround as reference_manager.py — see that module
# for the full explanation of why this is necessary on some systems
# (notably macOS with Python installed from python.org).
try:
    import certifi
    _SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    _SSL_CONTEXT = ssl.create_default_context()

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# Identifying this tool to NCBI's API, per their usage guidelines. No API
# key is used, so requests are subject to NCBI's default (unauthenticated)
# rate limit of ~3 requests/second, which is more than sufficient here.
TOOL_NAME = "spatial-tissue-atlas-rnaseq-portal"

# NCBI's E-utilities rate limit for unauthenticated requests (no API
# key) is ~3 requests/second -- but this is NOT simply a long-run
# average a client can freely burst against; firing several requests
# back-to-back in a tight loop (e.g. one esearch call per accession in
# a batch lookup -- see lookup_multiple_accessions below) can trigger a
# "429 Too Many Requests" response even when comfortably under the
# long-run average, since NCBI's rate limiter appears to also consider
# short-term burstiness. This was hit directly during real testing with
# a batch of 8 accessions looked up in quick succession. Two measures
# address this together:
#   1. A minimum delay is enforced between consecutive E-utilities
#      requests (see _rate_limit_wait below) -- proactively spacing
#      requests out rather than firing them as fast as Python can loop,
#      which is what caused the burst in the first place.
#   2. Any request that still gets a 429 despite the proactive spacing
#      is automatically retried with exponential backoff (see
#      _fetch_url_with_retry below) -- a 429 is an explicitly transient,
#      recoverable condition per NCBI's own semantics (unlike a 400/404,
#      which indicate a genuinely malformed/nonexistent request and
#      should NOT be retried), so silently giving up on the very first
#      429 -- as the previous version of this module did -- needlessly
#      fails a request that would very likely succeed moments later.
_MIN_SECONDS_BETWEEN_REQUESTS = 0.4  # a bit under 3/second, with headroom
_MAX_429_RETRIES = 5
_last_request_time = [0.0]  # mutable single-element list -- simple, GIL-safe module-level "shared" state


def _rate_limit_wait():
    """
    Block, if necessary, so that at least _MIN_SECONDS_BETWEEN_REQUESTS
    has elapsed since the last E-utilities request this process made --
    a simple, proactive throttle to avoid bursting requests at NCBI
    faster than its rate limit tolerates (see the rate-limit comment
    above for why this matters even when the long-run average request
    rate is comfortably under NCBI's stated 3/second limit).
    """
    elapsed = time.time() - _last_request_time[0]
    if elapsed < _MIN_SECONDS_BETWEEN_REQUESTS:
        time.sleep(_MIN_SECONDS_BETWEEN_REQUESTS - elapsed)
    _last_request_time[0] = time.time()


# fasterq-dump's default --disk-limit-tmp is conservative and can be
# exceeded by larger (e.g. human) RNA-seq runs, which need substantial
# temporary scratch space during extraction/decompression — this was
# never an issue with small bacterial genomes (e.g. E. coli) but shows
# up immediately with human datasets like "airway". Set generously high
# (500GB) so it practically never blocks a legitimate extraction; the
# check still exists as a safety net against genuinely pathological
# cases, just no longer tripped by normal large-genome RNA-seq runs.
FASTERQ_DISK_LIMIT_BYTES = 500_000_000_000


def _fetch_url(url, timeout=30):
    """
    Fetch a URL's raw bytes, using the certifi-backed SSL context.

    Every call is preceded by a proactive rate-limit wait (see
    _rate_limit_wait), and a "429 Too Many Requests" response is
    automatically retried with exponential backoff up to
    _MAX_429_RETRIES times before finally propagating the error to the
    caller -- see the module-level rate-limit comment above for the
    full rationale (a 429 is explicitly a transient/recoverable
    condition per NCBI's own semantics, unlike other HTTP errors like
    400/404 which indicate a genuinely malformed/nonexistent request
    and are correctly NOT retried here, propagating immediately as
    before).
    """
    attempt = 0
    while True:
        _rate_limit_wait()
        try:
            with urllib.request.urlopen(url, timeout=timeout, context=_SSL_CONTEXT) as response:
                return response.read()
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < _MAX_429_RETRIES:
                # Exponential backoff: 1s, 2s, 4s, 8s, 16s -- generous
                # enough to reliably clear a transient rate-limit
                # window without making a genuinely failed batch lookup
                # hang for an unreasonably long time.
                backoff_seconds = 2 ** attempt
                time.sleep(backoff_seconds)
                attempt += 1
                continue
            raise


def tools_available():
    """Check whether the SRA Toolkit (prefetch, fasterq-dump) is installed."""
    return shutil.which("prefetch") is not None, shutil.which("fasterq-dump") is not None


def _esearch_sra(term, retmax=200):
    """
    Convert a search term (accession, study, or BioProject) into a list
    of internal NCBI UIDs for the "sra" database, using esearch.

    Returns (uid_list, error_message_or_None).
    """
    params = {
        "db": "sra",
        "term": term,
        "retmode": "json",
        "retmax": str(retmax),
        "tool": TOOL_NAME,
    }
    url = f"{EUTILS_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"

    try:
        raw = _fetch_url(url)
    except (urllib.error.URLError, TimeoutError) as e:
        return None, f"Could not reach NCBI (esearch): {e}"
    except urllib.error.HTTPError as e:
        return None, f"NCBI rejected the search request (HTTP {e.code}): {e.reason}"

    import json
    try:
        data = json.loads(raw)
        return data["esearchresult"].get("idlist", []), None
    except (json.JSONDecodeError, KeyError) as e:
        return None, f"Could not parse NCBI's search response: {e}"


def _efetch_sra_full_xml(uid_list):
    """
    Fetch full XML metadata records for the given list of SRA UIDs using
    efetch (rettype=full, retmode=xml). This returns an
    EXPERIMENT_PACKAGE_SET containing one EXPERIMENT_PACKAGE per
    experiment, each of which may contain one or more sequencing runs.

    Returns (xml_root_or_None, error_message_or_None).
    """
    if not uid_list:
        return None, "No matching records found."

    params = {
        "db": "sra",
        "id": ",".join(uid_list),
        "rettype": "full",
        "retmode": "xml",
        "tool": TOOL_NAME,
    }
    url = f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"

    try:
        raw = _fetch_url(url, timeout=60)
    except (urllib.error.URLError, TimeoutError) as e:
        return None, f"Could not reach NCBI (efetch): {e}"
    except urllib.error.HTTPError as e:
        return None, f"NCBI rejected the fetch request (HTTP {e.code}): {e.reason}"

    try:
        root = ET.fromstring(raw)
        return root, None
    except ET.ParseError as e:
        return None, f"Could not parse NCBI's metadata response: {e}"


def _parse_sample_attributes(sample_el):
    """
    Extract SAMPLE_ATTRIBUTES (TAG/VALUE pairs) from a <SAMPLE> XML
    element. These are the depositor-provided characteristics fields —
    the same ones shown on a GEO sample page (e.g. "treatment",
    "genotype", "media", "strain", "tissue", "cell type") — and are what
    let us auto-build a usable metadata table for samples fetched from
    SRA, rather than requiring manual entry.

    Returns a dict of {tag: value}. Tag names are used as-is from NCBI
    (lowercased and with spaces replaced by underscores for convenience
    as column names later).
    """
    attributes = {}
    for attr in sample_el.findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
        tag_el = attr.find("TAG")
        value_el = attr.find("VALUE")
        if tag_el is not None and tag_el.text:
            tag = tag_el.text.strip().lower().replace(" ", "_")
            value = value_el.text.strip() if value_el is not None and value_el.text else ""
            attributes[tag] = value
    return attributes


def _find_geo_accession(sample_el):
    """
    Look for a GEO accession (GSM number) in a <SAMPLE> element's
    IDENTIFIERS/EXTERNAL_ID entries. GEO-derived SRA records commonly
    carry the original GSM accession here even when SAMPLE_ATTRIBUTES is
    empty — this lets us point the user directly at the original GEO
    record (which almost always has a full "characteristics" table) as
    a fallback when SRA's own metadata is sparse.

    Returns the GSM accession string, or None if not found.
    """
    for ext_id in sample_el.findall("IDENTIFIERS/EXTERNAL_ID"):
        namespace = (ext_id.get("namespace") or "").upper()
        if namespace == "GEO" and ext_id.text:
            return ext_id.text.strip()
    return None


def _parse_experiment_package(pkg):
    """
    Extract the fields we care about from a single <EXPERIMENT_PACKAGE>
    XML element, returning one dict per run it contains (an experiment
    can have more than one associated run, though 1:1 is most common).

    Each row includes the run-level fields (used for the Library
    Strategy check), a nested "attributes" dict of any SAMPLE_ATTRIBUTES
    found, and fallback descriptive fields (sample/experiment title, GEO
    accession) used when SAMPLE_ATTRIBUTES is empty — see the module
    docstring's "Known limitation" section for why these fallbacks
    exist.
    """
    exp = pkg.find("EXPERIMENT")
    design = exp.find("DESIGN/LIBRARY_DESCRIPTOR") if exp is not None else None

    strategy_el = design.find("LIBRARY_STRATEGY") if design is not None else None
    strategy = strategy_el.text if strategy_el is not None and strategy_el.text else "unknown"

    layout_el = design.find("LIBRARY_LAYOUT") if design is not None else None
    if layout_el is not None and layout_el.find("PAIRED") is not None:
        layout = "PAIRED"
    elif layout_el is not None and layout_el.find("SINGLE") is not None:
        layout = "SINGLE"
    else:
        layout = "unknown"

    exp_title = None
    if exp is not None:
        exp_title_el = exp.find("TITLE")
        if exp_title_el is not None and exp_title_el.text:
            exp_title = exp_title_el.text.strip()

    study = pkg.find("STUDY")
    bioproject = None
    sra_study = study.get("accession") if study is not None else None
    if study is not None:
        for ext_id in study.findall("IDENTIFIERS/EXTERNAL_ID"):
            if ext_id.get("namespace") == "BioProject":
                bioproject = ext_id.text

    sample = pkg.find("SAMPLE")
    organism = "unknown"
    sample_title = None
    geo_accession = None
    attributes = {}
    if sample is not None:
        sci_name_el = sample.find("SAMPLE_NAME/SCIENTIFIC_NAME")
        if sci_name_el is not None and sci_name_el.text:
            organism = sci_name_el.text
        title_el = sample.find("TITLE")
        if title_el is not None and title_el.text:
            sample_title = title_el.text.strip()
        attributes = _parse_sample_attributes(sample)
        geo_accession = _find_geo_accession(sample)

    # Fallback chain for a human-readable description when
    # SAMPLE_ATTRIBUTES is empty: prefer the sample-level title, then
    # fall back to the experiment-level title (GEO-derived records
    # sometimes place the original GEO sample title here instead).
    best_title = sample_title or exp_title

    rows = []
    for run in pkg.findall("RUN_SET/RUN"):
        size_bytes = run.get("size")
        size_mb = round(int(size_bytes) / (1024 * 1024), 1) if size_bytes else None

        rows.append({
            "Run": run.get("accession", "unknown"),
            "LibraryStrategy": strategy,
            "LibraryLayout": layout,
            "ScientificName": organism,
            "spots": run.get("total_spots", "—"),
            "size_MB": size_mb if size_mb is not None else "—",
            "SRAStudy": sra_study or "—",
            "BioProject": bioproject or "—",
            "SampleTitle": best_title or "",
            "GEOAccession": geo_accession or "",
            "attributes": attributes,
        })

    return rows


def lookup_accession(term):
    """
    Query NCBI for metadata about the given accession/search term (a run,
    experiment, study, or BioProject accession).

    Returns (success: bool, rows: list[dict] or None, message: str).
    Each row dict includes: Run, LibraryStrategy, LibraryLayout,
    ScientificName, spots, size_MB, SRAStudy, BioProject, SampleTitle,
    GEOAccession, and a nested "attributes" dict of depositor-provided
    sample characteristics (used to auto-build metadata).

    For looking up MULTIPLE individual accessions at once (e.g. a batch
    of SRR run numbers pasted or uploaded together), use
    lookup_multiple_accessions() instead -- see that function's
    docstring and the module docstring's "Multiple accessions at once"
    section for why a single combined term string isn't the right tool
    for that case.
    """
    term = term.strip()
    if not term:
        return False, None, "Please enter an accession or search term."

    uid_list, error = _esearch_sra(term)
    if error:
        return False, None, error
    if not uid_list:
        return False, None, (
            f"No results found for '{term}'. Double-check the accession "
            "is correct (e.g. SRR12345678, PRJNA123456, SRP123456)."
        )

    root, error = _efetch_sra_full_xml(uid_list)
    if error:
        return False, None, error

    rows = []
    for pkg in root.findall("EXPERIMENT_PACKAGE"):
        rows.extend(_parse_experiment_package(pkg))

    if not rows:
        return False, None, (
            f"Found matching records for '{term}', but couldn't extract "
            "run-level details from NCBI's response. This can happen "
            "with unusual/malformed records — try a more specific "
            "accession (e.g. a single SRR run) instead."
        )

    return True, rows, f"Found {len(rows)} run(s)."


def _split_accession_list_text(text):
    """
    Split a raw text blob containing one or more NCBI accessions into a
    clean list of individual accession strings. Accepts accessions
    separated by any combination of commas, whitespace (spaces, tabs),
    semicolons, and newlines -- e.g. "SRR1, SRR2\\nSRR3" and
    "SRR1 SRR2,SRR3" both parse into ["SRR1", "SRR2", "SRR3"], so a user
    can paste a list in essentially any reasonable format (one per line,
    comma-separated, space-separated, or a mix) without needing to know
    which exact format this tool expects.

    Returns a list of non-empty, stripped accession strings, in the
    order they first appeared, with exact duplicates removed
    (order-preserving) -- safe to call with an empty/whitespace-only
    string, returning [] in that case.
    """
    if not text or not text.strip():
        return []
    raw_tokens = re.split(r"[,\s;]+", text.strip())
    seen = set()
    result = []
    for tok in raw_tokens:
        tok = tok.strip()
        if tok and tok not in seen:
            seen.add(tok)
            result.append(tok)
    return result


def parse_accessions_from_file(file_or_path):
    """
    Extract a list of accession strings from an uploaded file -- .txt or
    .csv (read as plain text, then split using the exact same rule as a
    pasted text blob -- see _split_accession_list_text -- so a
    one-accession-per-line file, a comma-separated file, or any mix of
    the two all work without needing real CSV-dialect parsing) -- or
    .xlsx/.xls (every non-empty cell across every row/column is
    collected as a potential accession; deliberately permissive, so
    accessions pasted into a single column, spread across multiple
    columns, or with a stray header/title row all work -- a header cell
    like "Run" or "Accession" just becomes one harmless "accession" that
    later returns zero SRA results, rather than causing a parsing
    failure).

    file_or_path: either a Streamlit UploadedFile object (or any
    file-like object exposing .name and .read()) OR a plain path string
    -- matching the same dual-input convention used by
    bulk_rnaseq_workspace.py's other file inputs (e.g.
    _read_metadata_file), so this function works identically whether the
    accession list came from a browser upload or the server-side file
    browser.

    Returns (accession_list, error_message_or_None). accession_list is
    always a plain list (possibly empty) even when an error occurred, so
    callers can safely iterate it either way.
    """
    filename = (
        file_or_path.name if hasattr(file_or_path, "name") else os.path.basename(file_or_path)
    ).lower()

    try:
        if filename.endswith((".xlsx", ".xls")):
            import pandas as pd
            df = pd.read_excel(file_or_path, header=None)
            tokens = []
            for value in df.values.flatten():
                if value is None:
                    continue
                text = str(value).strip()
                if text and text.lower() != "nan":
                    tokens.extend(_split_accession_list_text(text))
            seen = set()
            result = []
            for tok in tokens:
                if tok not in seen:
                    seen.add(tok)
                    result.append(tok)
            return result, None
        else:
            # .txt, .csv, or anything else not recognized as Excel --
            # read as plain text and split the same way as a pasted
            # text blob. This handles a simple one-accession-per-line
            # .csv perfectly well without needing a real CSV parser,
            # since comma-splitting is already part of
            # _split_accession_list_text.
            if hasattr(file_or_path, "read"):
                raw = file_or_path.read()
                if isinstance(raw, bytes):
                    raw = raw.decode("utf-8", errors="ignore")
            else:
                with open(file_or_path, "r", errors="ignore") as f:
                    raw = f.read()
            return _split_accession_list_text(raw), None
    except Exception as e:
        return [], f"Could not read this file: {e}"


def lookup_multiple_accessions(terms):
    """
    Query NCBI for metadata about MULTIPLE accessions/search terms at
    once -- e.g. a batch of individual SRR run accessions pasted or
    uploaded together, as opposed to lookup_accession's single term
    (which is still the right tool for a single study/BioProject
    accession that itself expands to many runs).

    Each individual term is looked up via esearch separately (esearch's
    query syntax does not reliably support OR-ing a raw comma-separated
    list of accessions together the way one might expect), but every
    term's resulting UIDs are combined into ONE single efetch call --
    so a batch of, say, 8 run accessions still only costs 8 cheap
    esearch requests (well within NCBI's ~3/second unauthenticated rate
    limit for a batch of this size) plus exactly ONE metadata-fetch
    request, rather than 8 of each.

    terms: list of accession/search term strings (e.g. ["SRR1039508",
        "SRR1039509", ...]) -- duplicates and blank entries are safe to
        include; they're deduplicated automatically before any network
        request is made.

    Returns (success: bool, rows: list[dict] or None, message: str,
    not_found: list[str]).
        - rows has the exact same shape as lookup_accession's rows.
        - not_found lists any individual term that returned zero
          esearch results (e.g. a typo'd or invalid accession), so the
          caller can surface EXACTLY which entries in a pasted/uploaded
          batch didn't resolve, rather than only reporting a vague
          aggregate failure/success count.
    """
    cleaned_terms = []
    seen = set()
    for t in terms:
        t = (t or "").strip()
        if t and t not in seen:
            seen.add(t)
            cleaned_terms.append(t)

    if not cleaned_terms:
        return False, None, "No accessions were provided.", []

    all_uids = []
    uid_seen = set()
    not_found = []

    for term in cleaned_terms:
        uid_list, error = _esearch_sra(term)
        if error:
            return False, None, error, []
        if not uid_list:
            not_found.append(term)
            continue
        for uid in uid_list:
            if uid not in uid_seen:
                uid_seen.add(uid)
                all_uids.append(uid)

    if not all_uids:
        return False, None, (
            f"None of the {len(cleaned_terms)} accession(s) provided "
            "returned any results. Double-check they're correct (e.g. "
            "SRR12345678, PRJNA123456, SRP123456)."
        ), not_found

    root, error = _efetch_sra_full_xml(all_uids)
    if error:
        return False, None, error, not_found

    rows = []
    for pkg in root.findall("EXPERIMENT_PACKAGE"):
        rows.extend(_parse_experiment_package(pkg))

    if not rows:
        return False, None, (
            "Found matching records, but couldn't extract run-level "
            "details from NCBI's response."
        ), not_found

    n_resolved = len(cleaned_terms) - len(not_found)
    message = f"Found {len(rows)} run(s) from {n_resolved} of {len(cleaned_terms)} accession(s)."
    return True, rows, message, not_found


def is_rna_seq(row):
    """
    Check whether a run's metadata indicates it's actually RNA-Seq data,
    based on NCBI's own "Library Strategy" annotation for that run.
    """
    return row.get("LibraryStrategy", "").strip().upper() == "RNA-SEQ"


def is_paired_layout(row):
    return row.get("LibraryLayout", "").strip().upper() == "PAIRED"


def build_metadata_dataframe(rows, selected_runs=None):
    """
    Build a metadata table (as a list of dicts, ready to become a
    pandas DataFrame) from SRA lookup rows, using each sample's
    depositor-provided SAMPLE_ATTRIBUTES as columns, with a fallback to
    descriptive title text and the GEO accession when no structured
    attributes are available (see module docstring's "Known
    limitation").

    This lets a user who fetched samples via the SRA lookup feature skip
    manually typing out a metadata CSV — the same characteristics fields
    visible on the original GEO/SRA record (e.g. "treatment", "media",
    "genotype") become columns automatically, with the Run accession as
    the "sample" column (matching the FASTQ filenames fasterq-dump
    produces, e.g. SRR12345678_1.fastq.gz).

    If selected_runs is given (a list/set of Run accessions), only those
    rows are included — otherwise all rows are used.

    Returns a list of dicts (one per sample), suitable for
    pd.DataFrame(...). Column names beyond "sample" depend entirely on
    whatever attributes NCBI/the depositor actually provided, so this
    will vary by dataset — the caller should treat it as a helpful
    starting point that the user can review/edit, not a guaranteed-
    complete metadata table.
    """
    if selected_runs is not None:
        selected_runs = set(selected_runs)
        rows = [r for r in rows if r["Run"] in selected_runs]

    metadata_rows = []
    for row in rows:
        entry = {"sample": row["Run"]}
        if row.get("SampleTitle"):
            entry["sample_title"] = row["SampleTitle"]
        if row.get("GEOAccession"):
            entry["geo_accession"] = row["GEOAccession"]
        entry.update(row.get("attributes", {}))
        metadata_rows.append(entry)

    return metadata_rows


def has_any_descriptive_metadata(rows, selected_runs=None):
    """
    Check whether the given rows have *anything* usable beyond the bare
    Run accession — SAMPLE_ATTRIBUTES, a title, or a GEO accession.
    Used to decide whether to show the "nothing useful found" warning
    versus letting the user proceed with at least a title/GEO link to
    work from.
    """
    if selected_runs is not None:
        selected_runs = set(selected_runs)
        rows = [r for r in rows if r["Run"] in selected_runs]

    for row in rows:
        if row.get("attributes") or row.get("SampleTitle") or row.get("GEOAccession"):
            return True
    return False


def geo_url(geo_accession):
    """Build a direct link to a GEO sample's page, given its GSM accession."""
    return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_accession}"


def _gzip_and_remove(src_path):
    """Compress src_path to src_path + '.gz', then delete the original."""
    gz_path = src_path + ".gz"
    with open(src_path, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(src_path)
    return gz_path


def download_sra_run(accession, dest_dir, progress_callback=None, threads=4):
    """
    Download a single SRA run accession's FASTQ file(s) using the SRA
    Toolkit (prefetch + fasterq-dump), then gzip the output and place it
    directly in dest_dir (intended to be a project's fastq_dir, so the
    files are immediately picked up by the existing FASTQ matching logic
    in bulk_rnaseq_workspace.py — no separate handling needed).

    fasterq-dump's default output naming for paired-end data is
    "<accession>_1.fastq" / "<accession>_2.fastq", which already matches
    this tool's existing "_1"/"_2" pairing detection — no renaming
    required. Single-end data is named "<accession>.fastq".

    threads: fasterq-dump's own "-e"/"--threads" flag, controlling how
    many threads it uses internally while extracting/decompressing a
    single accession. Note this is a distinct axis of parallelism from
    download_sra_runs_parallel's max_workers, which controls how many
    *different accessions* download at once — the two multiply together
    (e.g. 3 concurrent accessions x 4 threads each = up to 12 threads
    active simultaneously), so it's worth lowering one if raising the
    other on a smaller machine.

    Implementation note: every path used here is converted to an
    *absolute* path up front, and no subprocess call uses a changed
    working directory (cwd=). An earlier version of this function passed
    a relative path to fasterq-dump while also running it with a changed
    cwd — that combination causes the relative path to be re-resolved
    against the new working directory, silently pointing at the wrong
    location and producing a confusing
    "failed to resolve accession ... 404" error even though prefetch had
    already succeeded. Using absolute paths throughout avoids this
    category of bug entirely.

    Implementation note (disk limit): fasterq-dump has a built-in safety
    check (--disk-limit-tmp) on how much temporary scratch space it's
    allowed to use during extraction/decompression, and its default is
    conservative enough that larger runs (e.g. human RNA-seq, such as
    the "airway" dataset) can exceed it even when the machine actually
    has plenty of free disk space — failing with "disk-limit exceeded!"
    before extraction even starts. This was never triggered by small
    bacterial genome runs (e.g. E. coli) in earlier testing, only
    appearing once larger eukaryotic datasets were used. We pass an
    explicit, generously high limit (FASTERQ_DISK_LIMIT_BYTES, 500GB) to
    avoid this for any realistic RNA-seq run size, while still leaving
    some safety check in place against truly pathological cases.

    Returns (success: bool, log: str).
    """
    dest_dir = os.path.abspath(dest_dir)
    os.makedirs(dest_dir, exist_ok=True)

    if progress_callback:
        progress_callback(f"Downloading {accession} from NCBI (prefetch)...")

    prefetch_cmd = ["prefetch", accession, "-O", dest_dir]
    try:
        result = subprocess.run(
            prefetch_cmd, capture_output=True, text=True, check=True, timeout=3600
        )
        log = result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, f"prefetch failed: {(e.stdout or '') + (e.stderr or '')}"
    except subprocess.TimeoutExpired:
        return False, f"prefetch timed out after 1 hour on '{accession}'."

    if progress_callback:
        progress_callback(f"Extracting FASTQ from {accession} (fasterq-dump)...")

    # prefetch downloads into dest_dir/<accession>/<accession>.sra by
    # default. We point fasterq-dump directly at that specific .sra file
    # (rather than just the containing directory or the bare accession
    # name) since this is the most unambiguous, reliable input — it
    # doesn't depend on fasterq-dump's own accession-resolution/caching
    # behavior working correctly in this environment.
    sra_subdir = os.path.join(dest_dir, accession)
    sra_file = os.path.join(sra_subdir, f"{accession}.sra")

    if os.path.exists(sra_file):
        fasterq_input = sra_file
    elif os.path.isdir(sra_subdir):
        fasterq_input = sra_subdir
    else:
        # Fall back to the bare accession name, letting fasterq-dump's
        # own SRA cache resolution find it — this matches prefetch's
        # default caching behavior if it placed the file somewhere other
        # than the expected dest_dir/<accession>/<accession>.sra layout.
        fasterq_input = accession

    fasterq_cmd = [
        "fasterq-dump", fasterq_input,
        "--split-files",
        "-O", dest_dir,
        "-e", str(threads),
        "--disk-limit-tmp", str(FASTERQ_DISK_LIMIT_BYTES),
    ]
    try:
        result = subprocess.run(
            fasterq_cmd, capture_output=True, text=True, check=True, timeout=3600,
        )
        log += result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, log + f"\nfasterq-dump failed: {(e.stdout or '') + (e.stderr or '')}"
    except subprocess.TimeoutExpired:
        return False, log + f"\nfasterq-dump timed out after 1 hour on '{accession}'."

    if progress_callback:
        progress_callback(f"Compressing FASTQ output for {accession}...")

    produced_files = []
    for candidate in (f"{accession}_1.fastq", f"{accession}_2.fastq", f"{accession}.fastq"):
        candidate_path = os.path.join(dest_dir, candidate)
        if os.path.exists(candidate_path):
            gz_path = _gzip_and_remove(candidate_path)
            produced_files.append(os.path.basename(gz_path))

    # Clean up the intermediate .sra download to save disk space, now
    # that we have the extracted (and compressed) FASTQ files.
    if os.path.isdir(sra_subdir):
        shutil.rmtree(sra_subdir, ignore_errors=True)

    if not produced_files:
        return False, log + "\nNo FASTQ output files were found after extraction."

    return True, f"Downloaded and compressed: {', '.join(produced_files)}"


def download_sra_runs_parallel(accessions, dest_dir, max_workers=3, on_run_complete=None, threads_per_run=4):
    """
    Download multiple SRA run accessions concurrently, using a thread
    pool. This is a meaningful speedup over downloading one accession at
    a time, since prefetch/fasterq-dump are almost entirely I/O-bound
    (waiting on network transfer from NCBI, not doing heavy local
    computation) — running several downloads at once lets one
    accession's network wait time overlap with another's, rather than
    the whole batch being fully serial.

    max_workers defaults to 3 as a reasonably safe concurrency level:
    high enough to meaningfully overlap network waits, but conservative
    enough to avoid overwhelming either the local network connection or
    NCBI's servers (NCBI does not publish a hard concurrent-connection
    limit for the SRA Toolkit specifically, but being a considerate,
    moderate default is good practice for a shared public resource).
    Increase this if your connection is fast and you're downloading many
    accessions; decrease it if downloads seem to be failing or timing
    out under load.

    threads_per_run is passed straight through to each individual
    download_sra_run() call's own "threads" parameter (fasterq-dump's
    "-e" flag). Note this multiplies with max_workers for total
    simultaneous thread usage (e.g. 3 concurrent accessions x 4 threads
    each = up to 12 threads active at once) — on smaller machines,
    lowering one when raising the other helps avoid over-committing the
    system.

    Note: since each download now allows fasterq-dump to use up to
    FASTERQ_DISK_LIMIT_BYTES of temporary scratch space, running several
    large (e.g. human) downloads in parallel multiplies actual peak disk
    usage accordingly — if you're downloading multiple large runs at
    once and are concerned about disk space, consider lowering
    max_workers for that batch.

    IMPORTANT (Streamlit-specific): this function itself does NOT touch
    any Streamlit UI elements (st.*) from within the worker threads —
    only individual accessions' plain success/log strings are returned.
    Streamlit is not thread-safe for UI updates from background threads,
    so any progress display must be driven from the *calling* code in
    the main thread (e.g. by iterating results as they complete and
    updating a progress bar there), not from inside this function or the
    on_run_complete callback if you choose to pass one that touches the
    UI directly.

    on_run_complete, if given, is called as on_run_complete(accession,
    success, message) from the *main* thread (not a worker thread) each
    time one accession finishes — this is safe to use for Streamlit UI
    updates like advancing a progress bar, since as_completed() below
    yields control back to the main thread for each callback invocation.

    Returns a list of (accession, success, message) tuples, in
    completion order (not necessarily the same order as the input list,
    since faster downloads finish first).
    """
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_accession = {
            executor.submit(download_sra_run, accession, dest_dir, None, threads_per_run): accession
            for accession in accessions
        }

        for future in as_completed(future_to_accession):
            accession = future_to_accession[future]
            try:
                success, message = future.result()
            except Exception as e:
                # Catch-all: a worker thread raising an unexpected
                # exception (rather than returning a (success, message)
                # tuple as designed) shouldn't crash the whole batch —
                # report it as a failure for that one accession instead.
                success, message = False, f"Unexpected error: {e}"

            results.append((accession, success, message))
            if on_run_complete:
                on_run_complete(accession, success, message)

    return results
