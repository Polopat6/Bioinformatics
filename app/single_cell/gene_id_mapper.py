"""
single_cell/gene_id_mapper.py

Thin re-export of the real gene_id_mapper.py (repo/app/gene_id_mapper.py)
-- exists only so single-cell modules (which resolve imports against
single_cell/'s own sys.path entry, not app/'s) can `import gene_id_mapper
as gim` and get the SAME shared implementation the bulk pipeline uses.

IMPORTANT: this canNOT be done with a plain `from gene_id_mapper import *`
-- since this file is itself named gene_id_mapper.py, Python's import
system registers "gene_id_mapper" in sys.modules as the CURRENTLY-LOADING
module (this file) before this file's own code even finishes running.
A same-named import then just refers back to this half-initialized
module instead of reaching the real one in app/, silently producing a
module with none of the real functions defined. Loading the real file
via importlib under a DIFFERENT internal name avoids this self-reference
entirely.

NOTE (2026-09-15): this file was found to have reverted to a plain
duplicate at some point after initially being fixed -- if gene ID
detection/conversion ever appears to silently regress again on the
single-cell side despite app/gene_id_mapper.py being correct, check
here FIRST before assuming the bug is back in the real module.
"""
import importlib.util
import os
import sys

_this_dir = os.path.dirname(os.path.abspath(__file__))
_real_path = os.path.join(_this_dir, "..", "gene_id_mapper.py")

_spec = importlib.util.spec_from_file_location("_gene_id_mapper_real", _real_path)
_real_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_real_module)

_current_module = sys.modules[__name__]
for _name in dir(_real_module):
    if not _name.startswith("__"):
        setattr(_current_module, _name, getattr(_real_module, _name))
