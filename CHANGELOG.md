# Changelog

## Unreleased

### Bug Fixes

- **`MolToXYZFile` wrote the same conformer to every `.xyz` file** (`base.py`): Missing
  `confId=conformer_id` argument meant all 50 conformer files were identical, making
  downstream xTB ranking meaningless. Fixed by passing the correct conformer ID.

- **`Monomer` crashed on construction for RER reactions** (`base.py`): `Initiator` was
  always instantiated in `Monomer.__init__`, even for RER where `initiator = ""`. Calling
  `get_mol("")` raised `PolymerBuildError` immediately. `Initiator` is now only created
  when `reaction.type == "ROR"`.

- **`Exception.__init__` packed args into a tuple** (`exceptions.py`): `super().__init__(args)`
  passed the whole `*args` tuple as a single argument, so `str(PolymerBuildError("message"))`
  rendered as `"('message',)"` instead of `"message"`. Fixed to `super().__init__(*args)`.

- **Failed external tool runs were silently ignored** (`base.py`): `subprocess.Popen(...).wait()`
  never checked the return code, so a crashed xTB, CREST, or CENSO process would be treated
  as success. Replaced with `subprocess.run(..., check=True)` at all five call sites, which
  raises `CalledProcessError` on non-zero exit.

- **Falsy energy check treated `0.0` as missing** (`parsers.py`): `if xtb_output["H"] and
  xtb_output["G"]` evaluates `0.0` as falsy, so entropy would not be computed if either
  energy happened to be exactly zero. Changed to explicit `is not None` checks in
  `parse_xtb_output` and `parse_censo_output`.

- **`sdfstream` could leak on exception** (`base.py`): The SDF file was opened with a bare
  `open()` and closed manually at the end of the conformer loop. An exception mid-loop would
  skip `sdfstream.close()`. Replaced with a `with` context manager.

### Design / Architecture

- **Removed logger-passing anti-pattern** (`parsers.py`, `chem.py`, `checkpoint.py`,
  `output.py`, `setup.py`): Every helper function previously required a `logger: Logger`
  parameter, polluting all signatures. Each module now owns a
  `logger = logging.getLogger(__name__)` at module level, consistent with standard Python
  practice. All call sites updated.

- **Eliminated double exception logging** (`base.py`, `parsers.py`, `chem.py`, `setup.py`):
  The pervasive pattern of `logger.exception(err); logger.exception("Failed to X"); raise`
  logged two entries (both with tracebacks) per exception, and caused cascading duplicates as
  exceptions bubbled up. Replaced with a single `logger.exception("Failed to X")` per catch
  site, only where additional context is added. Redundant try/except-log-raise wrappers with
  no added context were removed entirely.

- **Removed runtime config mutation** (`base.py`, `setup.py`): `_setup_pipeline()` stored
  `temp_config_dir` into `self._config["temp_config_dir"]`, turning the config dict into a
  hybrid config/state object. The temp directory path is now passed explicitly as a
  `temp_dir: Path` parameter to `create_censorc`, and derived directly from `self._args.config`
  at each call site.

- **Refactored repetitive checkpoint loaders** (`checkpoint.py`): Five near-identical
  `_load_*_results` functions (~275 lines total) reduced to a single
  `_load_stage(path, load_fn, label)` helper plus two small stage-specific parse helpers
  (~80 lines). Each stage now specifies only its file path and parse logic.

- **Removed no-op `Initiator.__init__` override** (`base.py`): The override did nothing
  except call `super().__init__()` with the same arguments. Deleted; `LarestMol.__init__`
  is now used directly.

### Minor

- **`os.makedirs` replaced with `Path.mkdir`** (`output.py`): `os.makedirs(path, exist_ok=False)`
  with a `FileExistsError` catch replaced by `path.mkdir(parents=True, exist_ok=True)`.
  Removed `import os`.

- **`config` variable shadowing in `parse_command_args`** (`parsers.py`): The loop that
  traversed the config rebound the parameter name `config` at each step. Renamed the
  traversal variable to `cfg`.

- **`summary_headings` moved outside the section loop** (`base.py`): The list was identical
  on every iteration of `for section in self._results`. Moved above the loop.
