# Implementation Plan: UMIErrorCorrect Modernization

This document outlines a comprehensive plan to modernize the UMIErrorCorrect codebase, improve code quality, and establish proper testing infrastructure.

## Executive Summary

The codebase is functional and well-structured for its purpose (bioinformatics UMI error correction), but shows signs of organic growth without modern Python tooling. Key issues include:

- **No type hints** throughout the codebase
- **No test suite** (only manual integration testing)
- **Legacy build system** (setup.py without pyproject.toml)
- **Security vulnerability** (shell=True in subprocess)
- **Code quality issues** (long functions, hardcoded values, inconsistent naming)

**Estimated effort:** The full modernization is a significant undertaking. Prioritize by phase.

---

## Phase 1: Modern Build System and Tooling (Foundation)

### 1.1 Create pyproject.toml with Hatchling

Replace `setup.py` and `setup.cfg` with a modern `pyproject.toml`:

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "umierrorcorrect"
dynamic = ["version"]
description = "Pipeline for analyzing barcoded amplicon sequencing data with Unique Molecular Identifiers (UMI)"
readme = "README.md"
license = "MIT"
requires-python = ">=3.9"
authors = [
    { name = "Tobias Osterlund", email = "tobias.osterlund@gu.se" }
]
classifiers = [
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
]
dependencies = [
    "pysam>=0.8.4",
    "scipy",
    "matplotlib",
]

[project.optional-dependencies]
dev = [
    "ruff",
    "pytest",
    "pytest-cov",
    "mypy",
]

[project.scripts]
run_umierrorcorrect = "umierrorcorrect.run_umierrorcorrect:main_cli"
preprocess = "umierrorcorrect.preprocess:main_cli"
run_mapping = "umierrorcorrect.run_mapping:main_cli"
umi_error_correct = "umierrorcorrect.umi_error_correct:main_cli"
get_consensus_statistics = "umierrorcorrect.get_consensus_statistics:main_cli"
call_variants = "umierrorcorrect.call_variants:main_cli"
filter_bam = "umierrorcorrect.filter_bam:main_cli"
filter_cons = "umierrorcorrect.filter_cons:main_cli"
downsampling_plots = "umierrorcorrect.downsampling_plots:main_cli"
fit_background_model = "umierrorcorrect.fit_background_model:main_cli"

[project.urls]
Homepage = "https://github.com/stahlberggroup/umierrorcorrect"
Documentation = "https://github.com/stahlberggroup/umierrorcorrect/wiki"
Repository = "https://github.com/stahlberggroup/umierrorcorrect"

[tool.hatch.version]
path = "umierrorcorrect/version.py"

[tool.hatch.build.targets.sdist]
include = [
    "/umierrorcorrect",
    "/test_data",
    "/doc",
]

[tool.hatch.build.targets.wheel]
packages = ["umierrorcorrect"]
```

### 1.2 Configure Ruff for Linting and Formatting

Add to `pyproject.toml`:

```toml
[tool.ruff]
target-version = "py39"
line-length = 120
src = ["umierrorcorrect"]

[tool.ruff.lint]
select = [
    "E",      # pycodestyle errors
    "W",      # pycodestyle warnings
    "F",      # Pyflakes
    "I",      # isort
    "B",      # flake8-bugbear
    "C4",     # flake8-comprehensions
    "UP",     # pyupgrade
    "ARG",    # flake8-unused-arguments
    "SIM",    # flake8-simplify
    "S",      # flake8-bandit (security)
    "PTH",    # flake8-use-pathlib
]
ignore = [
    "E501",   # line too long (handled by formatter)
    "S101",   # assert usage (ok in tests)
]
# Always auto-fix when possible
fixable = ["ALL"]

[tool.ruff.lint.per-file-ignores]
"tests/*" = ["S101", "ARG"]

[tool.ruff.format]
quote-style = "double"
indent-style = "space"
```

### 1.3 Set Up uv for Package Management

Create a `uv.lock` file and development workflow:

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Initialize project with uv
uv sync

# Development installation
uv sync --dev

# Run tools via uv
uv run ruff check --fix .
uv run ruff format .
uv run pytest
```

Add to documentation:

```markdown
## Development Setup

1. Install [uv](https://docs.astral.sh/uv/):
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

2. Clone and install:

   ```bash
   git clone https://github.com/stahlberggroup/umierrorcorrect.git
   cd umierrorcorrect
   uv sync --dev
   ```

3. Run tests:

   ```bash
   uv run pytest
   ```
```

### 1.4 Migration Tasks

| Task | Description | Breaking Change |
|------|-------------|-----------------|
| Refactor entry points | Change from `scripts` to `project.scripts` with proper `main_cli()` functions | No (same CLI names) |
| Remove setup.py | Delete after pyproject.toml is validated | No |
| Remove setup.cfg | Delete after pyproject.toml is validated | No |
| Update version.py | Ensure compatible with hatch version reading | No |

**Entry point refactoring example** (for each script):

```python
# Current (run_umierrorcorrect.py):
if __name__ == '__main__':
    args = parseArgs()
    main(args)

# New:
def main_cli():
    """CLI entry point."""
    args = parseArgs()
    main(args)

if __name__ == '__main__':
    main_cli()
```

---

## Phase 2: Test Suite Infrastructure

### 2.1 Test Directory Structure

```text
tests/
├── __init__.py
├── conftest.py              # Shared fixtures
├── test_data/               # Small test fixtures (symlink or copy)
│   ├── small_R1.fastq.gz    # Minimal FASTQ for unit tests
│   ├── small.bam            # Minimal BAM for unit tests
│   └── regions.bed
├── unit/
│   ├── __init__.py
│   ├── test_umi_cluster.py
│   ├── test_get_consensus.py
│   ├── test_get_cons_info.py
│   ├── test_group.py
│   ├── test_get_regions_from_bed.py
│   ├── test_check_args.py
│   └── test_handle_sequences.py
└── integration/
    ├── __init__.py
    ├── test_preprocess.py
    ├── test_run_mapping.py
    ├── test_umi_error_correct.py
    ├── test_pipeline.py      # End-to-end pipeline test
    └── test_call_variants.py
```

### 2.2 Pytest Configuration

Add to `pyproject.toml`:

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
python_functions = ["test_*"]
addopts = [
    "-v",
    "--tb=short",
    "--strict-markers",
]
markers = [
    "slow: marks tests as slow (deselect with '-m \"not slow\"')",
    "integration: marks tests as integration tests",
    "requires_bwa: marks tests that require bwa installed",
]
filterwarnings = [
    "ignore::DeprecationWarning:pysam.*",
]
```

### 2.3 Test Fixtures (conftest.py)

```python
import pytest
from pathlib import Path
import tempfile
import shutil

@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent / "test_data"

@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for test outputs."""
    tmpdir = tempfile.mkdtemp()
    yield Path(tmpdir)
    shutil.rmtree(tmpdir)

@pytest.fixture
def sample_umi_groups():
    """Sample UMI barcode groups for clustering tests."""
    return {
        "ACGTACGTACGT": ["read1", "read2", "read3"],
        "ACGTACGTACGA": ["read4"],  # 1 edit distance from first
        "GGGGGGGGGGGG": ["read5", "read6"],
    }

@pytest.fixture
def sample_reads():
    """Sample aligned reads for consensus testing."""
    # Return mock pysam AlignedSegment objects or dict representations
    pass
```

### 2.4 Priority Test Cases

**Unit tests for core algorithms (highest priority):**

| Module | Function | Test Cases |
|--------|----------|------------|
| `umi_cluster.py` | `hamming_distance()` | Identical strings, single diff, all diff, different lengths |
| `umi_cluster.py` | `edit_distance()` | Same as above plus insertions/deletions |
| `umi_cluster.py` | `cluster_barcodes()` | Various distance thresholds, edge cases |
| `umi_cluster.py` | `get_connected_components()` | Graph connectivity scenarios |
| `get_consensus3.py` | `getConsensus3()` | Uniform reads, mixed reads, low quality, indels |
| `get_cons_info.py` | `get_cons_info()` | Various coverage depths, variant positions |
| `get_regions_from_bed.py` | `read_bed()` | Valid BED, empty, malformed |
| `get_regions_from_bed.py` | `merge_regions()` | Overlapping, adjacent, disjoint regions |

**Integration tests (medium priority):**

| Test | Description |
|------|-------------|
| `test_preprocess.py` | FASTQ preprocessing with various UMI configurations |
| `test_pipeline.py` | End-to-end with test_data files |

### 2.5 Coverage Configuration

Add to `pyproject.toml`:

```toml
[tool.coverage.run]
source = ["umierrorcorrect"]
branch = true
omit = [
    "*/test_*.py",
    "*/__init__.py",
]

[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "if __name__ == .__main__.:",
    "raise NotImplementedError",
]
show_missing = true
fail_under = 60  # Start low, increase over time
```

---

## Phase 3: Critical Bug Fixes and Security

### 3.1 CRITICAL: Fix Command Injection Vulnerability

**Location:** `umierrorcorrect/umi_error_correct.py` lines 227-229

**Current (vulnerable):**
```python
command1 = ['sort tmp.txt | uniq -d']
p1 = subprocess.Popen(command1, shell=True, stdout=g)
```

**Fixed:**
```python
import shutil

# Use proper subprocess chaining without shell=True
sort_cmd = [shutil.which("sort") or "sort", tmp_file.name]
uniq_cmd = [shutil.which("uniq") or "uniq", "-d"]

p1 = subprocess.Popen(sort_cmd, stdout=subprocess.PIPE)
p2 = subprocess.Popen(uniq_cmd, stdin=p1.stdout, stdout=g)
p1.stdout.close()
p2.communicate()
```

### 3.2 Fix Resource Leak

**Location:** `umierrorcorrect/src/check_args.py` line 37-38

**Current:**
```python
devnull = open(os.devnull)
subprocess.Popen([name,'--version'], stdout=devnull, stderr=devnull)
```

**Fixed:**
```python
with open(os.devnull, 'w') as devnull:
    subprocess.run([name, '--version'], stdout=devnull, stderr=devnull, check=False)
```

### 3.3 Fix Unreachable Code

**Location:** `umierrorcorrect/src/check_args.py` lines 77-83

**Current:**
```python
except ValueError as e:
    raise(e + " Barcode length needs to be an integer")
    sys.exit(1)  # Unreachable
```

**Fixed:**
```python
except ValueError as e:
    raise ValueError(f"Barcode length must be an integer: {e}") from e
```

### 3.4 Use Temporary Files Properly

**Location:** `umierrorcorrect/umi_error_correct.py` lines 219, 228

**Current:**
```python
f = open('tmp.txt', 'w')
```

**Fixed:**
```python
import tempfile

with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
    tmp_file = f.name
    # ... use f ...
# Clean up later
os.unlink(tmp_file)
```

---

## Phase 4: Code Quality Improvements

### 4.1 Add Type Hints

**Priority files** (core algorithms):

```python
# umierrorcorrect/src/umi_cluster.py

from typing import Dict, List, Set, Tuple

def hamming_distance(s1: str, s2: str) -> int:
    """Calculate Hamming distance between two strings of equal length."""
    ...

def edit_distance(s1: str, s2: str) -> int:
    """Calculate Levenshtein edit distance between two strings."""
    ...

def cluster_barcodes(
    barcode_counts: Dict[str, int],
    distance_threshold: int = 1
) -> Dict[str, List[str]]:
    """Cluster barcodes by edit distance."""
    ...
```

Add mypy configuration to `pyproject.toml`:

```toml
[tool.mypy]
python_version = "3.9"
warn_return_any = true
warn_unused_configs = true
ignore_missing_imports = true
exclude = ["tests/", "build/"]

[[tool.mypy.overrides]]
module = "pysam.*"
ignore_missing_imports = true
```

### 4.2 Extract Constants

Create `umierrorcorrect/constants.py`:

```python
"""Constants used throughout the UMIErrorCorrect pipeline."""

# Default family sizes for consensus statistics
DEFAULT_FAMILY_SIZES: list[int] = [0, 1, 2, 3, 4, 5, 7, 10, 20, 30]

# Maximum reads per chunk for parallel processing
MAX_READS_PER_CHUNK: int = 100_000

# Default edit distance threshold for UMI clustering
DEFAULT_EDIT_DISTANCE_THRESHOLD: int = 1

# Default position threshold for grouping
DEFAULT_POSITION_THRESHOLD: int = 10

# Default consensus frequency threshold (percent)
DEFAULT_CONSENSUS_FREQUENCY: float = 60.0

# Default indel frequency threshold (percent)
DEFAULT_INDEL_FREQUENCY: float = 60.0
```

### 4.3 Refactor Long Functions

**Target: `get_consensus3.py:getConsensus3()` (170 lines)**

Split into:
- `_initialize_consensus_read()` - Set up initial state
- `_process_alignment_column()` - Handle single column
- `_apply_quality_thresholds()` - Filter by quality
- `_build_cigar_string()` - Construct CIGAR
- `getConsensus3()` - Orchestrate the above

**Target: `umi_error_correct.py:cluster_consensus_worker()` (70 lines)**

Split into:
- `_cluster_umis()` - Just the clustering logic
- `_generate_consensus_for_cluster()` - Consensus generation
- `_write_cluster_output()` - File I/O
- `cluster_consensus_worker()` - Orchestrate

### 4.4 Modernize Path Handling

Replace string concatenation with pathlib:

```python
# Before
output_file = output_path + '/' + sample_name + '_consensus_reads.bam'
basename = filename.split('/')[-1]

# After
from pathlib import Path

output_file = Path(output_path) / f"{sample_name}_consensus_reads.bam"
basename = Path(filename).name
```

### 4.5 Modernize String Formatting

Replace `.format()` with f-strings:

```python
# Before
parser = argparse.ArgumentParser(description="UmiErrorCorrect v. {}. Pipeline...".format(__version__))

# After
parser = argparse.ArgumentParser(description=f"UmiErrorCorrect v. {__version__}. Pipeline...")
```

### 4.6 Remove Python 2 Compatibility Code

Remove from all files:
```python
from __future__ import division  # Not needed in Python 3
```

### 4.7 Fix Typos in User-Facing Text

| File | Line | Current | Fixed |
|------|------|---------|-------|
| `run_umierrorcorrect.py` | 44 | "emove the original" | "Remove the original" |
| `umi_error_correct.py` | 54 | "emove the original" | "Remove the original" |
| `run_mapping.py` | 26 | "emove the original" | "Remove the original" |
| `umi_error_correct.py` | 491 | "0cluster umis" | "cluster UMIs" |

---

## Phase 5: Pre-commit Hooks

Create `.pre-commit-config.yaml`:

```yaml
repos:
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.8.0
    hooks:
      - id: ruff
        args: [--fix]
      - id: ruff-format

  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v5.0.0
    hooks:
      - id: trailing-whitespace
      - id: end-of-file-fixer
      - id: check-yaml
      - id: check-added-large-files
```

---

## Phase 6: Documentation Updates

### 6.1 Add Docstrings

Use Google-style docstrings:

```python
def cluster_barcodes(
    barcode_counts: Dict[str, int],
    distance_threshold: int = 1
) -> Dict[str, List[str]]:
    """Cluster UMI barcodes by edit distance.

    Groups barcodes that are within the specified edit distance
    threshold, with the most abundant barcode as the cluster center.

    Args:
        barcode_counts: Dictionary mapping barcode sequences to their counts.
        distance_threshold: Maximum edit distance for clustering (default: 1).

    Returns:
        Dictionary mapping cluster center barcodes to lists of member barcodes.

    Example:
        >>> counts = {"ACGT": 10, "ACGA": 2, "GGGG": 5}
        >>> clusters = cluster_barcodes(counts, distance_threshold=1)
        >>> clusters
        {"ACGT": ["ACGT", "ACGA"], "GGGG": ["GGGG"]}
    """
```

### 6.2 Update README.md

Add sections for:
- Development setup with uv
- Running tests
- Contributing guidelines
- Link to this implementation plan

---

## Implementation Priority Matrix

| Phase | Priority | Effort | Impact | Dependencies |
|-------|----------|--------|--------|--------------|
| Phase 1 (Build System) | HIGH | Medium | High | None |
| Phase 2 (Test Suite) | HIGH | High | High | Phase 1 |
| Phase 3 (Security Fixes) | CRITICAL | Low | Critical | None |
| Phase 4 (Code Quality) | MEDIUM | High | Medium | Phase 1, 2 |
| Phase 5 (Pre-commit Hooks) | MEDIUM | Low | Medium | Phase 1 |
| Phase 6 (Documentation) | LOW | Medium | Medium | Phase 4 |

## Recommended Implementation Order

1. **Immediate (Week 1):**
   - Phase 3.1: Fix command injection vulnerability
   - Phase 3.2-3.4: Fix other critical bugs
   - Phase 1.1-1.2: Set up pyproject.toml with hatchling and ruff

2. **Short-term (Weeks 2-3):**
   - Phase 1.3-1.4: Complete build system migration
   - Phase 2.1-2.3: Set up test infrastructure
   - Phase 5: Set up pre-commit hooks

3. **Medium-term (Weeks 4-6):**
   - Phase 2.4-2.5: Write priority unit tests
   - Phase 4.1: Add type hints to core modules
   - Phase 4.2: Extract constants

4. **Long-term (Ongoing):**
   - Phase 4.3-4.6: Refactoring and modernization
   - Phase 6: Documentation improvements
   - Increase test coverage to 80%+

---

## Files to Create

| File | Purpose |
|------|---------|
| `pyproject.toml` | Modern build configuration |
| `.pre-commit-config.yaml` | Pre-commit hooks |
| `tests/conftest.py` | Test fixtures |
| `tests/unit/test_umi_cluster.py` | Unit tests |
| `umierrorcorrect/constants.py` | Extracted constants |

## Files to Delete (after migration)

| File | Reason |
|------|--------|
| `setup.py` | Replaced by pyproject.toml |
| `setup.cfg` | Replaced by pyproject.toml |

## Files to Modify

| File | Changes |
|------|---------|
| All scripts in `umierrorcorrect/` | Add `main_cli()` entry points |
| `umierrorcorrect/umi_error_correct.py` | Security fixes, refactoring |
| `umierrorcorrect/src/check_args.py` | Bug fixes, resource leak |
| `README.md` | Development instructions |
| `CLAUDE.md` | Update build/test commands |

---

## Appendix: Current State Assessment

### Strengths
- Well-structured pipeline with clear separation of concerns
- Modular design allows running steps independently
- Published and validated algorithm (Clinical Chemistry 2022)
- Docker support for reproducibility
- Reasonable documentation for end users

### Weaknesses
- No automated testing (only manual integration test)
- No type hints (harder to maintain, no IDE support)
- Legacy build system (setup.py)
- Security vulnerability (shell=True)
- Code debt (long functions, hardcoded values, old patterns)
- No CI/CD pipeline
- Inconsistent code style

### Risk Assessment
- **High Risk:** Command injection vulnerability could allow arbitrary code execution
- **Medium Risk:** Lack of tests means regressions go undetected
- **Low Risk:** Code style issues affect maintainability but not functionality
