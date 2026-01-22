# Implementation Status & Architecture: UMIErrorCorrect

This document outlines the current architecture, utilized technologies, and implementation details of the UMIErrorCorrect package. It serves as a guide for developers to understand the codebase structure and design decisions.

## Code Guidelines

- **Parallel Processing**: Uses Python's `multiprocessing.Pool` for parallel region processing.
- **BAM/SAM Handling**: All alignment operations use the `pysam` library.
- **Region Detection**: Regions can be auto-detected from BAM or defined via BED file.
- **Package Management**: Orchestration using `uv`.
- **Build System**: PEP 621 compliant `pyproject.toml` with `hatchling` build backend.
- **CLI**: Unified CLI using `typer`.
- **Logging**: Structured logging with `loguru` and `rich`.
- **Linting & Formatting**: Enforced by `ruff`.
- **Validation**: Data validation using `pydantic`.
- **Type Hints**: Extensive usage of type annotations and `mypy` for static analysis.

---

## 1. Architecture Overview

The package has been modernized to use standard Python packaging and tooling.

### Directory Structure

```text
umierrorcorrect/
├── umierrorcorrect/        # Main package source
│   ├── core/               # Core libraries and utilities
│   │   ├── consensus.py    # Consensus generation logic
│   │   ├── constants.py    # Global constants and type aliases
│   │   ├── logging_config.py # Loguru configuration
│   │   ├── umi_cluster.py  # UMI clustering algorithms
│   │   └── ...
│   ├── models/             # Pydantic data models
│   ├── cli.py              # Main CLI entry point (Typer)
│   ├── batch.py            # Batch processing logic
│   ├── glue logic files... # (align.py, preprocess.py, etc.)
│   └── version.py
├── tests/                  # Test suite
│   ├── unit/
│   ├── integration/
│   └── conftest.py
├── pyproject.toml          # Project configuration
└── uv.lock                 # Dependency lock file
```

### Build & Dependency Management

- **Build Backend**: `hatchling`
- **Dependency Manager**: `uv`
- **Configuration**: All metadata and tool configuration (Ruff, Pytest, Coverage, Mypy) is consolidated in `pyproject.toml`.

### Command Line Interface

The package provides a single entry point `umierrorcorrect` which exposes multiple subcommands via `Typer`.

**Entry Point**: `umierrorcorrect.cli:main_cli`

**Subcommands**:

- `batch`: Run the complete pipeline (preprocessing -> mapping -> consensus -> stats -> variants).
- `preprocess`: UMI extraction and adapter trimming.
- `mapping`: BWA alignment.
- `consensus`: Generate consensus sequences (error correction).
- `stats`: Generate consensus statistics.
- `variants`: Call variants from consensus reads.
- `filter-bam`: Filter BAM files by consensus depth.
- `downsampling`: Analysis of downsampling.

---

## 2. Pipeline Implementation

The pipeline processes raw sequencing data into error-corrected consensus reads and variants.

### Pipeline Flow Diagrams

The preprocessing pipeline handles adapter trimming and UMI extraction.

#### Overview: Complete Pipeline

```text
┌─────────────────────────────────────────────────────────────────────────────┐
│                         UMI Error Correct Pipeline                          │
└─────────────────────────────────────────────────────────────────────────────┘

Input: FASTQ R1 (+ R2 for paired-end)
                    │
                    ▼
        ┌───────────────────────┐
        │   PREPROCESSING       │ ◄── Handles adapter trimming & UMI extraction
        │   (fastp or cutadapt) │
        └───────────────────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │   BWA MAPPING         │
        │   (align_bwa)         │
        └───────────────────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │   UMI ERROR CORRECT   │
        │   - UMI clustering    │
        │   - Consensus gen     │
        └───────────────────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │   CONSENSUS STATS     │
        └───────────────────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │   VARIANT CALLING     │
        └───────────────────────┘
                    │
                    ▼
Output: VCF, consensus BAM, statistics
```

#### Preprocessing Logic

Configuration flags determine whether `fastp` (default, faster) or `cutadapt` is used.

##### **Path A: With fastp (Default)**

- Quality filtering (Q20 default)
- Adapter trimming
- Read merging (optional)
- UMI extraction

##### **Path B: Without fastp (`--no-fastp`)**

- Adapter trimming via `cutadapt` (if enabled)
- No quality filtering
- UMI extraction

---

## 3. Core Components

### 3.1 UMI Clustering & Error Correction

**Location**: `umierrorcorrect/core/umi_cluster.py` & `umierrorcorrect/umi_error_correct.py`

- **Algorithm**: Adjacency-based clustering using Hamming or Levenshtein distance.
- **Optimization**: Accesses `pysam` objects directly and uses `multiprocessing` to handle genomic regions in parallel.
- **Security**: Subprocess calls (e.g., for sorting/uniq) use `shutil.which` and avoid `shell=True` to prevent command injection. Secure temporary files are used.

### 3.2 Consensus Generation

**Location**: `umierrorcorrect/core/consensus.py`

- Generates consensus sequences from clustered reads.
- Supports position-based grouping.
- Filters consensus bases based on frequency thresholds (default 60%).

### 3.3 Logging System

**Location**: `umierrorcorrect/core/logging_config.py`

- Implemented using `loguru`.
- **Console**: Uses `rich` for pretty formatting.
- **File**: Rotated, compressed log files saved to the output directory.
- **Context**: Captures function names and line numbers for debugging.

---

## 4. ctDNA Sample Analysis Data Models

**Location**: `umierrorcorrect/models/models.py`

The package supports personalized ctDNA (circulating tumor DNA) analysis through a hierarchical data model that organizes patients, their tracked mutations, and longitudinal samples.

### 4.1 Data Model Hierarchy

```text
SampleSheet (CSV parser)
└── patients: list[Patient]

Patient
├── patient_id: str
├── mutation_bed: Path
├── mutations: list[Mutation]      ← Auto-loaded from BED file
└── samples: list[Sample]
    └── results: list[MutationResult]  ← Populated after analysis
```

### 4.2 Model Classes

#### `Mutation`

Represents a patient-specific mutation to track, parsed from a mutation BED file.

| Field | Type | Description |
|-------|------|-------------|
| `name` | `str` | Mutation identifier (e.g., `KRAS_G12D`) |
| `chromosome` | `str` | Chromosome (e.g., `chr3`) |
| `position` | `int` | 0-based genomic position |
| `reference` | `str` | Reference allele |
| `alternate` | `str` | Alternate allele |

**Mutation BED format** (tab-separated):

```text
# chrom  start  end  name  ref  alt
chr14   22982082   22982083   AJUBA_c.184G>T   G   T
chr11   108250700  108250701  ATM_c.1236G>T    G   T
```

#### `MutationResult`

Stores analysis results for a tracked mutation in a specific sample.

| Field | Type | Description |
|-------|------|-------------|
| `mutation_name` | `str` | Links to `Mutation.name` |
| `vaf` | `float \| None` | Variant allele frequency (0.0–1.0) |
| `mm_count` | `int \| None` | Mutant molecule count |
| `wt_count` | `int \| None` | Wild-type molecule count |
| `total_reads` | `int \| None` | Total reads at position |
| `mm_per_ml` | `float \| None` | Mutant molecules per mL plasma |

#### `Sample`

Represents a single sequencing sample with FASTQ files and metadata.

| Field | Type | Description |
|-------|------|-------------|
| `name` | `str` | Sample identifier |
| `read1` | `Path` | Path to R1 FASTQ file |
| `read2` | `Path \| None` | Path to R2 FASTQ file (paired-end) |
| `sample_type` | `"validation" \| "ctdna" \| None` | Sample category |
| `replicate` | `int \| None` | Replicate number |
| `ng_input` | `float \| None` | DNA input in nanograms |
| `ml_plasma` | `float \| None` | Plasma volume in milliliters |
| `region_bed` | `Path \| None` | Target regions BED file |
| `collection_date` | `str \| None` | Collection date (YYYY-MM-DD) |
| `results` | `list[MutationResult]` | Analysis results (populated post-processing) |

#### `Patient`

Groups mutations and samples for a single patient.

| Field | Type | Description |
|-------|------|-------------|
| `patient_id` | `str` | Patient identifier |
| `mutation_bed` | `Path` | Path to mutation BED file |
| `mutations` | `list[Mutation]` | Parsed mutations (auto-loaded) |
| `samples` | `list[Sample]` | Associated samples |

**Helper methods**:

- `get_sample(name)` → Get sample by name
- `get_samples_by_type(type)` → Filter by sample type
- `get_samples_sorted_by_date()` → Chronological order for longitudinal analysis

#### `SampleSheet`

Parses and validates a CSV sample sheet, grouping samples by patient.

| Field | Type | Description |
|-------|------|-------------|
| `csv_path` | `Path` | Path to CSV file |
| `patients` | `list[Patient]` | Parsed patients with samples |
| `base_path` | `Path \| None` | Base path for resolving relative paths |

### 4.3 Sample Sheet CSV Format

**Required columns**: `patient_id`, `mutation_bed`, `sample_name`, `read1`

**Optional columns**: `read2`, `sample_type`, `replicate`, `ng_input`, `ml_plasma`, `region_bed`, `collection_date`

```csv
patient_id,mutation_bed,sample_name,read1,read2,sample_type,replicate,ng_input,ml_plasma,region_bed,collection_date
P001,/data/P001_mutations.bed,P001_T0,/data/P001_T0_R1.fq.gz,/data/P001_T0_R2.fq.gz,ctdna,1,17.5,3.0,/data/regions.bed,2024-01-15
P001,/data/P001_mutations.bed,P001_T1,/data/P001_T1_R1.fq.gz,/data/P001_T1_R2.fq.gz,ctdna,1,49.1,3.0,/data/regions.bed,2024-02-10
P002,/data/P002_mutations.bed,P002_T0,/data/P002_T0_R1.fq.gz,/data/P002_T0_R2.fq.gz,ctdna,1,18.9,2.5,/data/regions.bed,2024-01-20
```

**Validation rules**:

- All required columns must be present
- `patient_id` cannot be empty
- `mutation_bed` must be consistent for the same patient across rows
- File paths are validated on `Sample` instantiation
- `collection_date` must be in `YYYY-MM-DD` format

### 4.4 Usage Example

```python
from pathlib import Path
from umierrorcorrect.models.models import SampleSheet, MutationResult

# Load sample sheet (mutations auto-loaded from BED files)
sheet = SampleSheet(csv_path=Path("samples.csv"), base_path=Path("/data"))

print(f"Loaded {sheet.patient_count} patients, {sheet.sample_count} samples")

for patient in sheet.patients:
    print(f"\nPatient {patient.patient_id}:")
    print(f"  Tracking {len(patient.mutations)} mutations")

    for sample in patient.get_samples_sorted_by_date():
        print(f"  - {sample.name} ({sample.collection_date})")

        # After analysis, results are populated:
        for mutation in patient.mutations:
            result = sample.get_result(mutation.name)
            if result and result.vaf is not None:
                print(f"      {mutation.name}: VAF={result.vaf:.4f}")
```

---

## 5. Development & Testing Infrastructure

### 5.1 Development Setup

Using `uv` for fast, reproducible environments.

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Sync dependencies
uv sync --dev

# Run tests
uv run pytest
```

### 5.2 Test Suite

**Location**: `tests/`

- **Framework**: `pytest`
- **Unit Tests**: `tests/unit/` - Test individual functions and classes (formatting, clustering logic, args).
- **Integration Tests**: `tests/integration/` - Test full pipeline steps (preprocessing, alignment).
- **Fixtures**: Shared test data and artifacts defined in `tests/conftest.py`.

### 5.3 Code Quality Tools

All configured in `pyproject.toml`.

- **Linting**: `ruff check .`
- **Formatting**: `ruff format .`
- **Type Checking**: `mypy .`
- **Coverage**: `pytest --cov=umierrorcorrect`

---

## 6. Completed Modernization Tasks

The following improvements have been fully implemented:

1. **Refactored to `src` layout (implicit)**: Code moved to `umierrorcorrect/` package structure.
2. **Modern Packaging**: Replaced `setup.py` with `pyproject.toml`.
3. **Typer CLI**: Replaced `argparse` with `Typer` for a robust, self-documenting CLI.
4. **Security Fixes**: Removed vulnerable `shell=True` usage in subprocess calls.
5. **Resource Management**: Improved temporary file handling with `tempfile` and `contextlib`.
6. **Type Safety**: Added type hints to core modules (`constants.py`, `umi_cluster.py`, etc.).
7. **Constants**: Externalized magic numbers to `core/constants.py`.
