# Implementation Plan: ctDNA Analysis Module

## Overview

Add ctDNA-specific analysis capabilities to umierrorcorrect while keeping the generic UMI pipeline intact.

### Goals

1. Keep current generic UMI pipeline unchanged
2. Add `umierrorcorrect analysis --ctdna` command for post-hoc analysis
3. Add `--analysis-mode ctdna` flag to batch for integrated processing
4. Export results as TSV (reporting module deferred)
5. Family size configurable, default 3

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                   Generic UMI Pipeline (unchanged)                       │
│                                                                          │
│   umierrorcorrect batch → [preprocess → align → consensus → variants]   │
│                                ↓                                         │
│                      ProcessingResult objects                            │
│                      (paths to .cons, .vcf, .bam)                        │
└─────────────────────────────────────────────────────────────────────────┘
                                 ↓
                    (if --analysis-mode ctdna OR run analysis command)
                                 ↓
┌─────────────────────────────────────────────────────────────────────────┐
│                      ctDNA Analysis Module (NEW)                         │
│                                                                          │
│   Inputs:                                                                │
│     - ProcessingResult OR results directory from previous run            │
│     - Extended sample sheet (with patient_id, mutation_bed, etc.)        │
│                                                                          │
│   Processing:                                                            │
│     1. Load mutations from mutation BED file(s)                          │
│     2. Parse .cons files at mutation positions                           │
│     3. Extract counts at specified family size threshold                 │
│     4. Calculate VAF, mm_count, mm_per_ml                               │
│     5. Populate Sample.results with MutationResult objects               │
│                                                                          │
│   Outputs:                                                               │
│     - TSV file with mutation results per sample                          │
│     - TSV file with patient-level summary (if multiple samples)          │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## File Structure

```text
umierrorcorrect/
├── cli.py                      # Add 'analysis' command
├── batch.py                    # Add --analysis-mode option, hook ctdna after processing
├── models/
│   └── models.py               # Minor updates to SampleSheet parsing
│
└── ctdna/                      # NEW PACKAGE
    ├── __init__.py
    ├── analyzer.py             # CtdnaAnalyzer: main orchestrator
    ├── mutation_tracker.py     # Parse .cons files for specific mutations
    ├── sample_sheet.py         # Extended sample sheet parser (ctDNA-specific)
    └── export.py               # TSV export functions
```

---

## Implementation Steps

### Phase 1: Core ctDNA Module

#### 1.1 Create `umierrorcorrect/ctdna/__init__.py`

```python
"""ctDNA analysis module for UMI Error Correct."""
from umierrorcorrect.ctdna.analyzer import CtdnaAnalyzer
from umierrorcorrect.ctdna.mutation_tracker import MutationTracker

__all__ = ["CtdnaAnalyzer", "MutationTracker"]
```

#### 1.2 Create `umierrorcorrect/ctdna/mutation_tracker.py`

Responsible for extracting mutation data from .cons files. The mutations
were orignallay provided by the mutation bed file.

```python
"""Track mutations in consensus data."""
from pathlib import Path
from umierrorcorrect.models.models import Mutation, MutationResult

class MutationTracker:
    """Extract mutation information from .cons files."""

    def __init__(self, family_size: int = 3):
        self.family_size = family_size

    def track_mutation(
        self,
        cons_file: Path,
        mutation: Mutation,
        ml_plasma: float | None = None
    ) -> MutationResult:
        """
        Find mutation position in .cons file and calculate metrics.

        The .cons file has columns:
        - sample, contig, pos, annotation, ref, A, C, G, T, I, D, N,
          coverage, fsize, max_nonref_count, max_nonref_freq, max_nonref_allele

        Multiple rows exist per position (different family sizes).
        We select the row matching our family_size threshold.
        """
        # Parse cons file, find position, extract metrics
        # Return MutationResult with vaf, mm_count, total_reads, mm_per_ml
        ...

    def track_mutations(
        self,
        cons_file: Path,
        mutations: list[Mutation],
        ml_plasma: float | None = None
    ) -> list[MutationResult]:
        """Track multiple mutations in a single cons file."""
        return [self.track_mutation(cons_file, m, ml_plasma) for m in mutations]
```

#### 1.3 Create `umierrorcorrect/ctdna/sample_sheet.py`

Extended sample sheet parser for ctDNA workflows.

```python
"""Extended sample sheet parsing for ctDNA analysis."""
from pathlib import Path
from umierrorcorrect.models.models import Patient, Sample, Mutation, SampleSheet

def parse_ctdna_sample_sheet(csv_path: Path, base_path: Path | None = None) -> list[Patient]:
    """
    Parse extended sample sheet for ctDNA analysis.

    Required columns: sample_name, read1, patient_id, mutation_bed
    Optional columns: read2, sample_type, replicate, ng_input, ml_plasma,
                      region_bed, collection_date

    Returns list of Patient objects with their samples and mutations loaded.
    """
    # Use existing SampleSheet class from models.py but handle missing patient_id gracefully
    ...

def parse_ctdna_sample_sheet_for_existing_results(
    csv_path: Path,
    results_dir: Path
) -> list[Patient]:
    """
    Parse sample sheet and associate with existing results.

    For each sample, locate the corresponding .cons file in results_dir.
    Validates that expected output files exist.
    """
    ...
```

#### 1.4 Create `umierrorcorrect/ctdna/analyzer.py`

Main orchestrator for ctDNA analysis.

```python
"""ctDNA analysis orchestrator."""
from dataclasses import dataclass
from pathlib import Path
from umierrorcorrect.models.models import Patient, Sample, MutationResult
from umierrorcorrect.ctdna.mutation_tracker import MutationTracker
from umierrorcorrect.ctdna.export import export_results_tsv

@dataclass
class CtdnaAnalysisResult:
    """Result of ctDNA analysis for a sample."""
    sample_name: str
    patient_id: str
    cons_file: Path
    results: list[MutationResult]
    success: bool
    error_message: str | None = None

class CtdnaAnalyzer:
    """Orchestrate ctDNA analysis across patients and samples."""

    def __init__(self, family_size: int = 3):
        self.family_size = family_size
        self.tracker = MutationTracker(family_size=family_size)

    def analyze_sample(
        self,
        sample: Sample,
        mutations: list[Mutation],
        cons_file: Path
    ) -> CtdnaAnalysisResult:
        """Analyze a single sample for tracked mutations."""
        results = self.tracker.track_mutations(
            cons_file=cons_file,
            mutations=mutations,
            ml_plasma=sample.ml_plasma
        )
        # Populate sample.results
        sample.results = results
        return CtdnaAnalysisResult(...)

    def analyze_patient(self, patient: Patient, results_dir: Path) -> list[CtdnaAnalysisResult]:
        """Analyze all samples for a patient."""
        results = []
        for sample in patient.samples:
            cons_file = self._find_cons_file(sample.name, results_dir)
            result = self.analyze_sample(sample, patient.mutations, cons_file)
            results.append(result)
        return results

    def analyze_all(
        self,
        patients: list[Patient],
        results_dir: Path,
        output_dir: Path
    ) -> Path:
        """
        Analyze all patients and export results.
        Returns path to output TSV.
        """
        all_results = []
        for patient in patients:
            all_results.extend(self.analyze_patient(patient, results_dir))

        output_file = output_dir / "ctdna_analysis_results.tsv"
        export_results_tsv(all_results, patients, output_file)
        return output_file

    def _find_cons_file(self, sample_name: str, results_dir: Path) -> Path:
        """Locate .cons file for a sample in results directory."""
        # Search in results_dir/samples/{sample_name}/ or results_dir/
        ...
```

#### 1.5 Create `umierrorcorrect/ctdna/export.py`

TSV export functions.

```python
"""Export ctDNA analysis results."""
from pathlib import Path
from umierrorcorrect.models.models import Patient
from umierrorcorrect.ctdna.analyzer import CtdnaAnalysisResult

def export_results_tsv(
    results: list[CtdnaAnalysisResult],
    patients: list[Patient],
    output_file: Path
) -> None:
    """
    Export analysis results to TSV.

    Columns:
    - patient_id
    - sample_name
    - sample_type
    - collection_date
    - mutation_name
    - chromosome
    - position
    - ref
    - alt
    - vaf
    - mm_count (mutant molecule count)
    - total_reads
    - mm_per_ml (if ml_plasma available)
    - ml_plasma
    - ng_input
    - family_size_threshold
    """
    ...

def export_patient_summary_tsv(
    patients: list[Patient],
    output_file: Path
) -> None:
    """
    Export patient-level summary.

    One row per patient with aggregated metrics.
    """
    ...
```

---

### Phase 2: CLI Integration

#### 2.1 Add `analysis` command to `cli.py`

```python
@app.command()
def analysis(
    sample_sheet: Annotated[
        Path,
        typer.Option("--sample-sheet", "-s", help="Extended sample sheet CSV with patient/mutation info.")
    ],
    results_dir: Annotated[
        Path,
        typer.Option("--results-dir", "-r", help="Directory containing umierrorcorrect output.")
    ],
    output_dir: Annotated[
        Optional[Path],
        typer.Option("-o", "--output-dir", help="Output directory for analysis results.")
    ] = None,
    ctdna: Annotated[
        bool,
        typer.Option("--ctdna", help="Run ctDNA mutation tracking analysis.")
    ] = False,
    family_size: Annotated[
        int,
        typer.Option("-f", "--family-size", help="Family size threshold for analysis.")
    ] = 3,
) -> None:
    """Run secondary analysis on umierrorcorrect results.

    This command performs analysis on existing umierrorcorrect output.
    Currently supports ctDNA mutation tracking analysis.

    Examples:
        # Run ctDNA analysis on existing results
        umierrorcorrect analysis --ctdna -s samples.csv -r results/ -o analysis/

        # With custom family size threshold
        umierrorcorrect analysis --ctdna -s samples.csv -r results/ -f 5
    """
    if not ctdna:
        console.print("[red]Error:[/red] Please specify an analysis mode (--ctdna)")
        raise typer.Exit(1)

    if ctdna:
        from umierrorcorrect.ctdna import CtdnaAnalyzer
        from umierrorcorrect.ctdna.sample_sheet import parse_ctdna_sample_sheet_for_existing_results

        # Parse sample sheet and validate results exist
        patients = parse_ctdna_sample_sheet_for_existing_results(sample_sheet, results_dir)

        # Run analysis
        analyzer = CtdnaAnalyzer(family_size=family_size)
        output = output_dir or results_dir / "ctdna_analysis"
        output.mkdir(parents=True, exist_ok=True)

        result_file = analyzer.analyze_all(patients, results_dir, output)
        console.print(f"[green]Analysis complete![/green] Results: {result_file}")
```

#### 2.2 Add `--analysis-mode` to batch command

Add parameter to batch command:

```python
analysis_mode: Annotated[
    Optional[str],
    typer.Option("--analysis-mode", help="Run additional analysis after processing (e.g., 'ctdna').")
] = None,
analysis_family_size: Annotated[
    int,
    typer.Option("--analysis-family-size", help="Family size threshold for analysis.")
] = 3,
```

Add logic after `batch_process()` returns:

```python
# Run analysis if requested
if analysis_mode == "ctdna":
    from umierrorcorrect.ctdna import CtdnaAnalyzer
    from umierrorcorrect.ctdna.sample_sheet import parse_ctdna_sample_sheet

    if sample_sheet is None:
        console.print("[red]Error:[/red] --analysis-mode ctdna requires --sample-sheet with patient/mutation data")
        raise typer.Exit(1)

    # Re-parse sample sheet for ctDNA (gets patient/mutation info)
    patients = parse_ctdna_sample_sheet(sample_sheet, base_path=sample_sheet.parent)

    analyzer = CtdnaAnalyzer(family_size=analysis_family_size)
    analysis_output = output_dir / "ctdna_analysis"
    result_file = analyzer.analyze_all(patients, output_dir, analysis_output)

    console.print(f"[green]ctDNA analysis complete![/green] Results: {result_file}")
```

---

### Phase 3: Model Updates

#### 3.1 Update `models.py` - Make patient_id optional in SampleSheet

The current `SampleSheet._parse_csv()` assumes patient_id exists. Update to handle both cases:

```python
def _parse_csv(self) -> list[Patient]:
    """Parse the CSV file and return a list of Patient objects."""
    # ... existing code ...

    for row_num, row in enumerate(reader, start=2):
        patient_id = row.get("patient_id", "").strip()

        # If no patient_id, skip patient-level grouping
        # (samples will be orphaned - only valid for basic processing)
        if not patient_id:
            continue  # Or handle differently

        # ... rest of existing logic ...
```

#### 3.2 Add output tracking fields to Sample model (optional enhancement)

```python
class Sample(BaseModel):
    # ... existing fields ...

    # Output paths (populated after processing)
    output_dir: Path | None = None
    consensus_bam: Path | None = None
    consensus_tsv: Path | None = None  # .cons file
    vcf_file: Path | None = None
```

---

## Sample Sheet Format

### Basic (generic pipeline - unchanged)

```csv
sample_name,read1,read2
Sample1,/path/R1.fq.gz,/path/R2.fq.gz
```

### Extended (ctDNA analysis)

```csv
sample_name,read1,read2,patient_id,mutation_bed,sample_type,ml_plasma,collection_date,ng_input,region_bed
P001_T0,/path/R1.fq.gz,/path/R2.fq.gz,P001,/path/P001_mutations.bed,ctdna,2.5,2024-01-15,10.5,/path/regions.bed
P001_T1,/path/R1.fq.gz,/path/R2.fq.gz,P001,/path/P001_mutations.bed,ctdna,2.0,2024-03-20,8.0,/path/regions.bed
P001_VAL,/path/R1.fq.gz,/path/R2.fq.gz,P001,/path/P001_mutations.bed,validation,,,15.0,/path/regions.bed
```

---

## Output Format

### ctdna_analysis_results.tsv

```text
patient_id	sample_name	sample_type	collection_date	mutation_name	chromosome	position	ref	alt	vaf	mm_count	total_reads	mm_per_ml	ml_plasma	ng_input	family_size
P001	P001_T0	ctdna	2024-01-15	KRAS_G12D	chr12	25398284	G	A	0.0234	12	513	4.8	2.5	10.5	3
P001	P001_T0	ctdna	2024-01-15	TP53_R248W	chr17	7577538	C	T	0.0089	5	562	2.0	2.5	10.5	3
P001	P001_T1	ctdna	2024-03-20	KRAS_G12D	chr12	25398284	G	A	0.0156	8	513	4.0	2.0	8.0	3
```

---

## Usage Examples

### Standalone analysis on existing results

```bash
# Process samples first (generic pipeline)
umierrorcorrect batch --sample-sheet samples.csv -r genome.fa -o results/

# Run ctDNA analysis separately
umierrorcorrect analysis --ctdna --sample-sheet samples.csv --results-dir results/
```

### Integrated ctDNA analysis

```bash
# Process and analyze in one command
umierrorcorrect batch --sample-sheet samples.csv -r genome.fa -o results/ --analysis-mode ctdna
```

### Custom family size

```bash
umierrorcorrect analysis --ctdna -s samples.csv -r results/ --family-size 5
```

---

## Testing Plan

1. **Unit tests for MutationTracker**
   - Parse .cons file correctly
   - Handle missing positions gracefully
   - Calculate VAF correctly
   - Handle different family sizes

2. **Unit tests for sample sheet parsing**
   - Basic format (backward compatibility)
   - Extended format with all fields
   - Missing optional fields
   - Invalid data handling

3. **Integration tests**
   - Run analysis on test data output
   - Verify TSV output format
   - Test --analysis-mode ctdna in batch

4. **Test with existing test_data**
   - Use existing test FASTQ files
   - Create mock mutation BED file
   - Verify end-to-end workflow

---

## Implementation Order

1. **Phase 1.2**: `mutation_tracker.py` - Core parsing logic
2. **Phase 1.5**: `export.py` - TSV export
3. **Phase 1.4**: `analyzer.py` - Orchestration
4. **Phase 1.3**: `sample_sheet.py` - Extended parsing
5. **Phase 2.1**: CLI `analysis` command
6. **Phase 2.2**: CLI `batch --analysis-mode`
7. **Phase 3**: Model updates (if needed)

---

## Open Questions (Resolved)

| Question | Resolution |
| ---------- | ---------- |
| CLI structure | `analysis --ctdna` command + `batch --analysis-mode ctdna` |
| Family size | Configurable, default 3 |
| Report format | TSV for now, reporting module later |
| Integration | Both standalone and integrated modes supported |
