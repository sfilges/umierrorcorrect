"""Unit tests for umierrorcorrect.get_consensus_statistics module."""

from umierrorcorrect.core.constants import DEFAULT_FAMILY_SIZES
from umierrorcorrect.get_consensus_statistics import (
    RegionConsensusStats,
    calculate_target_coverage,
    get_overall_statistics,
)


class TestRegionConsensusStats:
    """Tests for RegionConsensusStats class."""

    def test_init_basic(self):
        """Test basic initialization."""
        fsizes = [1, 2, 3, 5, 10]
        stat = RegionConsensusStats("1", "chr1:100-200", "gene1", 5, fsizes)

        assert stat.regionid == "1"
        assert stat.pos == "chr1:100-200"
        assert stat.name == "gene1"
        assert stat.singletons == 5
        assert stat.family_sizes == []
        assert stat.fsizes == fsizes

    def test_init_singleton_counts(self):
        """Test that singletons are counted correctly at initialization."""
        fsizes = [1, 2, 3, 5, 10]
        stat = RegionConsensusStats("1", "chr1:100-200", "gene1", 5, fsizes)

        # Singletons should contribute to threshold 0 and 1
        assert stat.total_reads[0] == 5
        assert stat.umis[0] == 5
        assert stat.total_reads[1] == 5
        assert stat.umis[1] == 5

        # Higher thresholds should start at 0
        assert stat.total_reads[2] == 0
        assert stat.umis[2] == 0

    def test_init_no_singletons(self):
        """Test initialization with no singletons."""
        fsizes = [1, 2, 3]
        stat = RegionConsensusStats("1", "chr1:100-200", "", 0, fsizes)

        assert stat.singletons == 0
        assert stat.total_reads[0] == 0
        assert stat.umis[0] == 0

    def test_add_family_sizes_basic(self):
        """Test adding family sizes updates statistics correctly."""
        fsizes = [1, 2, 3, 5]
        stat = RegionConsensusStats("1", "chr1:100-200", "", 0, fsizes)

        # Add families of sizes [5, 3, 2, 2, 1]
        sizes = [5, 3, 2, 2, 1]
        stat.add_family_sizes(sizes, fsizes)

        # Total reads at threshold 0: 5+3+2+2+1 = 13
        assert stat.total_reads[0] == 13
        # UMIs at threshold 0: Should equal total reads (13) because fsize=0 implies "raw reads"
        assert stat.umis[0] == 13

        # At threshold 1: all families pass
        assert stat.total_reads[1] == 13
        assert stat.umis[1] == 5

        # At threshold 2: families [5,3,2,2] pass = 12 reads, 4 families
        assert stat.total_reads[2] == 12
        assert stat.umis[2] == 4

        # At threshold 3: families [5,3] pass = 8 reads, 2 families
        assert stat.total_reads[3] == 8
        assert stat.umis[3] == 2

        # At threshold 5: families [5] pass = 5 reads, 1 family
        assert stat.total_reads[5] == 5
        assert stat.umis[5] == 1

    def test_add_family_sizes_extends_list(self):
        """Test that family sizes are added to the list."""
        fsizes = [1, 2, 3]
        stat = RegionConsensusStats("1", "chr1:100-200", "", 0, fsizes)

        stat.add_family_sizes([5, 3], fsizes)
        assert stat.family_sizes == [5, 3]

        stat.add_family_sizes([2, 1], fsizes)
        assert stat.family_sizes == [5, 3, 2, 1]

    def test_add_family_sizes_multiple_calls(self):
        """Test that multiple add_family_sizes calls accumulate correctly."""
        fsizes = [1, 2, 3]
        stat = RegionConsensusStats("1", "chr1:100-200", "", 0, fsizes)

        stat.add_family_sizes([5, 3], fsizes)
        stat.add_family_sizes([2, 1], fsizes)

        # Total: 5+3+2+1 = 11 reads
        # UMIs at threshold 0 should equal total reads (11)
        assert stat.total_reads[0] == 11
        assert stat.umis[0] == 11

    def test_write_stats_format(self):
        """Test that write_stats produces correct format."""
        fsizes = [1, 2, 3]
        stat = RegionConsensusStats("1", "chr1:100-200", "gene1", 2, fsizes)
        stat.add_family_sizes([5, 3], fsizes)

        output = stat.write_stats()
        lines = output.split("\n")

        # Should have one line per threshold: 0, 1, 2, 3
        assert len(lines) == 4

        # Check first line (threshold 0)
        parts = lines[0].split("\t")
        assert parts[0] == "1"  # regionid
        assert parts[1] == "chr1:100-200"  # pos
        assert parts[2] == "gene1"  # name
        assert parts[3] == "0"  # threshold


class TestCalculateTargetCoverage:
    """Tests for calculate_target_coverage function."""

    def test_does_not_mutate_input(self):
        """Test that the function does not mutate the input fsizes list."""
        fsizes = list(DEFAULT_FAMILY_SIZES)[1:]
        original_fsizes = fsizes.copy()

        stat = RegionConsensusStats("1", "chr1:100-200", "gene1", 5, fsizes)
        stat.add_family_sizes([10, 5, 3], fsizes)

        calculate_target_coverage([stat])

        # fsizes should be unchanged
        assert fsizes == original_fsizes

    def test_named_vs_unnamed_regions(self):
        """Test that only named regions contribute to target coverage."""
        fsizes = list(DEFAULT_FAMILY_SIZES)[1:]

        # Named region
        named = RegionConsensusStats("1", "chr1:100-200", "gene1", 0, fsizes)
        named.add_family_sizes([10], fsizes)

        # Unnamed region
        unnamed = RegionConsensusStats("2", "chr1:300-400", "", 0, fsizes)
        unnamed.add_family_sizes([10], fsizes)

        output = calculate_target_coverage([named, unnamed])
        lines = output.split("\n")

        # First line is threshold 0
        parts = lines[0].split("\t")
        target_reads = int(parts[1])
        all_reads = int(parts[2])

        # At threshold 0, target_reads should be total reads in named regions (10)
        # and all_reads should be total reads in all regions (20)
        assert target_reads == 10
        assert all_reads == 20

    def test_empty_stats(self):
        """Test with empty statistics list."""
        output = calculate_target_coverage([])
        lines = output.split("\n")

        # Should still produce output for all thresholds
        assert len(lines) == len(DEFAULT_FAMILY_SIZES)


class TestGetOverallStatistics:
    """Tests for get_overall_statistics function."""

    def test_aggregates_across_regions(self):
        """Test that statistics are correctly aggregated across regions."""
        fsizes = list(DEFAULT_FAMILY_SIZES)[1:]

        stat1 = RegionConsensusStats("1", "chr1:100-200", "gene1", 2, fsizes)
        stat1.add_family_sizes([5, 3], fsizes)

        stat2 = RegionConsensusStats("2", "chr1:300-400", "gene2", 3, fsizes)
        stat2.add_family_sizes([4, 2], fsizes)

        overall = get_overall_statistics([stat1, stat2])

        # Total reads at threshold 0:
        # stat1: 5+3 + 2 singletons = 10
        # stat2: 4+2 + 3 singletons = 9
        # overall: 19
        assert overall.total_reads[0] == 19

        # UMIs at threshold 0:
        # stat1: 5+3 + 2 singletons = 10
        # stat2: 4+2 + 3 singletons = 9
        # overall: 19
        assert overall.umis[0] == 19

    def test_overall_regionid(self):
        """Test that overall statistics have correct regionid."""
        fsizes = list(DEFAULT_FAMILY_SIZES)[1:]
        stat = RegionConsensusStats("1", "chr1:100-200", "gene1", 0, fsizes)

        overall = get_overall_statistics([stat])

        assert overall.regionid == "All"
        assert overall.pos == "all_regions"
        assert overall.name == ""

    def test_empty_regions(self):
        """Test with empty region list."""
        overall = get_overall_statistics([])

        assert overall.total_reads[0] == 0
        assert overall.umis[0] == 0
