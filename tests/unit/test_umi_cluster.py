"""Unit tests for umierrorcorrect.src.umi_cluster module."""

import pytest
from umierrorcorrect.src.umi_cluster import (
    cluster_barcodes,
    create_substring_matrix,
    get_connected_components,
    hamming_distance,
    merge_clusters,
    umi_cluster,
)


class TestHammingDistance:
    """Tests for hamming_distance function."""

    def test_identical_strings(self):
        """Identical strings should have distance 0."""
        assert hamming_distance("ACGT", "ACGT") == 0
        assert hamming_distance("AAAAAAAAAAAA", "AAAAAAAAAAAA") == 0

    def test_single_difference(self):
        """Strings differing by one character should have distance 1."""
        assert hamming_distance("ACGT", "ACGA") == 1
        assert hamming_distance("ACGT", "TCGT") == 1
        assert hamming_distance("AAAAAAAAAAAA", "AAAAAAAAAAAC") == 1

    def test_multiple_differences(self):
        """Test strings with multiple differences."""
        assert hamming_distance("ACGT", "TGCA") == 4
        assert hamming_distance("AAAA", "CCCC") == 4
        assert hamming_distance("ACGT", "ACCC") == 2

    def test_all_different(self):
        """Completely different strings should have max distance."""
        assert hamming_distance("AAAA", "TTTT") == 4

    def test_unequal_lengths_raises(self):
        """Strings of different lengths should raise AssertionError."""
        with pytest.raises(AssertionError):
            hamming_distance("ACGT", "ACGTA")

        with pytest.raises(AssertionError):
            hamming_distance("AC", "ACGT")

    def test_empty_strings(self):
        """Empty strings should have distance 0."""
        assert hamming_distance("", "") == 0


class TestUmiClusterClass:
    """Tests for umi_cluster class."""

    def test_init(self):
        """Test umi_cluster initialization."""
        umi = umi_cluster("ACGT", 10)
        assert umi.centroid == "ACGT"
        assert umi.count == 10

    def test_add_count(self):
        """Test adding counts to a cluster."""
        umi = umi_cluster("ACGT", 10)
        umi.add_count(5)
        assert umi.count == 15

    def test_change_centroid(self):
        """Test changing the centroid."""
        umi = umi_cluster("ACGT", 10)
        umi.change_centroid("TGCA")
        assert umi.centroid == "TGCA"


class TestCreateSubstringMatrix:
    """Tests for create_substring_matrix function."""

    def test_edit_distance_1(self):
        """Test substring matrix for edit distance threshold 1 (2 substrings)."""
        barcodes = {"ACGTACGT": 10, "ACGTACGA": 5}
        result = create_substring_matrix(barcodes, 1)

        assert len(result) == 2  # Two substring dictionaries
        substr_dict1, substr_dict2 = result

        # First half of ACGTACGT is ACGT
        assert "ACGT" in substr_dict1
        assert "ACGTACGT" in substr_dict1["ACGT"]

    def test_edit_distance_2(self):
        """Test substring matrix for edit distance threshold 2 (3 substrings)."""
        barcodes = {"ACGTACGTAC": 10, "ACGTACGTAG": 5}
        result = create_substring_matrix(barcodes, 2)

        assert len(result) == 3  # Three substring dictionaries


class TestClusterBarcodes:
    """Tests for cluster_barcodes function."""

    def test_no_similar_barcodes(self):
        """Barcodes with no similarity should not be clustered."""
        barcodes = {
            "AAAAAAAAAAAA": 10,
            "CCCCCCCCCCCC": 5,
        }
        adj_matrix = cluster_barcodes(barcodes, edit_distance_threshold=1)

        # No edges in adjacency matrix - barcodes are too different
        assert len(adj_matrix) == 0

    def test_similar_barcodes(self):
        """Barcodes within edit distance should be connected."""
        barcodes = {
            "AAAAAAAAAAAA": 10,
            "AAAAAAAAAAAC": 5,  # 1 edit distance
        }
        adj_matrix = cluster_barcodes(barcodes, edit_distance_threshold=1)

        # Higher count barcode should be the key, lower count should be in its list
        assert "AAAAAAAAAAAA" in adj_matrix
        assert "AAAAAAAAAAAC" in adj_matrix["AAAAAAAAAAAA"]

    def test_edit_distance_threshold(self):
        """Test that edit distance threshold is respected."""
        barcodes = {
            "AAAAAAAAAAAA": 10,
            "AAAAAAAAAACC": 5,  # 2 edit distance
        }

        # With threshold 1, they should not cluster
        adj_matrix_1 = cluster_barcodes(barcodes, edit_distance_threshold=1)
        assert len(adj_matrix_1) == 0

        # With threshold 2, they should cluster
        adj_matrix_2 = cluster_barcodes(barcodes, edit_distance_threshold=2)
        assert "AAAAAAAAAAAA" in adj_matrix_2

    def test_higher_count_is_centroid(self):
        """Barcode with higher count should always be the centroid."""
        barcodes = {
            "AAAAAAAAAAAA": 5,
            "AAAAAAAAAAAC": 10,  # Higher count
        }
        adj_matrix = cluster_barcodes(barcodes, edit_distance_threshold=1)

        # Higher count (AAAAAAAAAAAC) should be the key
        assert "AAAAAAAAAAAC" in adj_matrix
        assert "AAAAAAAAAAAA" in adj_matrix["AAAAAAAAAAAC"]

    def test_small_barcode_set_uses_combinations(self, sample_barcode_counts):
        """Small barcode sets (<30) use itertools.combinations."""
        # sample_barcode_counts has 5 barcodes, so should use combinations
        adj_matrix = cluster_barcodes(sample_barcode_counts, edit_distance_threshold=1)

        # ACGTACGTACGT (10) should cluster with ACGTACGTACGA (2) and ACGTACGTACGG (3)
        assert "ACGTACGTACGT" in adj_matrix
        assert "ACGTACGTACGA" in adj_matrix["ACGTACGTACGT"]
        assert "ACGTACGTACGG" in adj_matrix["ACGTACGTACGT"]

    def test_large_barcode_set_uses_substrings(self, sample_barcode_counts_large):
        """Large barcode sets (>30) use substring optimization."""
        adj_matrix = cluster_barcodes(sample_barcode_counts_large, edit_distance_threshold=1)

        # Known cluster: AAAAAAAAAAAA (50), AAAAAAAAAAAC (10), AAAAAAAAAAAG (5)
        assert "AAAAAAAAAAAA" in adj_matrix
        neighbors = adj_matrix["AAAAAAAAAAAA"]
        assert "AAAAAAAAAAAC" in neighbors or "AAAAAAAAAAAG" in neighbors


class TestGetConnectedComponents:
    """Tests for get_connected_components function."""

    def test_single_barcode(self):
        """Single barcode should form its own cluster."""
        barcodes = {"AAAAAAAAAAAA": 10}
        adj_matrix = {}
        clusters = get_connected_components(barcodes, adj_matrix)

        assert len(clusters) == 1
        assert clusters[0] == ["AAAAAAAAAAAA"]

    def test_connected_barcodes(self):
        """Connected barcodes should form a single cluster."""
        barcodes = {
            "AAAAAAAAAAAA": 10,
            "AAAAAAAAAAAC": 5,
        }
        adj_matrix = {"AAAAAAAAAAAA": ["AAAAAAAAAAAC"]}
        clusters = get_connected_components(barcodes, adj_matrix)

        assert len(clusters) == 1
        assert "AAAAAAAAAAAA" in clusters[0]
        assert "AAAAAAAAAAAC" in clusters[0]

    def test_disconnected_barcodes(self):
        """Disconnected barcodes should form separate clusters."""
        barcodes = {
            "AAAAAAAAAAAA": 10,
            "CCCCCCCCCCCC": 5,
        }
        adj_matrix = {}
        clusters = get_connected_components(barcodes, adj_matrix)

        assert len(clusters) == 2

    def test_highest_count_first_in_cluster(self):
        """Highest count barcode should be first (centroid) in cluster."""
        barcodes = {
            "AAAAAAAAAAAA": 10,  # Highest count
            "AAAAAAAAAAAC": 5,
            "AAAAAAAAAAAG": 3,
        }
        adj_matrix = {"AAAAAAAAAAAA": ["AAAAAAAAAAAC", "AAAAAAAAAAAG"]}
        clusters = get_connected_components(barcodes, adj_matrix)

        assert len(clusters) == 1
        # First element should be the highest count barcode
        assert clusters[0][0] == "AAAAAAAAAAAA"


class TestMergeClusters:
    """Tests for merge_clusters function."""

    def test_single_barcode_cluster(self):
        """Single barcode cluster should not change counts."""
        barcodes = {"AAAAAAAAAAAA": 10}
        clusters = [["AAAAAAAAAAAA"]]
        umis = merge_clusters(barcodes, clusters)

        assert umis["AAAAAAAAAAAA"].count == 10
        assert umis["AAAAAAAAAAAA"].centroid == "AAAAAAAAAAAA"

    def test_merged_cluster_counts(self):
        """Merged cluster should sum counts from all members."""
        barcodes = {
            "AAAAAAAAAAAA": 10,
            "AAAAAAAAAAAC": 5,
            "AAAAAAAAAAAG": 3,
        }
        clusters = [["AAAAAAAAAAAA", "AAAAAAAAAAAC", "AAAAAAAAAAAG"]]
        umis = merge_clusters(barcodes, clusters)

        # Centroid (first in cluster) should have sum of all counts
        assert umis["AAAAAAAAAAAA"].count == 18  # 10 + 5 + 3
        assert umis["AAAAAAAAAAAA"].centroid == "AAAAAAAAAAAA"

        # Other members should point to same umi_cluster object
        assert umis["AAAAAAAAAAAC"] is umis["AAAAAAAAAAAA"]
        assert umis["AAAAAAAAAAAG"] is umis["AAAAAAAAAAAA"]

    def test_multiple_clusters(self):
        """Multiple separate clusters should be handled correctly."""
        barcodes = {
            "AAAAAAAAAAAA": 10,
            "AAAAAAAAAAAC": 5,
            "CCCCCCCCCCCC": 8,
            "CCCCCCCCCCCT": 2,
        }
        clusters = [
            ["AAAAAAAAAAAA", "AAAAAAAAAAAC"],
            ["CCCCCCCCCCCC", "CCCCCCCCCCCT"],
        ]
        umis = merge_clusters(barcodes, clusters)

        assert umis["AAAAAAAAAAAA"].count == 15  # 10 + 5
        assert umis["CCCCCCCCCCCC"].count == 10  # 8 + 2

        # Different clusters should have different umi_cluster objects
        assert umis["AAAAAAAAAAAA"] is not umis["CCCCCCCCCCCC"]


class TestIntegration:
    """Integration tests for the full clustering pipeline."""

    def test_full_clustering_pipeline(self, sample_barcode_counts):
        """Test the full clustering pipeline from barcodes to merged clusters."""
        edit_distance_threshold = 1

        # Step 1: Create adjacency matrix
        adj_matrix = cluster_barcodes(sample_barcode_counts, edit_distance_threshold)

        # Step 2: Get connected components
        clusters = get_connected_components(sample_barcode_counts, adj_matrix)

        # Step 3: Merge clusters
        umis = merge_clusters(sample_barcode_counts, clusters)

        # Verify results
        # ACGTACGTACGT cluster should have: 10 + 2 + 3 = 15
        assert umis["ACGTACGTACGT"].count == 15
        assert umis["ACGTACGTACGA"].centroid == "ACGTACGTACGT"
        assert umis["ACGTACGTACGG"].centroid == "ACGTACGTACGT"

        # GGGGGGGGGGGG cluster should have: 5 + 1 = 6
        assert umis["GGGGGGGGGGGG"].count == 6
        assert umis["GGGGGGGGGGGA"].centroid == "GGGGGGGGGGGG"

    def test_empty_barcode_dict(self):
        """Empty barcode dictionary should return empty results."""
        barcodes = {}
        adj_matrix = cluster_barcodes(barcodes, edit_distance_threshold=1)
        assert adj_matrix == {}
