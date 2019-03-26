import pandas as pd
import numpy as np
import hicsuntdracones.hicmatrix


def test_matrix_generation():
    hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix.csv")
    assert type(hic_matrix) == hicsuntdracones.hicmatrix.HiCMatrix
    assert hic_matrix.hic_matrix_df.shape == (12, 14)
    assert hic_matrix.number_of_bins == 12
    assert hic_matrix.chromosomes == ["chr1", "chr2"]


def test_normalize_by_columns_sum():
    hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    assert hic_matrix._calc_column_sum_median() == 20.0
    hic_matrix.normalize_by_columns_sum()
    pd.testing.assert_frame_equal(
        pd.DataFrame.from_dict(
            {"HiCMatrix": [
                "chr1-0", "chr1-10000", "chr2-0", "chr2-10000"],
             "Regions": [
                 "chr1-0", "chr1-10000", "chr2-0", "chr2-10000"],
             "chr1-0": [0.50, 0.15, 0.0, 0.35],
             "chr1-10000": [0.15, 0.65, 0.0, 0.20],
             "chr2-0": [0.0, 0.0, 0.0, 0.0],
             "chr2-10000": [0.35, 0.20, 0.00, 0.45]}),
        hic_matrix.hic_matrix_df)


def test_bins():
    hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix.csv")
    assert type(hic_matrix.bins()) == pd.Series
    pd.testing.assert_series_equal(hic_matrix.bins(), pd.Series(
        ["chr1-0", "chr1-10000", "chr1-20000", "chr1-30000", "chr1-40000",
         "chr2-0", "chr2-10000", "chr2-20000", "chr2-30000", "chr2-40000",
         "chr2-50000", "chr2-60000"], name="bins"))


def test_select_keep_pattern():
    hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    sub_hic_matrix = hic_matrix.select(keep_pattern="chr1")
    pd.testing.assert_frame_equal(
        pd.DataFrame.from_dict(
            {"HiCMatrix": [
                "chr1-0", "chr1-10000"],
             "Regions": [
                 "chr1-0", "chr1-10000"],
             "chr1-0": [10.0, 3.0],
             "chr1-10000": [3.0, 13.0]}),
        sub_hic_matrix.hic_matrix_df)


def test_select_keep_pattern_inplace():
    hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    hic_matrix.select(keep_pattern="chr1", inplace=True)
    pd.testing.assert_frame_equal(
        pd.DataFrame.from_dict(
            {"HiCMatrix": [
                "chr1-0", "chr1-10000"],
             "Regions": [
                 "chr1-0", "chr1-10000"],
             "chr1-0": [10.0, 3.0],
             "chr1-10000": [3.0, 13.0]}),
        hic_matrix.hic_matrix_df)


def test_select_remove_pattern():
    hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    sub_hic_matrix = hic_matrix.select(remove_pattern="chr2")
    pd.testing.assert_frame_equal(
        pd.DataFrame.from_dict(
            {"HiCMatrix": [
                "chr1-0", "chr1-10000"],
             "Regions": [
                 "chr1-0", "chr1-10000"],
             "chr1-0": [10.0, 3.0],
             "chr1-10000": [3.0, 13.0]}),
        sub_hic_matrix.hic_matrix_df)


def test_select_keep_vs_remove():
    """Test if the remove pattern is stronger than the keep pattern.
    """
    hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    sub_hic_matrix = hic_matrix.select(
        keep_pattern="chr1", remove_pattern="chr1-0")
    pd.testing.assert_frame_equal(
        pd.DataFrame.from_dict(
            {"HiCMatrix": ["chr1-10000"],
             "Regions": ["chr1-10000"],
             "chr1-10000": [13.0]}),
        sub_hic_matrix.hic_matrix_df)


def test_matrix_values():
    hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    data = {"chr1-0": [10.0, 3.0, 0.0, 7.0],
             "chr1-10000": [3.0, 13.0, 0.0, 4.0],
             "chr2-0": [0.0, 0.0, 0.0, 0.0],
             "chr2-10000": [7.0, 4.0, 0.0, 9.0]}
    pd.testing.assert_frame_equal(
        pd.DataFrame.from_dict(data, orient='index',
                               columns=["chr1-0", "chr1-10000", "chr2-0", "chr2-10000"]),
        hic_matrix.matrix_values())


def test__div_matrix():
    hic_matrix_1 = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    hic_matrix_2 = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    div_matrix = hic_matrix_1._div_by(hic_matrix_2, 0)
    print(hic_matrix_1.matrix_values())
    print(hic_matrix_2.matrix_values())
    print(div_matrix)
    data = {"HiCMatrix": [
                "chr1-0", "chr1-10000", "chr2-0", "chr2-10000"],
             "Regions": [
                 "chr1-0", "chr1-10000", "chr2-0", "chr2-10000"],
             "chr1-0": [1.0, 1.0, np.nan, 1.0],
             "chr1-10000": [1.0, 1.0, np.nan, 1.0],
             "chr2-0": [np.nan, np.nan, np.nan, np.nan],
             "chr2-10000": [1.0, 1.0, np.nan, 1.0]}
    pd.testing.assert_frame_equal(
        pd.DataFrame.from_dict(data, orient='index',
                               columns=["chr1-0", "chr1-10000", "chr2-0", "chr2-10000"]),
        div_matrix)


def test__div_matrix_2():
    hic_matrix_1 = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    hic_matrix_2 = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    div_matrix = hic_matrix_1._div_by(hic_matrix_2, 1)
    data = {"HiCMatrix": [
                "chr1-0", "chr1-10000", "chr2-0", "chr2-10000"],
             "Regions": [
                 "chr1-0", "chr1-10000", "chr2-0", "chr2-10000"],
             "chr1-0": [1.0, 1.0, 1.0, 1.0],
             "chr1-10000": [1.0, 1.0, 1.0, 1.0],
             "chr2-0": [1.0, 1.0, 1.0, 1.0],
             "chr2-10000": [1.0, 1.0, 1.0, 1.0]}
    pd.testing.assert_frame_equal(
        pd.DataFrame.from_dict(data, orient='index',
                               columns=["chr1-0", "chr1-10000", "chr2-0", "chr2-10000"]),
        div_matrix)


def test_div_matrix():
    hic_matrix_1 = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    hic_matrix_2 = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    diff_hic_matrix = hic_matrix_1.div_by(hic_matrix_2, pseudocount=1.0)
    data = {"HiCMatrix": [
                "chr1-0", "chr1-10000", "chr2-0", "chr2-10000"],
             "Regions": [
                 "chr1-0", "chr1-10000", "chr2-0", "chr2-10000"],
             "chr1-0": [1.0, 1.0, 1.0, 1.0],
             "chr1-10000": [1.0, 1.0, 1.0, 1.0],
             "chr2-0": [1.0, 1.0, 1.0, 1.0],
             "chr2-10000": [1.0, 1.0, 1.0, 1.0]}
    pd.testing.assert_frame_equal(
        pd.DataFrame.from_dict(data, orient='index',
                               columns=["chr1-0", "chr1-10000", "chr2-0", "chr2-10000"]),
        diff_hic_matrix.hic_matrix_df)


def test_div_matrix_inplace():
    hic_matrix_1 = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    hic_matrix_2 = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix_small.txt")
    hic_matrix_1.div_by(hic_matrix_2, pseudocount=1.0, inplace=True)
    data = {"HiCMatrix": [
                "chr1-0", "chr1-10000", "chr2-0", "chr2-10000"],
             "Regions": [
                 "chr1-0", "chr1-10000", "chr2-0", "chr2-10000"],
             "chr1-0": [1.0, 1.0, 1.0, 1.0],
             "chr1-10000": [1.0, 1.0, 1.0, 1.0],
             "chr2-0": [1.0, 1.0, 1.0, 1.0],
             "chr2-10000": [1.0, 1.0, 1.0, 1.0]}
    pd.testing.assert_frame_equal(
        pd.DataFrame.from_dict(data, orient='index',
                               columns=["chr1-0", "chr1-10000", "chr2-0", "chr2-10000"]),
        hic_matrix_1.hic_matrix_df)


def test__get_counts_for_chroms():
    hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix.csv")
    chroms_dists_and_countings = hic_matrix.calc_distance_dependent_decay(
        10000)
    expected_chroms_dists_and_countings = {
        "chr1": {
            "countings": np.array(
                [54,  7, 96, 42, 11,
                 97, 48, 86, 84,
                 91, 95, 49,
                 57, 14,
                 70]),
            "dists": np.array([
                0, 10000, 20000, 30000, 40000,
                0, 10000, 20000, 30000,
                0, 10000, 20000,
                0, 10000,
                0])},
        "chr2" : {
            "countings": np.array(
                [87, 93, 93, 20, 61, 52, 18,
                 92, 94, 87, 42, 34, 97,
                 54, 11, 23, 20, 50,
                 49, 48, 5, 12,
                 81, 76, 13,
                 14, 85,
                 16]),
            "dists": np.array([
                0, 10000, 20000, 30000, 40000, 50000, 60000,
                0, 10000, 20000, 30000, 40000, 50000,
                0, 10000, 20000, 30000, 40000,
                0, 10000, 20000, 30000,
                0, 10000, 20000,
                0, 10000,
                0])}}
    np.testing.assert_array_equal(
        chroms_dists_and_countings["chr1"]["countings"],
        expected_chroms_dists_and_countings["chr1"]["countings"])
    np.testing.assert_array_equal(
        chroms_dists_and_countings["chr1"]["dists"],
        expected_chroms_dists_and_countings["chr1"]["dists"])
    np.testing.assert_array_equal(
        chroms_dists_and_countings["chr2"]["countings"],
        expected_chroms_dists_and_countings["chr2"]["countings"])
    np.testing.assert_array_equal(
        chroms_dists_and_countings["chr2"]["dists"],
        expected_chroms_dists_and_countings["chr2"]["dists"])
