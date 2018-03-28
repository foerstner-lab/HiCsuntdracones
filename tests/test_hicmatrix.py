import pandas as pd
import hicsuntdracones.hicmatrix


def test_matrix_generation():
    hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix.txt")
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
