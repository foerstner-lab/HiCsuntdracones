import hicsuntdracones.hicmatrix


def test_matrix_generation():
    hic_matrix = hicsuntdracones.hicmatrix.HiCMatrix(
        "tests/fixtures/homer_matrix.txt")
    assert type(hic_matrix) == hicsuntdracones.hicmatrix.HiCMatrix
    assert hic_matrix.hic_matrix_df.shape == (12, 14)
    assert hic_matrix.number_of_bins == 12
