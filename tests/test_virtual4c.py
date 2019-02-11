import hicsuntdracones.virtual4c as v


def test_virtual4c():
    obj = v.Virtual4C("tests/fixtures/50000_testmatrix.txt",
                      "tests/fixtures/Testgenome.gff",
                      65,
                      "x",
                      "myout_")
    obj.read_gff_file()
    obj.read_matrix_file_and_add_as_features()
    obj.extract_features_overlapping_bins()
    obj.generate_wiggle_file()
