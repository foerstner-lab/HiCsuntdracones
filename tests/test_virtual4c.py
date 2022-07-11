import hicsuntdracones.virtual4c as v
import hashlib


def test_virtual4c():
    v4c = v.Virtual4C("tests/fixtures/50000_testmatrix.txt",
                      "tests/fixtures/Testgenome.gff",
                      65,
                      "x",
                      "myout_")

def test_read_gff_file():
    v4c = v.Virtual4C("tests/fixtures/50000_testmatrix.txt",
                      "tests/fixtures/Testgenome.gff",
                      65,
                      "x",
                      "myout_")

    v4c.read_gff_file()

    target_files_and_hashs = {}
    generated_files_and_hashs = {}

    # get_hash

    md5_read_gff_file_v4c_target = _get_hash("tests/fixtures/myout_target_virtual_4c.txt")
    md5_read_gff_file_v4c_generated = _get_hash("./myout___")

    # add to dictionary

    target_files_and_hashs["v4c"] = md5_read_gff_file_v4c_target
    generated_files_and_hashs["v4c"] = md5_read_gff_file_v4c_generated

    assert target_files_and_hashs == generated_files_and_hashs

def _get_hash(path):
    hash = hashlib.md5(open(path, 'rb').read()).hexdigest()
    return hash

    v4c.read_matrix_file_and_add_as_features()
    v4c.extract_features_overlapping_bins()
    v4c.generate_wiggle_file()
