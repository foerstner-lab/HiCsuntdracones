import pandas as pd
import hicsuntdracones.hicpro2homer
import os
import hashlib

def test_hicpro2homer():
    tmp_output = "tests/fixtures/homer_matrix.txt"
    hicsuntdracones.hicpro2homer.hicpro2homer(
        "tests/fixtures/hicpro_abs.bed",
        "tests/fixtures/hicpro.matrix",
        tmp_output)
    output_homer_table = pd.read_csv(tmp_output, sep="\t")
    wanted_homer_table = pd.read_csv("tests/fixtures/homer_matrix.txt", sep="\t")
    pd.testing.assert_frame_equal(output_homer_table, wanted_homer_table)
    os.remove(tmp_output)


'''

    md5_outer_homer_table = hashlib.md5(open(tmp_output).read()).hexdigest()
    md5_wanted_homer_table = hashlib.md5(open("tests/fixtures/homer_matrix.txt", sep="\t").read()).hexdigest()
    assert md5_outer_homer_table == md5_wanted_homer_table
'''
