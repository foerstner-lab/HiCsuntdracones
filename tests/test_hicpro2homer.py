import pandas as pd
import hicsuntdracones.hicpro2homer
import os


def test_hicpro2homer():
    tmp_output = "/tmp/homer_matrix.txt"
    hicsuntdracones.hicpro2homer.hicpro2homer(
        "tests/fixtures/hicpro_abs.bed",
        "tests/fixtures/hicpro.matrix",
        tmp_output)
    output_homer_table = pd.read_table(tmp_output)
    wanted_homer_table = pd.read_table("tests/fixtures/homer_matrix.txt")
    pd.testing.assert_frame_equal(output_homer_table, wanted_homer_table)
    os.remove(tmp_output)
