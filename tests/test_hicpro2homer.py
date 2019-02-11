import pandas as pd
import hicsuntdracones.hicpro2homer
import os


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
