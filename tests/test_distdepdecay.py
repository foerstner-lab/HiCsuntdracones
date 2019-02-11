import hicsuntdracones.distdepdecay as d

def test_distdepdecay():
    obj = d.DistDepDecayOutputGenerator("tests/fixtures/50000_testmatrix.txt",
                                        65,
                                        "myout_")
    obj.write_table_file()
    obj.plot_bin_averages()
    obj.plot_bin_averages_with_error_bars()
    obj.plot_all_values()

