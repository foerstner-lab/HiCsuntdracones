import hicsuntdracones.distdepdecay as d
import hashlib
import shutil
import os
'''

def test_distdepdecay():
    dddo_generator = d.DistDepDecayOutputGenerator("tests/fixtures/50000_testmatrix.txt",
                                        65,
                                        "myout_")
    #dddo_generator.write_table_file()
    #dddo_generator.plot_bin_averages()
    #dddo_generator.plot_bin_averages_with_error_bars()
    #dddo_generator.plot_all_values()

    table_file_content_target = open("tests/fixtures/ddd_example.csv").read()
    table_file_content_generated = open("./myout__distance_dependent_decay.csv").read()


    #md5_code_target = hashlib.md5(str("tests/fixtures/ddd_example.csv").encode('utf-8')).hexdigest()
    #md5_code_generated = hashlib.md5(str("./myout__distance_dependent_decay.csv").encode('utf-8')).hexdigest()
    md5_code_target1 = hashlib.md5(open("tests/fixtures/ddd_example.csv", 'rb').read()).hexdigest()
    md5_code_target2 = hashlib.md5(open("./myout__distance_dependent_decay.csv", 'rb').read()).hexdigest()
    assert table_file_content_target == table_file_content_generated
    assert md5_code_target1 == md5_code_target2
    # md5_code_target = hashlib.md5("tests/fixtures/ddd_example.csv").encode('utf-8').hexdigest()
'''


def test_write_table():
    dddo_generator = d.DistDepDecayOutputGenerator("tests/fixtures/50000_testmatrix.txt",
                                                   65,
                                                   "myout_")
    dddo_generator.write_table_file()
    write_table_file_content_target = open("tests/fixtures/ddd_example.csv").read()
    write_table_file_content_generated = open("./myout__distance_dependent_decay.csv").read()


    md5_writetable_code_target = hashlib.md5(open("tests/fixtures/ddd_example.csv", 'rb').read()).hexdigest()
    md5_writetable_code_generated = hashlib.md5(open("./myout__distance_dependent_decay.csv", 'rb').read()).hexdigest()

    assert write_table_file_content_target == write_table_file_content_generated
    assert md5_writetable_code_target == md5_writetable_code_generated


def test_plot_bin_averages():
    dddo_generator = d.DistDepDecayOutputGenerator("tests/fixtures/50000_testmatrix.txt",
                                                   65,
                                                   "myout_")
    dddo_generator.plot_bin_averages()
    plot_bin_target = open("tests/fixtures/ddd_example.csv").read()
    plot_bin_generated = open("./myout__distance_dependent_decay.csv").read()


    md5_plot_bin_target = hashlib.md5(open("tests/fixtures/ddd_example.csv", 'rb').read()).hexdigest()
    md5_plot_bin_generated= hashlib.md5(open("./myout__distance_dependent_decay.csv", 'rb').read()).hexdigest()

    assert plot_bin_target == plot_bin_generated
    assert md5_plot_bin_target == md5_plot_bin_generated


def test_plot_bin_averages_with_error_bars():
    dddo_generator = d.DistDepDecayOutputGenerator("tests/fixtures/50000_testmatrix.txt",
                                                   65,
                                                   "myout_")
    dddo_generator.plot_bin_averages_with_error_bars()
    plot_bin_averages_with_error_bars_target = open("tests/fixtures/ddd_example.csv").read()
    plot_bin_averages_with_error_bars_generated = open("./myout__distance_dependent_decay.csv").read()


    md5_plot_bin_averages_with_error_bars_target = hashlib.md5(open("tests/fixtures/ddd_example.csv", 'rb').read()).hexdigest()
    md5_plot_bin_averages_with_error_bars_generated = hashlib.md5(open("./myout__distance_dependent_decay.csv", 'rb').read()).hexdigest()

    assert plot_bin_averages_with_error_bars_target == plot_bin_averages_with_error_bars_generated
    assert md5_plot_bin_averages_with_error_bars_target == md5_plot_bin_averages_with_error_bars_generated


def test_plot_all_values():
    dddo_generator = d.DistDepDecayOutputGenerator("tests/fixtures/50000_testmatrix.txt",
                                                   65,
                                                   "myout_")
    dddo_generator.plot_all_values()
    plot_all_values_target = open("tests/fixtures/ddd_example.csv").read()
    plot_all_values_generated = open("./myout__distance_dependent_decay.csv").read()


    md5_plot_all_values_target = hashlib.md5(open("tests/fixtures/ddd_example.csv", 'rb').read()).hexdigest()
    md5_plot_all_values_generated = hashlib.md5(open("./myout__distance_dependent_decay.csv", 'rb').read()).hexdigest()

    assert plot_all_values_target == plot_all_values_generated
    assert md5_plot_all_values_target == md5_plot_all_values_generated

    if os.path.exists("./myout__distance_dependent_decay.csv"):
        os.remove("./myout__distance_dependent_decay.csv")

#def _get_hash(path):
 #   hash = hashlib.md5(open(path, 'rb').read()).hexdigest()
  #  return hash

