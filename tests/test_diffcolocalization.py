#import pytest
import hicsuntdracones.diffcolocalization as d

#@pytest.mark.skip(reason="not a proper test needs to be fixed")
def test_colo():
    diffcolo = d.DiffColocalizationTester("tests/fixtures/50000_testmatrix.txt",
                                 "tests/fixtures/50000_testmatrix.txt",
                                 "tests/fixtures/Testgenome.gff",
                                 "CDS",
                                 65,
                                 "myout_",
                                 2,
                                 10,
                                 True,
                                 True)

    #diffcolo.read_gff_file()



def test_read_gff_file():
    diffcolo = d.DiffColocalizationTester("tests/fixtures/50000_testmatrix.txt",
                                          "tests/fixtures/50000_testmatrix.txt",
                                          "tests/fixtures/Testgenome.gff",
                                          "CDS",
                                          65,
                                          "myout_",
                                          2,
                                          10,
                                          True,
                                          True)
    diffcolo.read_gff_file()

    diffcolo.generate_both_interaction_matrices()
    diffcolo.extract_features_overlapping_bins()
    diffcolo.extract_random_bins()
    diffcolo.compile_interaction_valus_of_bins()
    diffcolo.calculate_interaction_ratios()
    diffcolo.perform_t_test()
    diffcolo.write_countings_for_file()
    diffcolo.plot_distribution()
    '''
    plot_all_values_target = open("tests/fixtures/ddd_example.csv").read()
    plot_all_values_generated = open("./myout__distance_dependent_decay.csv").read()


    md5_plot_all_values_target = hashlib.md5(open("tests/fixtures/ddd_example.csv", 'rb').read()).hexdigest()
    md5_plot_all_values_generated = hashlib.md5(open("./myout__distance_dependent_decay.csv", 'rb').read()).hexdigest()

    assert plot_all_values_target == plot_all_values_generated
    assert md5_plot_all_values_target == md5_plot_all_values_generated
'''
