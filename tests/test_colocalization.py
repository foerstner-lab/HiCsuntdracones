import hicsuntdracones.colocalization as c

def test_colo():
    obj = c.ColocalizationTester("tests/fixtures/50000_testmatrix.txt",
                                                              "tests/fixtures/Testgenome.gff",
                                                              "CDS",
                                                              65,
                                                              "myout_",
                                                              2,
                                                              10,
                                                              True,
                                                              True)
    obj.read_gff_file()
    obj.generate_interaction_matrix()
    obj.extract_features_overlapping_bins()
    obj.extract_random_bins()
    obj.compile_interaction_valus_of_bins()
    obj.perform_t_test()
    obj.write_countings_for_file()
    obj.plot_distribution()
