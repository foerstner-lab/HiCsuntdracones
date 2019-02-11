import hicsuntdracones.diffcolocalization as d


def test_colo():
    obj = d.DiffColocalizationTester("tests/fixtures/50000_testmatrix.txt",
                                 "tests/fixtures/50000_testmatrix.txt",
                                 "tests/fixtures/Testgenome.gff",
                                 "CDS",
                                 65,
                                 "myout_",
                                 2,
                                 10,
                                 True,
                                 True)

    obj.read_gff_file()
    obj.generate_both_interaction_matrices()
    obj.extract_features_overlapping_bins()
    obj.extract_random_bins()
    obj.compile_interaction_valus_of_bins()
    obj.calculate_interaction_ratios()
    obj.perform_t_test()
    obj.write_countings_for_file()
    obj.plot_distribution()
