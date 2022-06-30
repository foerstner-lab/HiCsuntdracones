
import hicsuntdracones.colocalization as c
import hashlib
import os
import random

def test_colo():
    colo = c.ColocalizationTester("tests/fixtures/50000_testmatrix.txt",
                                                              "tests/fixtures/Testgenome.gff",
                                                              "CDS",
                                                              65,
                                                              "myout_",
                                                              2,
                                                              10,
                                                              True,
                                                              True)


def test_read_gff_file():
    random.seed(1)
    colo = c.ColocalizationTester("tests/fixtures/50000_testmatrix.txt",
                                  "tests/fixtures/Testgenome.gff",
                                  "CDS",
                                  65,
                                  "myout_",
                                  2,
                                  10,
                                  True,
                                  True)
    colo.read_gff_file()
    colo.generate_interaction_matrix()
    colo.extract_features_overlapping_bins()
    colo.extract_random_bins()
    colo.compile_interaction_valus_of_bins()
    colo.perform_t_test()
    colo.write_countings_for_file()
    colo.plot_distribution()

    #Creating Hashs for 2 out of 3 files (countings, t-test)

    md5_countings_reference = hashlib.md5(open("tests/fixtures/myout__countings.txt", 'rb').read()).hexdigest()
    md5_countings_generated = hashlib.md5(open("./myout__countings.txt", 'rb').read()).hexdigest()

    md5_ttest_reference = hashlib.md5(open("tests/fixtures/myout__t-test_results.txt", 'rb').read()).hexdigest()
    md5_ttest_generated = hashlib.md5(open("./myout__t-test_results.txt", 'rb').read()).hexdigest()

    #myout_histogram has a random base and cannot be tests with hashs
    #testing histogram for file size and type

    histo_output = "./myout__histograms.pdf"
    filesize_histo = os.path.getsize(histo_output)

    if histo_output.endswith('.pdf'):
        assert True
    else:
        assert False

    assert md5_countings_reference == md5_countings_generated
    assert md5_ttest_reference == md5_ttest_generated
    assert filesize_histo == 17173

    '''
    if os.path.exists("./myout__countings.txt"):
        os.remove("./myout__countings.txt")

    if os.path.exists("./myout__t-test_results.txt"):
        os.remove("./myout__t-test_results.txt")

    if os.path.exists("./myout__histograms.pdf"):
        os.remove("./myout__histograms.pdf")
    '''
