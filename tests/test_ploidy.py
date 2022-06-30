import hicsuntdracones.ploidy as p
import pytest

@pytest.mark.skip(reason="not a proper test needs to be fixed")

def test_ploidy():
    obj = p.Ploidy("",
                   "",
                   "",
                   "",
                   "",
                   0)
    obj.main()
