from augur.proximity import (
    get_valid_range,
    to_numpy_array,
)


class TestGetValidRange:
    def test_start_end(self):
        seq = to_numpy_array("NNNNNATGNN")
        #                     0123456789
        start,end = get_valid_range(seq)
        assert start == 5
        assert end == 8
