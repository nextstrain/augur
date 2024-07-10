from unittest import mock
from augur.export_v2 import is_numeric, format_number, attr_confidence

class TestIsNumeric():
    """This tests _our_ is_numeric function, it's not a philosophical discussion about numbers!"""
    def test_strings_are_not_numeric(self):
        assert not is_numeric("abc")
        assert not is_numeric("1")
        assert not is_numeric("1.0")
    def test_numbers_are_numeric(self):
        assert is_numeric(1)
        assert is_numeric(1.1)
        assert is_numeric(1.04e-2)

class TestFormatNumber():
    def test_integers_are_unchanged(self):
        assert 123 == format_number(123)
        assert 0 == format_number(0)
        assert -123 == format_number(-123)

    def test_positive_float_sig_figs(self):
        assert "2010.0" == str(format_number(2010.0))
        assert "2010.0" == str(format_number(2010.0000000001))
        assert "2010.123" == str(format_number(2010.1234000000001))
        assert "2010.0" == str(format_number(2010.0001))
        assert "0.123" == str(format_number(0.123456))

    def test_negative_float_sig_figs(self):
        assert "-2010.0" == str(format_number(-2010.0))
        assert "-2010.0" == str(format_number(-2010.0000000001))
        assert "-2010.123" == str(format_number(-2010.1234000000001))
        assert "-2010.0" == str(format_number(-2010.0001))
        assert "-0.123" == str(format_number(-0.123456))

    def test_exponential_notation_sig_figs(self):
        assert "-1e-12" == str(format_number(-1.000088900581841e-12))
        assert "-1.01e-12" == str(format_number(-1.0088900581841e-12))
        assert "-8.89e-15" == str(format_number(-0.0088900581841e-12))
        
        assert "1e-12" == str(format_number(1.000088900581841e-12))
        assert "1.01e-12" == str(format_number(1.0088900581841e-12))
        assert "8.89e-15" == str(format_number(0.0088900581841e-12))

    def test_exponential_notation_conversion(self):
        assert ("1.23e-10") == str(format_number(0.000_000_000_123_456))


def no_op(*args):
    pass

class TestConfidenceExtraction():
    def test_no_confidence(self):
        assert {} == attr_confidence("name", {"key": "something"}, "key")

    def test_no_confidence_but_entropy(self):
        assert {} == attr_confidence("name", {"key": "something", "key_entropy": 2.0}, "key")

    @mock.patch("augur.export_v2.warn", no_op)
    def test_confidence_but_no_entropy(self):
        # array confidences don't need entropy (and are represented internally as tuples)
        assert {'confidence': (1,2)} == attr_confidence("name", {"key": 1.5, "key_confidence": [1,2]}, "key")
        # but dict confidences do...
        assert {} == attr_confidence("name", {"key": "foo", "key_confidence": {"foo": 0.99}}, "key")

    @mock.patch("augur.export_v2.warn", no_op)
    def test_array_confidence(self):
        # Input JSON array -> python tuple -> exported as JSON Array
        assert {'confidence': (1,2)} == attr_confidence("name", {"key": 1.5, "key_confidence": [1,2]}, "key")
        # invalid confidence. Warnings suppressed.
        assert {} == attr_confidence("name", {"key": 1.5, "key_confidence": [1,2,3]}, "key")


    @mock.patch("augur.export_v2.warn", no_op)
    def test_dict_confidence(self):
        conf = {"x": 0.888, "y": 1.1e-2}
        entropy = 0.123
        assert {'confidence': conf, "entropy": entropy} == \
            attr_confidence("name", {"key": "x", "key_confidence": conf, "key_entropy": entropy}, "key")

        # values must all be numbers
        assert {} == \
            attr_confidence("name", {"key": "x", "key_confidence": {**conf, "c": "d"}, "key_entropy": entropy}, "key")    

    @mock.patch("augur.export_v2.warn", no_op)
    def test_unknown_confidence(self):
        assert {} == attr_confidence("name", {"key": "x", "key_confidence": "0.999"}, "key")