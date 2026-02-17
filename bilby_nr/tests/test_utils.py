from bilby_nr import utils
import pytest

def test_convert_waveform_input_ValueError():
    with pytest.raises(ValueError):
        utils.convert_waveform_input(None)


def test_convert_waveform_input_list():
    out = utils.convert_waveform_input(["A", "B", "C"])
    assert out == ["A", "B", "C"]


def test_convert_waveform_input_string():
    variation = [
        "[A, B, C]",
        "['A', 'B', 'C']",
        "'A', 'B', 'C'",
        "A, B, C",
        "A,B,C",
    ]
    for opt in variation:
        out = utils.convert_waveform_input(opt)
        assert out == ["A", "B", "C"]
