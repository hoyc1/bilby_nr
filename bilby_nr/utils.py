# Licensed under an MIT style license -- see LICENSE.md

from bilby_pipe.utils import strip_quotes

__author__ = ["Charlie Hoy <charlie.hoy@port.ac.uk>"]


def convert_waveform_input(string):
    """Convert string inputs into a standard form for the waveform list

    Parameters
    ----------
    string: str
        A string representation to be converted
    """
    if string is None:
        raise VaueError("No input provided")
    if isinstance(string, list):
        string = ",".join(string)
    # Remove square brackets
    string = string.replace("[", "").replace("]", "")
    # Remove added quotes
    string = strip_quotes(string)
    # Replace multiple spaces with a single space
    string = " ".join(string.split())
    # Spaces can be either space or comma in input, convert to comma
    string = string.replace(" ,", ",").replace(", ", ",").replace(" ", ",")
    waveforms = string.split(",")
    waveforms.sort()
    return waveforms


def convert_waveform_list_from_input(waveform_list):
    """Convert a waveform list from the input into a list of strings

    Parameters
    ----------
    waveform_list: str or list
        The input waveform list

    Returns
    -------
    waveform_list: list
        The converted waveform list
    """
    if isinstance(waveform_list, list):
        return waveform_list
    return convert_waveform_input(waveform_list)
