# Licensed under an MIT style license -- see LICENSE.md

from bilby_pipe.input import Input as _Input
from bilby_pipe.main import MainInput as _MainInput, main
from bilby_pipe.data_analysis import DataAnalysisInput as _DAInput
from bilby_pipe.utils import logger, BilbyPipeError
import inspect

__author__ = ["Charlie Hoy <charlie.hoy@port.ac.uk>"]


class Input(_Input):
    """Superclass of input handlers inherited from bilby_pipe.input.Input
    """
    @property
    def waveform_approximant(self):
        return self._waveform_approximant

    @waveform_approximant.setter
    def waveform_approximant(self, waveform_approximant):
        """Set the waveform approximant.

        If a list of waveform approximants is provided, the code will check if
        the source model is compatible with sampling over multiple models. If it
        is, the code will set the `_sample_multiple_models` flag to True. If
        not, it will raise an error.

        Parameters
        ----------
        waveform_approximant: str or list
            The waveform approximant to use.
        """
        from .utils import convert_waveform_list_from_input
        self._sample_multiple_models = False
        self._waveform_approximant = waveform_approximant
        if waveform_approximant is not None:
            models = convert_waveform_list_from_input(waveform_approximant)
            if isinstance(models, list) and len(models) > 1:
                source_model = self.bilby_frequency_domain_source_model
                source_module = inspect.getmodule(source_model).__name__
                if source_module == "bilby_nr.source":
                    logger.warning(
                        f"A list has been provided for the waveform "
                        f"approximant. Sampling over the models: "
                        f"{','.join(models)}."
                    )
                    self._sample_multiple_models = True
                elif isinstance(models, list):
                    raise BilbyPipeError(
                        f"A waveform approximant list has been provided but "
                        f"a frequency domain source model inside the "
                        f"`bilby_nr` package has not been specified. Unable to "
                        f"sample over multiple models."
                    )
            elif isinstance(models, list):
                # recover the original behaviour
                models = models[0]
            self._waveform_approximant = models

    def get_default_waveform_arguments(self):
        """Return the default waveform arguments.

        If the code is sampling over multiple models, this method will return
        the arguments for the first model in the list, and add the full list
        to the dictionary under the key `waveform_approximant_list`.

        Returns
        -------
        wfa: dict
            The default waveform arguments.
        """
        wfa = super().get_default_waveform_arguments()
        if self._sample_multiple_models:
            # store just the first approximant to prevent issues downstream
            # and add a new waveform argument which is used by the
            # bilby_nr source models
            wfa["waveform_approximant"] = self.waveform_approximant[0]
            wfa["waveform_approximant_list"] = self.waveform_approximant
        return wfa


class MainInput(Input, _MainInput):
    """An object to hold all the inputs to bilby_pipe. Inherited from
    bilby_pipe.main.MainInput
    """
    def __init__(self, *args, **kwargs):
        self.frequency_domain_source_model = args[0].frequency_domain_source_model
        super().__init__(*args, **kwargs)


class DataAnalysisInput(Input, _DAInput):
    """Handles user-input for the data analysis script. Inhertied from
    bilby_pipe.data_analysis.DataAnalysisInput
    """
    def __init__(self, *args, **kwargs):
        self.frequency_domain_source_model = args[0].frequency_domain_source_model
        super().__init__(*args, **kwargs)


def create_parser(top_level=False):
    """Extends the BilbyArgParser for bilby_pipe to include additional
    options for sampling over multiple models

    Parameters
    ----------
    top_level:
        If true, parser is to be used at the top-level with requirement
        checking etc, else it is an internal call and will be ignored.
    """
    from bilby_pipe.parser import create_parser
    parser = create_parser(top_level=top_level)
    parser.set_defaults(
        main_input_class="bilby_nr.bilby_pipe.MainInput",
        analysis_input_class="bilby_nr.bilby_pipe.DataAnalysisInput",
    )
    return parser
