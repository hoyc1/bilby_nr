from bilby_nr.bilby_pipe import Input
from bilby.gw.detector.networks import InterferometerList
from bilby_pipe.utils import BilbyPipeError
import pytest


class TestInput(object):
    def test_original_waveform_approximant_behaviour(self):
        inputs = Input(None, None)
        inputs.waveform_approximant = "A"
        assert inputs.waveform_approximant == "A"

    def test_waveform_approximant_list_with_correct_source_model(self):
        variations = [
            ["A", "B"], "[A, B]", "['A', 'B']", "A, B"
        ]
        for opt in variations:
            inputs = Input(None, None)
            inputs.frequency_domain_source_model = (
                "bilby_nr.source.multi_model_binary_black_hole"
            )
            inputs.waveform_approximant = opt
            assert inputs.waveform_approximant == ["A", "B"]

    def test_waveform_approximant_list_with_incorrect_source_model(self):
        inputs = Input(None, None)
        inputs.frequency_domain_source_model = "lal_binary_black_hole"
        with pytest.raises(BilbyPipeError):
            inputs.waveform_approximant = ["A", "B"]

    def test_waveform_generator(self):
        for num, approx in enumerate(["A", "[A, B]"]):
            inputs = Input(None, None)
            inputs.detectors = ["H1"]
            inputs.interferometers = InterferometerList(["H1"])
            inputs.conversion_function = "noconvert"
            inputs.waveform_generator_class = "WaveformGenerator"
            inputs.detectors = ["H1"]
            inputs.waveform_generator_class_ctor_args = None
            inputs.likelihood_type = "GravitationalWaveTransient"
            if num == 0:
                inputs.frequency_domain_source_model = "lal_binary_black_hole"
            else:
                inputs.frequency_domain_source_model = (
                    "bilby_nr.source.multi_model_binary_black_hole"
                )
            inputs.waveform_approximant = approx
            inputs.catch_waveform_errors = False
            inputs.reference_frequency = 20
            inputs.minimum_frequency = 20
            inputs.maximum_frequency = 1024
            inputs.pn_spin_order = -1
            inputs.pn_tidal_order = -1
            inputs.pn_phase_order = -1
            inputs.pn_amplitude_order = 0
            inputs.mode_array = None
            inputs.waveform_arguments_dict = None
            args = inputs.waveform_generator.waveform_arguments
            assert args['waveform_approximant'] == "A"
            if num == 0:
                assert 'waveform_approximant_list' not in args.keys()
            else:
                assert args['waveform_approximant_list'] == ["A", "B"]
