from bilby_nr.bilby_pipe import Input
from bilby.gw.detector.networks import InterferometerList
from bilby_pipe.utils import BilbyPipeError
import os
import pytest
import tempfile

tmpdir = tempfile.TemporaryDirectory(prefix=".", dir=".").name


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


def TestMainInput(object):
    def setup_method(self):
        from bilby_nr.bilby_pipe import (
            create_parser, MainInput
        )
        from bilby_pipe.main import (
            create_parser as original_parser,
            MainInput as OriginalInput
        )
        from bilby_pipe.main import parse_args

        self.outdir = tmpdir
        self.default_args_list = [
            "--ini",
            os.path.join(os.path.dirname(__file__), "test_config.ini"),
            "--outdir",
            self.outdir,
        ]
        self.my_parser = create_parser()
        self.my_inputs = MainInput(
            *parse_args(self.default_args_list, self.my_parser), test=True
        )
        self.original_parser = original_parser()
        self.original_inputs = OriginalInput(
            *parse_args(self.default_args_list, self.original_parser),
            test=True
        )

    def teardown_method(self):
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)

    def test_original_behaviour(self):
        keys = vars(self.original_inputs).keys()
        for key in keys:
            if key in ["known_args"]:
                continue
            _orig = getattr(self.original_inputs, key)
            _new = getattr(self.my_inputs, key)
            assert _orig == _new

    def test_new_options(self):
        from bilby_nr.bilby_pipe import (
            create_parser, MainInput
        )
        from bilby_pipe.main import parse_args

        self.outdir = tmpdir
        self.default_args_list = [
            "--ini",
            os.path.join(
                os.path.dirname(__file__),
                "test_config_with_additional_arguments.ini"
            ),
            "--outdir",
            self.outdir,
        ]
        self.my_parser = create_parser()
        self.my_inputs = MainInput(
            *parse_args(self.default_args_list, self.my_parser), test=True
        )
        assert sorted(self.my_inputs.waveform_approximant) == sorted([
            "IMRPhenomXPHMST", "IMRPhenomTPHM", "SEOBNRv5PHM"
        ])
        wvf_args = self.my_inputs.get_default_waveform_arguments()
        assert sorted(wvf_args["waveform_approximant_list"]) == sorted([
            "IMRPhenomXPHMST", "IMRPhenomTPHM", "SEOBNRv5PHM"
        ])
        assert wvf_args["waveform_approximant"] == sorted([
            "IMRPhenomXPHMST", "IMRPhenomTPHM", "SEOBNRv5PHM"
        ])[0]


class TestDataAnalysisInput(object):
    def setup_method(self):
        from bilby_nr.bilby_pipe import (
            create_parser, DataAnalysisInput
        )
        from bilby_pipe.data_analysis import (
            create_analysis_parser as original_parser,
            DataAnalysisInput as OriginalInput
        )
        from bilby_pipe.main import parse_args

        self.outdir = tmpdir
        self.default_args_list = [
            "--ini",
            os.path.join(os.path.dirname(__file__), "test_config.ini"),
            "--outdir",
            self.outdir,
        ]
        self.my_parser = create_parser()
        self.my_inputs = DataAnalysisInput(
            *parse_args(self.default_args_list, self.my_parser), test=True
        )
        self.original_parser = original_parser()
        self.original_inputs = OriginalInput(
            *parse_args(self.default_args_list, self.original_parser),
            test=True
        )

    def teardown_method(self):
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir) 

    def test_original_behaviour(self):
        keys = vars(self.original_inputs).keys()
        for key in keys:
            if key in ["known_args"]:
                continue
            _orig = getattr(self.original_inputs, key)
            _new = getattr(self.my_inputs, key)
            assert _orig == _new

    def test_new_options(self):
        from bilby_nr.bilby_pipe import (
            create_parser, DataAnalysisInput
        )
        from bilby_pipe.main import parse_args

        self.outdir = tmpdir
        self.default_args_list = [
            "--ini",
            os.path.join(
                os.path.dirname(__file__),
                "test_config_with_additional_arguments.ini"
            ),
            "--outdir",
            self.outdir,
        ]
        self.my_parser = create_parser()
        self.my_inputs = DataAnalysisInput(
            *parse_args(self.default_args_list, self.my_parser), test=True
        )
        assert sorted(self.my_inputs.waveform_approximant) == sorted([
            "IMRPhenomXPHMST", "IMRPhenomTPHM", "SEOBNRv5PHM"
        ])
        wvf_args = self.my_inputs.get_default_waveform_arguments()
        assert sorted(wvf_args["waveform_approximant_list"]) == sorted([
            "IMRPhenomXPHMST", "IMRPhenomTPHM", "SEOBNRv5PHM"
        ])
        assert wvf_args["waveform_approximant"] == sorted([
            "IMRPhenomXPHMST", "IMRPhenomTPHM", "SEOBNRv5PHM"
        ])[0]
