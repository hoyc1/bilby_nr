import numpy as np
from bilby.core.utils import create_frequency_series
import pytest


class TestMultiModelModel(object):
    def setup_method(self):
        self.parameters = dict(
            mass_1 = 100,
            mass_2 = 50,
            luminosity_distance = 100,
            a_1 = 0.6,
            tilt_1 = np.pi / 3,
            phi_12 = np.pi / 2,
            a_2 = 0.2,
            tilt_2 = np.pi / 10,
            phi_jl = np.pi,
            theta_jn = np.pi / 3,
            phase = 0.,
        )
        self.waveform_kwargs = dict(
            waveform_approximant="IMRPhenomPv2",
            reference_frequency=50.0,
            minimum_frequency=20.0,
            catch_waveform_errors=True,
            match_interpolant="bilby_nr.match.match_from_interpolant",
            waveform_approximant_list=["IMRPhenomXPHMST", "IMRPhenomTPHM"]
        )
        self.frequency_array = create_frequency_series(2048, 4)

    def test_multi_model_binary_black_hole_no_interp(self):
        from bilby_nr.source import multi_model_binary_black_hole
        _wvf_args = self.waveform_kwargs.copy()
        _wvf_args.pop("match_interpolant")
        self.parameters.update(_wvf_args)
        pols = multi_model_binary_black_hole(
            self.frequency_array, **self.parameters
        )
        assert isinstance(pols, dict)

    def test_multi_model_binary_black_hole_use_best(self):
        from bilby_nr.source import multi_model_binary_black_hole
        for opt in ["True", True, 1]:
            _wvf_args = self.waveform_kwargs.copy()
            _wvf_args["use_best_match"] = opt
            self.parameters.update(_wvf_args)
            pols = multi_model_binary_black_hole(
                self.frequency_array, **self.parameters
            )
            assert isinstance(pols, dict)

    def test_multi_model_binary_black_hole_interp(self):
        from bilby_nr.source import multi_model_binary_black_hole
        self.parameters.update(self.waveform_kwargs)
        pols = multi_model_binary_black_hole(
            self.frequency_array, **self.parameters
        )
        assert isinstance(pols, dict)

    def test_multi_model_binary_black_hole_fake_interp(self):
        from bilby_nr.source import multi_model_binary_black_hole
        _wvf_args = self.waveform_kwargs.copy()
        _wvf_args["match_interpolant"] = "fake.fake_interpolant"
        with pytest.raises(ValueError):
            pols = multi_model_binary_black_hole(
                self.frequency_array, **self.parameters
            )
        _wvf_args = self.waveform_kwargs.copy()
        _wvf_args["match_interpolant"] = "bilby_nr.match.no_match_interp"
        with pytest.raises(ValueError):
            self.parameters.update(_wvf_args)
            pols = multi_model_binary_black_hole(
                self.frequency_array, **self.parameters
            )

    def test_multi_model_binary_black_hole_incorrect_weights(self):
        from bilby_nr.source import _multi_model_binary_black_hole
        with pytest.raises(ValueError):
            pols = _multi_model_binary_black_hole(
                [None, None], self.waveform_kwargs["waveform_approximant_list"],
                self.frequency_array, **self.parameters
            )
        with pytest.raises(ValueError):
            pols = _multi_model_binary_black_hole(
                [1.], self.waveform_kwargs["waveform_approximant_list"],
                self.frequency_array, **self.parameters
            )


def test_weights_from_matches():
    from bilby_nr.source import _weights_from_matches
    matches = np.array([0.3, 0.9, 0.8])
    for opt in ["True", True, 1]:
        weights = _weights_from_matches(matches, use_best=opt)
        np.testing.assert_almost_equal(weights, [0., 1., 0.])

    _true = 1 / ((1 - matches)**4)
    _true /= np.sum(_true)
    for opt in ["False", False, 0]:
        weights = _weights_from_matches(matches, use_best=False)
        np.testing.assert_almost_equal(weights, _true)
    
