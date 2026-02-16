import numpy as np


class TestMatchInterpolants(object):
    def setup_method(self):
        self.parameters = dict(
            mass_1 = 100,
            mass_2 = 50,
            a_1 = 0.6,
            tilt_1 = np.pi / 3,
            phi_12 = np.pi / 2,
            a_2 = 0.2,
            tilt_2 = np.pi / 10,
            phi_jl = np.pi,
            theta_jn = np.pi / 3,
            phase = 0.,
        )
        self.waveform_approximant = "IMRPhenomXPHMST"

    def test_pade_pade_match_interpolant_for_float_input(self):
        from bilby_nr.match import match_from_interpolant

        match = match_from_interpolant(
            self.waveform_approximant,
            self.parameters["mass_1"],
            self.parameters["mass_2"],
            self.parameters["a_1"],
            self.parameters["tilt_1"],
            self.parameters["phi_12"],
            self.parameters["a_2"],
            self.parameters["tilt_2"],
            self.parameters["phi_jl"],
            self.parameters["theta_jn"],
            self.parameters["phase"],
            interp="pade_pade"
        )
        np.testing.assert_almost_equal(match, 0.9926632521327353)

    def test_pade_pade_match_interpolant_for_array_input(self):
        from bilby_nr.match import match_from_interpolant

        _parameters = {key: [item] for key, item in self.parameters.items()}
        match = match_from_interpolant(
            self.waveform_approximant,
            _parameters["mass_1"],
            _parameters["mass_2"],
            _parameters["a_1"],
            _parameters["tilt_1"],
            _parameters["phi_12"],
            _parameters["a_2"],
            _parameters["tilt_2"],
            _parameters["phi_jl"],
            _parameters["theta_jn"],
            _parameters["phase"],
            interp="pade_pade"
        )
        np.testing.assert_almost_equal(match, [0.9926632521327353])
