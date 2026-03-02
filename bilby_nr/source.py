# Licensed under an MIT style license -- see LICENSE.md

import numpy as np
from bilby.gw import source
import ast
from .utils import convert_waveform_list_from_input

__author__ = ["Charlie Hoy <charlie.hoy@port.ac.uk>"]


def multi_model_binary_black_hole(
    frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
):
    """Source model for a binary black hole with multiple models.

    Parameters
    ----------
    frequency_array: np.ndarray
        The frequency array
    mass_1: float
        The mass of the primary black hole
    mass_2: float
        The mass of the secondary black hole
    luminosity_distance: float
        The luminosity distance
    a_1: float
        The dimensionless spin magnitude of the primary black hole
    tilt_1: float
        The tilt angle of the primary black hole spin
    phi_12: float
        The difference in azimuthal angle between the two spins
    a_2: float
        The dimensionless spin magnitude of the secondary black hole
    tilt_2: float
        The tilt angle of the secondary black hole spin
    phi_jl: float
        The azimuthal angle of the total angular momentum
    theta_jn: float
        The angle between the total angular momentum and the line of sight
    phase: float
        The phase of the gravitational wave
    kwargs: dict
        Additional keyword arguments. We support additional arguments beyond
        those allowed in bilby. For example:

            - waveform_approximant_list: a list of waveform approximants
              you wish to sample over. If the user provides the waveform
              approximant 'IMRPhenomXPHMST', we populate the required flags
              to use the Spin-Taylor (ST) variant of IMRPhenomXPHM.
            - match_interpolant: the interpolant you wish to use to estimate the
              match for a given region in the parameter space
            - use_best_match: always use the model with the best match to
              evaluate the likelihood
            - match_to_weight: a string that can be evaluated to map an array
              of matches to a series of weights. The model will then be
              chosen probabilistically based on the weights.

    Returns
    -------
    polarizations: dict
        The polarizations
    """
    waveform_approximant_list = kwargs.pop("waveform_approximant_list", None)
    if waveform_approximant_list is None:
        raise ValueError(
            "Please provide a list of waveforms to sample over via the "
            "waveform_approximant_list waveform argument"
        )
    waveform_approximant_list = convert_waveform_list_from_input(
        waveform_approximant_list
    )
    match_interpolant = kwargs.pop("match_interpolant", None)
    if match_interpolant is not None:
        return _multi_model_match_informed_binary_black_hole(
            match_interpolant, waveform_approximant_list, frequency_array,
            mass_1, mass_2, luminosity_distance, a_1, tilt_1, phi_12, a_2,
            tilt_2, phi_jl, theta_jn, phase, **kwargs
        )
    return _multi_model_binary_black_hole(
        np.ones(len(waveform_approximant_list)) / len(waveform_approximant_list),
        waveform_approximant_list, frequency_array, mass_1, mass_2,
        luminosity_distance, a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn,
        phase, **kwargs
    )


def _multi_model_match_informed_binary_black_hole(
    interpolant, waveform_approximant_list, frequency_array, mass_1, mass_2,
    luminosity_distance, a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn,
    phase, **kwargs
):
    """Source model for a binary black hole with multiple models. The model
    used to generate the GW polarizations is chosen based on the match to
    numerical relativity simulations. The match is calculated based on the
    specified interpolant

    Parameters
    ----------
    interpolant: str
        The interpolant you wish to use to estimate the match
    waveform_approximant_list: list
        list of waveform approximants you wish to include
    frequency_array: np.ndarray
        The frequency array
    mass_1: float
        The mass of the primary black hole
    mass_2: float
        The mass of the secondary black hole
    luminosity_distance: float
        The luminosity distance
    a_1: float
        The dimensionless spin magnitude of the primary black hole
    tilt_1: float
        The tilt angle of the primary black hole spin
    phi_12: float
        The difference in azimuthal angle between the two spins
    a_2: float
        The dimensionless spin magnitude of the secondary black hole
    tilt_2: float
        The tilt angle of the secondary black hole spin
    phi_jl: float
        The azimuthal angle of the total angular momentum
    theta_jn: float
        The angle between the total angular momentum and the line of sight
    phase: float
        The phase of the gravitational wave
    kwargs: dict
        Additional keyword arguments

    Returns
    -------
    polarizations: dict
        The polarizations
    """
    import importlib
    try:
        _split = interpolant.split(".")
        _module = ".".join(_split[:-1])
        _function = _split[-1]
        module = importlib.import_module(_module)
        method = getattr(module, _function)
    except Exception as e:
        raise ValueError(f"Unable to import interpolant function because: {e}")

    _matches = np.array(
        [
            method(
                wvf, mass_1, mass_2, a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl,
                theta_jn, phase, 
            ) for wvf in waveform_approximant_list
        ]
    )
    # protect against negative matches
    _matches[_matches < 0.] = 0.
    use_best = kwargs.pop("use_best_match", False)
    mapping = kwargs.pop("match_to_weight", None)
    if isinstance(use_best, str):
        use_best = ast.literal_eval(use_best)
    _weights = _weights_from_matches(_matches, use_best=use_best, mapping=mapping)
    return _multi_model_binary_black_hole(
        _weights, waveform_approximant_list, frequency_array, mass_1, mass_2,
        luminosity_distance, a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl,
        theta_jn, phase, **kwargs
    )


def _multi_model_binary_black_hole(
    weights, waveform_approximant_list, frequency_array, mass_1, mass_2,
    luminosity_distance, a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn,
    phase, **kwargs
):
    """Source model for a binary black hole with multiple models. The model
    is chosen based on the weights provided.

    Parameters
    ----------
    weights: list
        The interpolant you wish to use to estimate the match
    waveform_approximant_list: list
        list of waveform approximants you wish to include. Must equal the length
        of weights
    frequency_array: np.ndarray
        The frequency array
    mass_1: float
        The mass of the primary black hole
    mass_2: float
        The mass of the secondary black hole
    luminosity_distance: float
        The luminosity distance
    a_1: float
        The dimensionless spin magnitude of the primary black hole
    tilt_1: float
        The tilt angle of the primary black hole spin
    phi_12: float
        The difference in azimuthal angle between the two spins
    a_2: float
        The dimensionless spin magnitude of the secondary black hole
    tilt_2: float
        The tilt angle of the secondary black hole spin
    phi_jl: float
        The azimuthal angle of the total angular momentum
    theta_jn: float
        The angle between the total angular momentum and the line of sight
    phase: float
        The phase of the gravitational wave
    kwargs: dict
        Additional keyword arguments

    Returns
    -------
    polarizations: dict
        The polarizations
    """
    if not np.any(weights):
        raise ValueError(
            "All weights are non-numeric. Please provide a numeric weight for "
            "each approximant in the list"
        )
    elif len(weights) != len(waveform_approximant_list):
        raise ValueError(
            "Please provide a weight for each approximant in the list"
        )
    else:
        waveform_approximant = np.random.choice(
            waveform_approximant_list, p=weights
        )
    kwargs["waveform_approximant"] = waveform_approximant
    # only use gwsignal for reviewed waveforms. This should be changed when
    # bilby updates their review statement
    if waveform_approximant in ["SEOBNRv5HM", "SEOBNRv5PHM"]:
        return source.gwsignal_binary_black_hole(
            frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
            phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
        )
    if waveform_approximant == "IMRPhenomXPHMST":
        kwargs.update(
            {
                "waveform_approximant": "IMRPhenomXPHM",
                "PhenomXPrecVersion": 320,
                "PhenomXPFinalSpinMod": 2,
                "PhenomXHMReleaseVersion": 122022
            }
        )
    return source.lal_binary_black_hole(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
    )


def _weights_from_matches(matches, use_best=False, mapping=None):
    """Calculate a weight based on the match to numerical relativity. If
    mapping is None, we use the recommendation from
    https://www.nature.com/articles/s41550-025-02579-7,
    i.e. weight = 1 / (1 - match)**4

    Parameters
    ----------
    matches: np.ndarray
        array of matches to calculate the weight for
    use_best: bool, optional
        if True, return a weight of 1 for the highest match and 0 for all other
        models. Default False
    mapping: str, optional
        a string that can be evaluated to map an array of matches to a series of
        weights. The string must contain the variable 'matches' and should
        contain the right hand side of the equation 'weights = f(matches)'.
        For example, you could provide '1 / ((1 - matches)**4)'. The weights
        will be rescaled to be between 0 and 1.
    """
    if use_best:
        weights = np.zeros(len(matches))
        weights[np.argmax(matches)] = 1.
        return weights
    if mapping is None:
        weights = 1 / ((1 - matches)**4)
    else:
        if "matches" not in mapping:
            raise ValueError(
                f"{mapping} must be a string that can be evaluated and contain "
                f"a variable called: 'match'"
            )
        if "=" in mapping:
            raise ValueError(
                f"{mapping} must not contain '='. It should contain the "
                f"right hand side of the equation 'weights = f(matches)'. For "
                f"example, you could provide '1 / ((1 - matches)**4)'"
            )
        try:
            weights = eval(mapping)
        except Exception as e:
            raise ValueError(
                f"Unable to generate weights from matches for the string "
                f"{mapping}."
            )
    weights /= np.sum(weights)
    weights[np.isnan(weights)] = 0.
    return weights
