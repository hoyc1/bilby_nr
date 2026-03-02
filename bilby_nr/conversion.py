# Licensed under an MIT style license -- see LICENSE.md

import copy
import numpy as np
from bilby.gw import conversion
from bilby.gw.conversion import _generate_all_cbc_parameters as _base_generate

__author__ = ["Charlie Hoy <charlie.hoy@port.ac.uk>"]


def chi_par_chi_perp_from_mass_spin(
    mass_1, mass_2, a_1, tilt_1, a_2, tilt_2, phi_12, phi_jl, theta_jn, phase,
    reference_frequency
):
    """Calculate and return the parallel and perpendicular spin, as defined in
    Eqs. 11 and 12 of https://www.nature.com/articles/s41550-025-02579-7

    Parameters
    ----------
    mass_1: float
        The mass of the primary black hole
    mass_2: float
        The mass of the secondary black hole
    a_1: float
        The dimensionless spin magnitude of the primary black hole
    tilt_1: float
        The tilt angle of the primary black hole spin
    a_2: float
        The dimensionless spin magnitude of the secondary black hole
    tilt_2: float
        The tilt angle of the secondary black hole spin
    phi_12: float
        The difference in azimuthal angle between the two spins
    phi_jl: float
        The azimuthal angle of the total angular momentum
    theta_jn: float
        The angle between the total angular momentum and the line of sight
    phase: float
        The phase of the gravitational wave
    reference_frequency: float
        The reference frequency

    Returns
    -------
    chi_par: float
        The parallel spin
    chi_perp: float
        The perpendicular spin
    """
    _chi_par = chi_par(mass_1, mass_2, a_1, tilt_1, a_2, tilt_2)
    _chi_perp = chi_perp(
        mass_1, mass_2, a_1, tilt_1, a_2, tilt_2, phi_12, phi_jl, theta_jn,
        phase, reference_frequency
    )
    return _chi_par, _chi_perp


def chi_par(mass_1, mass_2, a_1, tilt_1, a_2, tilt_2):
    """Calculate and return the parallel spin, as defined in Eq.12 of
    https://www.nature.com/articles/s41550-025-02579-7

    Parameters
    ----------
    mass_1: float
        The mass of the primary black hole
    mass_2: float
        The mass of the secondary black hole
    a_1: float
        The dimensionless spin magnitude of the primary black hole
    tilt_1: float
        The tilt angle of the primary black hole spin
    a_2: float
        The dimensionless spin magnitude of the secondary black hole
    tilt_2: float
        The tilt angle of the secondary black hole spin

    Returns
    -------
    chi_par: float
        The parallel spin
    """
    mass_ratio = mass_2 / mass_1
    numerator = (
        a_1 * np.cos(tilt_1) + mass_ratio**2 * a_2 * np.cos(tilt_2)
    )
    denominator = (1 + mass_ratio)**2
    return numerator / denominator


def chi_perp(
    mass_1, mass_2, a_1, tilt_1, a_2, tilt_2, phi_12, phi_jl, theta_jn, phase,
    reference_frequency
):
    """Calculate and return the perpependicular spin, as defined in Eq. 11 of
    https://www.nature.com/articles/s41550-025-02579-7

    Parameters
    ----------
    mass_1: float
        The mass of the primary black hole
    mass_2: float
        The mass of the secondary black hole
    a_1: float
        The dimensionless spin magnitude of the primary black hole
    tilt_1: float
        The tilt angle of the primary black hole spin
    a_2: float
        The dimensionless spin magnitude of the secondary black hole
    tilt_2: float
        The tilt angle of the secondary black hole spin
    phi_12: float
        The difference in azimuthal angle between the two spins
    phi_jl: float
        The azimuthal angle of the total angular momentum
    theta_jn: float
        The angle between the total angular momentum and the line of sight
    phase: float
        The phase of the gravitational wave
    reference_frequency: float
        The reference frequency

    Returns
    -------
    chi_perp: float
        The perpendicular spin
    """
    total_mass = mass_1 + mass_2
    spins = conversion.bilby_to_lalsimulation_spins(
        theta_jn, phi_jl, tilt_1, tilt_2, phi_12, a_1, a_2, mass_1, mass_2,
        reference_frequency, phase
    )
    S1_perp = mass_1**2 * np.array([spins[1], spins[2]])
    S2_perp = mass_2**2 * np.array([spins[4], spins[5]])
    S_perp = np.sum([S1_perp, S2_perp], axis=0)
    S_perp_mag = np.linalg.norm(S_perp)
    chi_perp = S_perp_mag / total_mass**2
    return chi_perp


def _generate_all_cbc_parameters(
    sample, defaults, base_conversion, likelihood=None, priors=None, npool=1
):
    """Extension of bilby.gw.conversion.generate_all_bbh_parameters to allow
    for multiple models

    Parameters
    ==========
    sample: dict, pandas.DataFrame
        Samples to fill in with extra parameters, this may be either an
        injection or posterior samples.
    likelihood: bilby.gw.likelihood.GravitationalWaveTransient, optional
        GravitationalWaveTransient used for sampling, used for waveform and
        likelihood.interferometers.
    priors: dict, optional
        Dictionary of prior objects, used to fill in non-sampled parameters.
    npool: int, optional
        Number of processes to use for the conversion. Default 1
    """
    _chosen_model = None
    if "log_likelihood" in sample.keys() and likelihood is not None:
        waveform_approximant_list = \
            likelihood.waveform_generator.waveform_arguments.get(
                "waveform_approximant_list", None
            )
        if waveform_approximant_list is not None:
            _chosen_model = determine_waveform_approximant_from_likelihood(
                sample, waveform_approximant_list, likelihood
            )

    if _chosen_model is None:
        # assume a default
        _chosen_model = "IMRPhenomTPHM"
    _likelihood = copy.deepcopy(likelihood)
    if _chosen_model is not None:
        if likelihood is not None:
            _likelihood.waveform_generator.waveform_arguments.update(
                {
                    "waveform_approximant_list": [_chosen_model],
                    "waveform_approximant": _chosen_model
                }
            )
        sample["waveform_approximant"] = _chosen_model
        defaults["waveform_approximant"] = _chosen_model
    return _base_generate(
        sample, defaults=defaults,
        base_conversion=base_conversion,
        likelihood=_likelihood, priors=priors, npool=npool
    )


def determine_waveform_approximant_from_likelihood(
    sample, waveform_approximant_list, likelihood
):
    """Determine the waveform approximant that evaluated the likelihood during
    the inference.

    Parameters
    ----------
    sample: dict, pandas.DataFrame
        the sample you wish to evaluate the likelihood for
    waveform_approximant_list: list
        list of waveform approximants you wish to consider
    likelihood: bilby.gw.likelihood.GravitationalWaveTransient
        likelihood object that was used during the sampling
    """
    _chosen_model = None
    if "log_likelihood" not in sample.keys():
        raise ValueError(
            "Unable to determine which model ws used to evaluate the "
            "likelihood as there is no 'log_likelihood' entry in the provided "
            "samples."
        )
    original = sample["log_likelihood"]
    if not isinstance(original, (float, int, np.number)):
        if isinstance(sample, dict):
            original = original[0]
        else:
            original = original.values[0]

    for model in waveform_approximant_list:
        _likelihood = copy.deepcopy(likelihood)
        logl = _likelihood_for_given_model(sample, model, _likelihood)
        if np.isclose(logl, original):
            _chosen_model = model
            break
    if _chosen_model is None:
        raise ValueError(
           "Unable to find a model that returns the same log likelihood as "
           "that stored in the sample."
        )
    return _chosen_model


def _likelihood_for_given_model(sample, waveform_approximant, likelihood):
    """Evaluate the likelihood for a given waveform approximant

    Parameters
    ----------
    sample: dict, pandas.DataFrame
        the sample you wish to evaluate the likelihood for
    waveform_approximant: str
        the waveform approximant you wish to evaluate the likelihood for
    likelihood: bilby.gw.likelihood.GravitationalWaveTransient
        likelihood object that was used during the sampling
    """
    _lkl = copy.deepcopy(likelihood)
    # remove cache of last sample by adding a small shift to all parameters
    _sample = sample.copy()
    for key in _sample.keys():
        _sample[key] += np.random.random() * 1e-5
    _lkl.log_likelihood_ratio(_sample)
    # now check the sample provided
    _lkl.waveform_generator.waveform_arguments["waveform_approximant_list"] = [
        waveform_approximant
    ]
    _lkl.waveform_generator.waveform_arguments["waveform_approximant"] = \
        waveform_approximant
    if isinstance(sample, dict):
        _sample = {
            key: item for key, item in sample.items() if key not in
            ["log_likelihood", "log_prior"]
        }
        if not isinstance(list(_sample.values())[0], (float, int, np.number)):
            _sample = {key: item[0] for key, item in _sample.items()}
    else:
        _sample = {
            key: item.values[0] for key, item in sample.items() if key not in
            ["log_likelihood", "log_prior"]
        }
            
    if _lkl.distance_marginalization:
        _sample["luminosity_distance"] = _lkl._ref_dist
    if _lkl.time_marginalization:
        _sample["geocent_time"] = _lkl.interferometers.start_time
        if "time_jitter" not in _sample.keys():
            raise ValueError(
                "Please provide a jitter time as time marginalization was "
                "used during the sampling"
            )
    if _lkl.phase_marginalization:
        _sample["phase"] = 0.

    logl = _lkl.log_likelihood_ratio(parameters=_sample)
    return logl


conversion._generate_all_cbc_parameters = _generate_all_cbc_parameters
