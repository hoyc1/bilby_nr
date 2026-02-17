# Licensed under an MIT style license -- see LICENSE.md

import numpy as np

__author__ = ["Charlie Hoy <charlie.hoy@port.ac.uk>"]


def match_from_interpolant(*args, interp="pade_pade"):
    """Return an estimate for the match based on an interpolant fit

    Parameters
    ----------
    *args: list
        The arguments to pass to the interpolant
    interp: str
        The name of the interpolant to use. Default is "pade_pade"

    Returns
    -------
    match: float
        The estimated match
    """
    if interp not in interpolant_map.keys():
        raise ValueError(
            f"Unable to evaluate the mismatch interpolant using {interp}. "
            f"Only allowed interpolants are: {','.join(interpolant_map.keys())}"
        )
    return interpolant_map[interp](*args)


def match_from_pade_pade_interpolant(
    waveform_approximant, mass_1, mass_2, a_1, tilt_1, phi_12, a_2, tilt_2,
    phi_jl, theta_jn, phase,
):
    """Return an estimate for the match based on a Pade-Pade fit, see
    https://www.nature.com/articles/s41550-025-02579-7

    Parameters
    ----------
    waveform_approximant: str
        The waveform approximant to use
    mass_1: float
        The mass of the primary black hole
    mass_2: float
        The mass of the secondary black hole
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

    Returns
    -------
    match: float
        The estimated match
    """
    from .interp.pade_pade import match_interpolant
    if isinstance(mass_1, (np.ndarray, list)):
        matches = np.zeros(len(mass_1))
        for ii in range(len(mass_1)):
            matches[ii] = match_interpolant(
                waveform_approximant, mass_1[ii], mass_2[ii], a_1[ii],
                tilt_1[ii], phi_12[ii], a_2[ii], tilt_2[ii], phi_jl[ii],
                theta_jn[ii], phase[ii]
            )
        return matches
    return match_interpolant(
        waveform_approximant, mass_1, mass_2, a_1, tilt_1, phi_12, a_2,
        tilt_2, phi_jl, theta_jn, phase
    )


interpolant_map = {
    "pade_pade": match_from_pade_pade_interpolant,
}
