# Licensed under an MIT style license -- see LICENSE.md

import numpy as np
from bilby.gw.conversion import (
    bilby_to_lalsimulation_spins, component_masses_to_symmetric_mass_ratio
)

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
    spins = bilby_to_lalsimulation_spins(
        theta_jn, phi_jl, tilt_1, tilt_2, phi_12, a_1, a_2, mass_1, mass_2,
        reference_frequency, phase
    )
    S1_perp = mass_1**2 * np.array([spins[1], spins[2]])
    S2_perp = mass_2**2 * np.array([spins[4], spins[5]])
    S_perp = np.sum([S1_perp, S2_perp], axis=0)
    S_perp_mag = np.linalg.norm(S_perp)
    chi_perp = S_perp_mag / total_mass**2
    return chi_perp
