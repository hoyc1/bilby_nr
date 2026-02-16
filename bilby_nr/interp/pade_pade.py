# Licensed under an MIT style license -- see LICENSE.md

import numpy as np
import os
from ..conversion import (
    chi_par_chi_perp_from_mass_spin, component_masses_to_symmetric_mass_ratio
)

__author__ = ["Charlie Hoy <charlie.hoy@port.ac.uk>"]


def Cijkl(i, j, k, l, fit_coeffs):
    """Return coefficients of the Pade Pade fit as defined in Eq.X of
    https://www.nature.com/articles/s41550-025-02579-7 given data structure
    """
    ind, = np.where(
        (fit_coeffs["i"] == i) & (fit_coeffs["j"] == j) &
        (fit_coeffs["k"] == k) & (fit_coeffs["l"] == l)
    )
    if len(ind) > 1:
        raise ValueError("Duplicated entry")
    return fit_coeffs["Cijkl"][ind[0]]


def ArcTan(x, y):
    """Convert Mathematica definition of ArcTan to numpy definition"""
    return np.arctan2(y, x)


def match_interpolant(
    waveform_approximant, mass_1, mass_2, a_1, tilt_1, phi_12, a_2, tilt_2,
    phi_jl, theta_jn, phase
):
    """Evaluate match interpolant as defined in Eq.X of
    https://www.nature.com/articles/s41550-025-02579-7

    Parameters
    ----------
    waveform_approximant: str
        Name of waveform approximant you wish to evaluate the interpolant for.
        Currently allowed waveform approximants include
        ["IMRPhenomXPHMST", "IMRPhenomTPHM", "SEOBNRv5PHM"]
    mass_1: float
        Detector-frame primary mass of the binary black hole
    mass_2: float
        Detector-frame secondary mass of the binary black hole
    a_1: float
        Magnitude of the spin vector associated with the primary component in
        the binary black hole
    tilt_1: float
        Polar angle of the spin vector associated with the primary component in
        the binary black hole
    phi_12: float
        Azimuthal angle between the primary spin and secondary spin vector
    a_2: float
        Magnitude of the spin vector associated with the secondary component
        in the binary black hole
    tilt_2: float
        Polar angle of the spin vector associated with the secondary component
        in the binary black hole
    phi_jl: float
        Azimuthal angle between the total angular momentum (J) and the orbital
        angular momentum (L)
    theta_jn: float
        Inclination angle of the binary: the angle between the total angular
        momentum (J) and the line
        of sight (N)
    phase: float
        The phase of the binary black hole

    Returns
    -------
    match: float
        An approximate match for a given model
    """
    return 1 - mismatch_interpolant(
        waveform_approximant, mass_1, mass_2, a_1, tilt_1, phi_12, a_2, tilt_2,
        phi_jl, theta_jn, phase
    )


def mismatch_interpolant(
    waveform_approximant, mass_1, mass_2, a_1, tilt_1, phi_12, a_2, tilt_2,
    phi_jl, theta_jn, phase
):
    """Evaluate mismatch interpolant as defined in Eq.X of
    https://www.nature.com/articles/s41550-025-02579-7

    Parameters
    ----------
    waveform_approximant: str
        Name of waveform approximant you wish to evaluate the interpolant for.
        Currently allowed waveform approximants include
        ["IMRPhenomXPHMST", "IMRPhenomTPHM", "SEOBNRv5PHM"]
    mass_1: float
        Detector-frame primary mass of the binary black hole
    mass_2: float
        Detector-frame secondary mass of the binary black hole
    a_1: float
        Magnitude of the spin vector associated with the primary component in
        the binary black hole
    tilt_1: float
        Polar angle of the spin vector associated with the primary component in
        the binary black hole
    phi_12: float
        Azimuthal angle between the primary spin and secondary spin vector
    a_2: float
        Magnitude of the spin vector associated with the secondary component
        in the binary black hole
    tilt_2: float
        Polar angle of the spin vector associated with the secondary component
        in the binary black hole
    phi_jl: float
        Azimuthal angle between the total angular momentum and the orbital
        angular momentum
    theta_jn: float
        Inclination angle of the binary: the angle between the total angular
        momentum and the line of sight
    phase: float
        The phase of the binary black hole

    Returns
    -------
    mismatch: float
        An approximate mismatch for a given model
    """
    allowed_models = ["IMRPhenomXPHMST", "IMRPhenomTPHM", "SEOBNRv5PHM"]
    if waveform_approximant not in allowed_models:
        raise ValueError(
            f"Unable to evaluate interpolant for waveform model "
            f"{waveform_approximant}. Please provide either 'IMRPhenomXPHMST', "
            f"'IMRPhenomTPHM' or 'SEOBNRv5PHM'."
        )
    if mass_1 < mass_2:
        raise ValueError(
            "Secondary mass must be smaller than the primary mass of the binary"
        )

    Mtot = mass_1 + mass_2
    eta = component_masses_to_symmetric_mass_ratio(mass_1, mass_2)
    # interpolant uses f_ref = 10.70629431812844 Hz. We use the same
    # reference frequency here
    chi_par, chi_perp = chi_par_chi_perp_from_mass_spin(
        mass_1, mass_2, a_1, tilt_1, a_2, tilt_2, phi_12, phi_jl, theta_jn,
        phase, 10.70629431812844
    )
    # Read in file containing fitting coefficients
    if waveform_approximant == "IMRPhenomXPHMST":
        identifier = "XPHMST"
    elif waveform_approximant == "IMRPhenomTPHM":
        identifier = "TPHM"
    else:
        identifier = "SEOB"

    filename = os.path.join(
        os.path.dirname(__file__),
        f"NatureAstronomy.XXX.YYY.2025.{identifier}.txt"
    )
        
    with open(filename, "r") as f:
        lines = f.readlines()
        lines = [l.strip() for l in lines]
        original_variables = lines[1].split("\t")
        changed_variables = lines[2].split("\t")

    # Evaluate the fitting co-efficients and return the mismatch
    fit_coeffs = np.genfromtxt(filename, skip_header=3, names=True)
    Sqrt = np.sqrt
    x, y, z, v = original_variables[:4]
    X, Y, Z, V = changed_variables[2:6]

    for var in [x, y, z, v]:
        X = X.replace(var.split("=")[0], f"({var.split('=')[1]})")
        Y = Y.replace(var.split("=")[0], f"({var.split('=')[1]})")
        Z = Z.replace(var.split("=")[0], f"({var.split('=')[1]})")
        V = V.replace(var.split("=")[0], f"({var.split('=')[1]})")

    X = eval(X.split("=")[1])
    Y = eval(Y.split("=")[1])
    Z = eval(Z.split("=")[1])
    V = eval(V.split("=")[1])

    log10_mismatch = np.sum([
        (X)**(i) * (Y)**(j) * np.sum([
            Cijkl(i, j, k, l, fit_coeffs) * (Z)**(k) * (V)**(l)
            for l in np.arange(0, np.max(fit_coeffs["l"]) - 1)
            for k in np.arange(0, np.max(fit_coeffs["k"]) + 1)
        ]) / np.sum([
            np.abs(Cijkl(i, j, k, l + 2, fit_coeffs)) * (Z)**(k) * (V)**(l)
            for l in np.arange(0, np.max(fit_coeffs["l"]) - 1)
            for k in np.arange(0, np.max(fit_coeffs["k"]) + 1)
        ])
        for j in np.arange(0, np.max(fit_coeffs["j"]) + 1)
        for i in np.arange(0, np.max(fit_coeffs["i"]) + 1)
    ])
    mismatch = 10**log10_mismatch
    return mismatch
