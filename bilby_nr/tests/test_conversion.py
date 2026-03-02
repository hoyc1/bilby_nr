from bilby_nr.conversion import (
    determine_waveform_approximant_from_likelihood
)
import pytest


def test_determine_waveform_approximant_from_likelihood():
    import bilby
    from gwpy.timeseries import TimeSeries

    ifo_list = bilby.gw.detector.InterferometerList([])
    for det in ["H1", "L1"]:
        ifo = bilby.gw.detector.get_empty_interferometer(det)
        data = TimeSeries.fetch_open_data(
            det, 1126259462.4 - 2., 1126259462.4 + 2.
        )
        ifo.strain_data.set_from_gwpy_timeseries(data)
        ifo_list.append(ifo)

    waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        waveform_arguments={
            "waveform_approximant": "IMRPhenomTPHM",
            "reference_frequency": 50,
            "minimum_frequency": 20,
        },
    )
    likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
        ifo_list, waveform_generator
    )
    sample = {
        "mass_1": 36, "mass_2": 32, "a_1": 0.8, "a_2": 0.7, "tilt_1": 2.4,
        "tilt_2": 0.8, "phi_12": 5.8, "phi_jl": 0.25, "theta_jn": 2.5,
        "ra": 2.2, "dec": -1.22, "geocent_time": 1126259462.408404,
        "phase": 4., "psi": 0.7, "luminosity_distance": 500
    }
    logl = likelihood.log_likelihood_ratio(sample)
    with pytest.raises(ValueError):
        model = determine_waveform_approximant_from_likelihood(
            sample, ["IMRPhenomPv2", "IMRPhenomTPHM", "IMRPhenomXPHM"],
            likelihood
        )
    sample["log_likelihood"] = logl
    model = determine_waveform_approximant_from_likelihood(
        sample, ["IMRPhenomPv2", "IMRPhenomTPHM", "IMRPhenomXPHM"],
        likelihood
    )
    assert model == "IMRPhenomTPHM"

    # try now with marginalization
    priors = bilby.gw.prior.BBHPriorDict()
    for marg in ["distance", "time", "phase"]:
        kwargs = {}
        _sample = sample.copy()
        if marg == "distance":
            kwargs["distance_marginalization"] = True
        elif marg == "time":
            kwargs["time_marginalization"] = True
            _sample["geocent_time"] = likelihood.interferometers.start_time
            _sample["time_jitter"] = 2.
        else:
            kwargs["phase_marginalization"] = True
            _sample["phase"] = 0.
        likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
            ifo_list, waveform_generator, priors=priors, **kwargs
        )
        logl = likelihood.log_likelihood_ratio(_sample)
        _sample["log_likelihood"] = logl
        model = determine_waveform_approximant_from_likelihood(
            _sample, ["IMRPhenomPv2", "IMRPhenomTPHM", "IMRPhenomXPHM"],
            likelihood
        )
        assert model == "IMRPhenomTPHM"
