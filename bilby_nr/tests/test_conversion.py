import bilby
from bilby.gw import conversion
from gwpy.timeseries import TimeSeries
from bilby_nr.conversion import (
    determine_waveform_approximant_from_likelihood
)
import pytest


def _setup(model="IMRPhenomTPHM"):
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
            "waveform_approximant": model,
            "reference_frequency": 50,
            "minimum_frequency": 20,
        },
    )
    likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
        ifo_list, waveform_generator
    )
    return likelihood, ifo_list, waveform_generator


def test_generate_all_bbh_parameters():
    samples = {
        "mass_1": [36, 35, 37],
        "mass_2": [32, 31, 33],
        "a_1": [0.8, 0.5, 0.7],
        "a_2": [0.5, 0.2, 0.8],
        "tilt_1": [2.4, 1.2, 4.5],
        "tilt_2": [0.8, 0.6, 2.1],
        "phi_12": [5.8, 4.0, 2.1],
        "phi_jl": [0.25, 0.4, 0.1],
        "theta_jn": [2.5, 2.4, 2.6],
        "ra": [2.2, 2.1, 2.3],
        "dec": [-1.22, -1.2, -1.25],
        "geocent_time": [1126259462.408404, 1126259462.408404, 1126259462.408404],
        "phase": [4., 5., 2.],
        "psi": [0.7, 0.5, 0.4],
        "luminosity_distance": [500, 505, 502],
    }
    waveform_approximants = ["IMRPhenomXPHM", "IMRPhenomTPHM", "IMRPhenomPv2"]
    logls = []
    for num in range(len(samples["mass_1"])):
        likelihood, _, _ = _setup(model=waveform_approximants[num])
        _sample = {key: item[num] for key, item in samples.items()}
        logls.append(likelihood.log_likelihood_ratio(_sample))
    samples["log_likelihood"] = logls
    likelihood.waveform_generator.waveform_arguments.update(
        {"waveform_approximant_list": waveform_approximants}
    )
    models = []
    for num in range(len(samples["mass_1"])):
        _sample = {key: item[num] for key, item in samples.items()}
        sample = conversion.generate_all_bbh_parameters(_sample, likelihood)
        models.append(sample["waveform_approximant"])
    assert models == waveform_approximants


def test_determine_waveform_approximant_from_likelihood():
    import bilby
    from gwpy.timeseries import TimeSeries
    from pandas import DataFrame

    likelihood, ifo_list, waveform_generator = _setup()
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

    # try passing a pandas dataframe instead
    _sample = DataFrame(sample, index=[0])
    model = determine_waveform_approximant_from_likelihood(
        _sample, ["IMRPhenomPv2", "IMRPhenomTPHM", "IMRPhenomXPHM"],
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
