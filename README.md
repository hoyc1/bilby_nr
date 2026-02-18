# bilby NR

## Development status

[![Coverage report](https://hoyc1.github.io/bilby_nr/coverage-badge.svg)](https://hoyc1.github.io/bilby_nr/coverage.xml) [![Pipeline Status](https://github.com/hoyc1/bilby_nr/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/hoyc1/bilby_nr/actions/workflows/test.yml)

This Python package incorporates model accuracy into gravitational wave Bayesian analyses via the bilby Python package.

## Installation

`bilby_nr` is currently available via PyPI and can be installed with:

```bash
pip install bilby_nr
```

For full installation instructions, see [our documentation](https://hoyc1.github.io/bilby_nr/installation.html).

## Usage in bilby_pipe

The functionality in `bilby_nr` can be used with `bilby_pipe` as you would with any other frequency domain source model. It simply requires the following options to be specified in your configuration file:

```ini
analysis_executable_parser=bilby_nr.bilby_pipe.create_parser
waveform-approximant=IMRPhenomXPHMST,IMRPhenomTPHM,SEOBNRv5PHM
frequency-domain-source-model = bilby_nr.source.multi_model_binary_black_hole
waveform-arguments-dict={'match_interpolant': 'bilby_nr.match.match_from_pade_pade_interpolant'}
```

## Citing

If you find `bilby_nr` useful in your work please cite the following papers:

```bibtex
@article{Hoy:2024vpc,
    author = "Hoy, Charlie and Akcay, Sarp and Mac Uilliam, Jake and Thompson, Jonathan E.",
    title = "{Incorporation of model accuracy in gravitational wave Bayesian inference}",
    eprint = "2409.19404",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    reportNumber = "LIGO-P2400393",
    doi = "10.1038/s41550-025-02579-7",
    journal = "Nature Astron.",
    volume = "9",
    number = "8",
    pages = "1256--1267",
    year = "2025"
}

@article{Hoy:2022tst,
    author = "Hoy, Charlie",
    title = "{Accelerating multimodel Bayesian inference, model selection, and systematic studies for gravitational wave astronomy}",
    eprint = "2208.00106",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    reportNumber = "LIGO-P2200228",
    doi = "10.1103/PhysRevD.106.083003",
    journal = "Phys. Rev. D",
    volume = "106",
    number = "8",
    pages = "083003",
    year = "2022"
}
```
