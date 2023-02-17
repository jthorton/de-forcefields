# DE-Force Fields
Transferable Double Exponential non-bonded potential for condensed phase simulations of small molecules.


[![CI](https://github.com/jthorton/de-forcefields/actions/workflows/CI.yaml/badge.svg)](https://github.com/jthorton/de-forcefields/actions/workflows/CI.yaml)

This repository contains general transferable force fields that use the Double Exponential functional form first proposed by [Brooks](http://aip.scitation.org/doi/10.1063/1.5107505). 
These force fields use the SMIRKS Native Open Force Field (SMIRNOFF) format. 
By convention these files use the `.offxml` file extension. 
The SMIRNOFF format has a [specification](https://openforcefield.github.io/standards/standards/smirnoff/) and is discussed in a [JCTC publication](https://doi.org/10.1021/acs.jctc.8b00640) and associated [pre-print](https://www.biorxiv.org/content/10.1101/286542v3).

The [OpenFF Toolkit](https://github.com/openforcefield/openff-toolkit) (version >0.10.6, <0.11.0) provides a reference implementation of the SMIRNOFF format. 
In particular, the `ForceField` class is used to load SMIRNOFF-format force fields and the `create_openmm_system` method enables the parametrization of small molecules into OpenMM objects.
[Smirnoff-plugins](https://github.com/openforcefield/smirnoff-plugins) provides a framework to extend the SMIRNOFF specification with custom force field functional forms such as the Double Exponential form used here, using a plugin system. 
See `smirnoff-plugins` for a list of the currently supported potentials.

Detailed usage examples can be found in the OpenFF Toolkit repository.

Each force field is currently available in two forms -- both with and without bond constraints to hydrogen. The default version of each force field (i.e. de-1.0.0.offxml) is suitable for typical molecular dynamics simulations with constrained bonds to hydrogen. 
The "unconstrained" version of each force field (i.e. de_unconstrained-1.0.0.offxml) should be used when single-point energies are a major concern (e.g. geometry optimizations) and when comparing the force field to QM data.

| Filename                              | DOI  | FF line    | Release Date | Major format changes? |
|---------------------------------------|------|------------|--------------|-----------------------|
| `de-force-1.0.0.offxml`               | TODO | DE-Force-1 | Feb 17, 2023 | No                    |
| `de-force_unconstrained-1.0.0.offxml` | TODO | DE-Force-1 | Feb 17, 2023 | No                    |


## Installation

```shell
conda install -c conda-forge de-forcefields
```


## Use

Installing this package exposes an entry point that makes the `deforcefields/offxml/` directory easily accessible by other packages in the same Python installation. 
If the [OpenFF Toolkit](https://github.com/openforcefield/openff-toolkit) is installed, it will automatically detect and use this entry point:

```python3
>>> from openff.toolkit.typing.engines.smirnoff import ForceField
>>> ff = ForceField('de-force-1.0.0.offxml', load_plugins=True)
```

Otherwise, the entry point can be accessed by querying the `openforcefield.smirnoff_forcefield_directory` entry point group.

```python3
>>> from pkg_resources import iter_entry_points
>>> for entry_point in iter_entry_points(group='openforcefield.smirnoff_forcefield_directory'):
...     print(entry_point.load()())
```


# History

Force fields in the `DE-Force-1` lines are descended from [OpenFF-2.0.0 Sage](https://doi.org/10.5281/zenodo.5214478).

## General versioning guidelines

Force fields moving forward will be called `name-X.Y.Z`

* `X` denotes some major change in functional form or fitting strategy.
* `Y` is the parameterization epoch / generation, or a minor change that can affect energy.
* `Z` is a bugfix version -- e.g. something we've caught and corrected.


## Versions
- `DE-Force-1` : Proof of concept general transferable Double Exponential force field fit using [Sage training data](https://doi.org/10.26434/chemrxiv-2022-n2z1c-v2). 
