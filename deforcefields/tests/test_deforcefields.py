"""
Test loading DE-Force fields via the plugin interface through the toolkit.
"""

import openmm
import pytest
from openff.interchange.drivers.openmm import get_openmm_energies
from openff.toolkit import ForceField, Molecule
from openmm import unit


@pytest.mark.parametrize(
    "forcefield",
    [
        pytest.param("de-force-1.0.2.offxml", id="No constraints"),
        pytest.param("de-force_unconstrained-1.0.2.offxml", id="Constraints"),
    ],
)
def test_load_de_ff(forcefield):
    """
    Load the DE FF and create an OpenMM system.
    """

    ff = ForceField(forcefield, load_plugins=True)
    ethanol = Molecule.from_smiles("CCO")

    system = ff.create_interchange(topology=ethanol.to_topology()).to_openmm(
        combine_nonbonded_forces=False,
    )

    forces = {force.__class__.__name__: force for force in system.getForces()}

    # make sure scaled 1-4 interactions are present
    assert "CustomBondForce" in forces
    assert "CustomNonbondedForce" in forces

    # make sure the sigma and epsilon are dummy parameters
    nonbonded: openmm.NonbondedForce = forces["NonbondedForce"]
    for i in range(nonbonded.getNumParticles()):
        _, sigma, epsilon = nonbonded.getParticleParameters(i)
        assert sigma == 0 * unit.nanometer
        assert epsilon == 0 * unit.kilocalorie_per_mole

    # check the custom force parameters are non-zero
    custom_force: openmm.CustomNonbondedForce = forces["CustomNonbondedForce"]
    for i in range(custom_force.getNumParticles()):
        epsilon, _ = custom_force.getParticleParameters(i)
        # no units with our custom force
        assert epsilon > 0

    assert custom_force.getNonbondedMethod() == openmm.NonbondedForce.NoCutoff


def test_openmm_energies_not_crazy():
    # Could pick any ligand that's supported by both FFs
    molecule = Molecule.from_smiles("NC(=O)c1cccc2c1OCCO2")
    molecule.generate_conformers(n_conformers=1)
    topology = molecule.to_topology()

    sage_energies = get_openmm_energies(
        ForceField("openff-2.1.0.offxml").create_interchange(topology),
        combine_nonbonded_forces=False,
    )

    de_energies = get_openmm_energies(
        ForceField(
            "de-force-1.0.2.offxml",
            load_plugins=True,
        ).create_interchange(topology),
        combine_nonbonded_forces=False,
    )

    # Energies should differ a little bit, but be similar
    assert 0.1 < sage_energies["vdW"] / de_energies["vdW"] < 10

    # Electrostatics methods are the same, so energies should be as well
    assert sage_energies["Electrostatics"].m == pytest.approx(
        de_energies["Electrostatics"].m
    )


def test_fails_unsupported_chemistry():
    # No double exponential parameters for H bonded to C#,
    # so this should error instead of assigning blank/bogus parameters
    with pytest.raises(KeyError, match="atom indices .*2"):
        ForceField(
            "de-force-1.0.2.offxml",
            load_plugins=True,
        ).create_interchange(
            Molecule.from_mapped_smiles("[H:3][C:1]#[C:2][H:4]").to_topology()
        ).to_openmm(combine_nonbonded_forces=False)
