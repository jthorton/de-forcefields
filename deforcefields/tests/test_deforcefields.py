"""
Test loading DE-Force fields via the plugin interface through the toolkit.
"""
import openmm
import pytest
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from openmm import unit


@pytest.mark.parametrize(
    "forcefield",
    [
        pytest.param("de-force-1.0.0.offxml", id="No constraints"),
        pytest.param("de-force_unconstrained-1.0.0.offxml", id="Constraints"),
    ],
)
def test_load_de_ff():
    """
    Load the DE FF and create an OpenMM system.
    """

    ff = ForceField("de-force-1.0.0.offxml", load_plugins=True)
    ethanol = Molecule.from_smiles("CCO")

    system = ff.create_openmm_system(topology=ethanol.to_topology())

    forces = {force.__class__.__name__: force for force in system.getForces()}

    # make sure scaled 1-4 interactions are present
    assert "CustomBondForce" in forces
    assert "CustomNonbondedForce" in forces

    # make sure the sigma and epsilon are dummy parameters
    nonbonded: openmm.NonbondedForce = forces["NonbondedForce"]
    for i in range(nonbonded.getNumParticles()):
        _, sigma, epsilon = nonbonded.getParticleParameters(i)
        assert sigma == 1 * unit.angstroms
        assert epsilon == 0 * unit.kilocalorie_per_mole

    # check the custom force parameters are non-zero
    custom_force: openmm.CustomNonbondedForce = forces["CustomNonbondedForce"]
    for i in range(custom_force.getNumParticles()):
        epsilon, _ = custom_force.getParticleParameters(i)
        # no units with our custom force
        assert epsilon > 0
