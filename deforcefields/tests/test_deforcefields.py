"""
Test loading DE-Force fields via the plugin interface through the toolkit.
"""
import pytest
from openff.toolkit.typing.engines.smirnoff import ForceField
from openmm import openmm, unit
from smirnoff_plugins.utilities.openmm import (
    evaluate_energy,
    evaluate_water_energy_at_distances,
)


@pytest.mark.parametrize(
    "forcefield",
    [
        pytest.param("de-force-1.0.0.offxml", id="No constraints"),
        pytest.param("de-force_unconstrained-1.0.0.offxml", id="Constraints"),
    ],
)
def test_load_de_ff(forcefield, ethanol_with_charges):
    """
    Load the DE FF and create an OpenMM system.
    """

    ff = ForceField(forcefield, load_plugins=True)

    system = ff.create_openmm_system(topology=ethanol_with_charges.to_topology())

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


@pytest.mark.parametrize(
    "forcefield, ref_energy",
    [
        pytest.param("de-force-1.0.0.offxml", 13.601144438830156, id="No constraints"),
        pytest.param(
            "de-force_unconstrained-1.0.0.offxml", 13.605201859835375, id="Constraints"
        ),
    ],
)
def test_energy_no_sites(forcefield, ref_energy, ethanol_with_charges):
    """
    Test calculating the single point energy of ethanol using constrained and unconstrained DE-FF with pre-computed
    partial charges from openeye.
    """

    ff = ForceField(forcefield, load_plugins=True)

    off_top = ethanol_with_charges.to_topology()
    openmm_top = off_top.to_openmm()
    system = ff.create_openmm_system(
        topology=off_top, charge_from_molecules=[ethanol_with_charges]
    )
    energy = evaluate_energy(
        system=system, topology=openmm_top, positions=ethanol_with_charges.conformers[0]
    )
    assert ref_energy == pytest.approx(energy)


def test_energy_sites():
    """
    Test calculating the energy for a system with two waters with virtual sites at set distances.
    """

    ff = ForceField("de-force-1.0.0.offxml", load_plugins=True)
    energies = evaluate_water_energy_at_distances(force_field=ff, distances=[2, 3, 4])
    ref_values = [1005.0846252441406, 44.696786403656006, 10.453390896320343]
    for i, energy in enumerate(energies):
        assert energy == pytest.approx(ref_values[i])
