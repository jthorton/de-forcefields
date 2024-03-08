import os.path

import pytest
from openff.toolkit.topology import Molecule


def get_data(relative_path: str) -> str:
    """
    Get the file path to some data in the package for testing.

    Args:
        relative_path: The relative path of the file to be loaded from deforcefields/data

    Returns:
        The absolute path to the requested file in deforcefields/data.
    """
    from pkg_resources import resource_filename

    file_name = resource_filename("deforcefields", os.path.join("data", relative_path))
    if not os.path.exists(file_name):
        raise ValueError(
            f"{relative_path} does not exist. If you have just added it, you'll have to re-install."
        )
    return file_name


@pytest.fixture()
def ethanol_with_charges() -> Molecule:
    """
    Return and OpenFF Molecule model of ethanol with `am1bccelf10` charges computed with openeye
    """
    ethanol = Molecule.from_file(get_data("ethanol.sdf"))
    return ethanol


@pytest.fixture()
def water() -> Molecule:
    """
    Return an OpenFF Molecule model of water with no charges.
    """
    return Molecule.from_file(get_data("water.sdf"))
