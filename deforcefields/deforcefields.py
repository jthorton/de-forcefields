"""
Entry point function to add extra installed force fields
"""

from openff.utilities.utilities import get_data_dir_path


def get_forcefield_paths() -> list[str]:
    """
    Returns:
         The paths to the directories containing force field files.
    """

    return [get_data_dir_path(relative_path="offxml", package_name="deforcefields")]
