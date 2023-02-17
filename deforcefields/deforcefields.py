"""
Entry point function to add extra installed force fields
"""


from typing import List

from pkg_resources import resource_filename


def get_forcefield_paths() -> List[str]:
    """
    Returns:
         The paths to the directories containing force field files.
    """

    return [resource_filename("deforcefields", "offxml")]
