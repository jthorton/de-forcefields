from . import _version

versions = _version.get_versions()
__version__ = versions["version"]
__git_revision = versions["full-revisionid"]
del versions
