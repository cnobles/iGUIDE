from semantic_version import Version
from pkg_resources import get_distribution

__version__ = str(Version.coerce(get_distribution('iguide').version))
