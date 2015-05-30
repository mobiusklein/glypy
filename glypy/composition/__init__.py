import pkg_resources
__all__ = ["composition"]
__package__ = "glypy.composition"

import composition
from .composition import Composition, calculate_mass

pkg_resources.declare_namespace('glypy.composition')
