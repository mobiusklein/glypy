import pkg_resources
__all__ = ["composition"]
__package__ = "pygly2.composition"

import composition
from .composition import Composition, calculate_mass

pkg_resources.declare_namespace('pygly2.composition')
