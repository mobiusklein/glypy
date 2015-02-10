__all__ = ["composition", "io", "structure", "utils"]

from .structure.named_structures import monosaccharides, glycans
from .structure import Glycan, Monosaccharide, Substituent, Link
from .composition import Composition

import pkg_resources

pkg_resources.declare_namespace('pygly2')