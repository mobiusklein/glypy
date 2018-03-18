from matplotlib.path import Path
from matplotlib.textpath import TextPath
import matplotlib.patches as patches
from matplotlib.transforms import IdentityTransform
import matplotlib
from matplotlib.colors import rgb2hex

from glypy.structure import Modification, Stem, SuperClass
from glypy.utils import where
from glypy.utils.enum import Enum
from glypy.io.nomenclature import identity
from glypy.algorithms import similarity


from .cfg_symbols import (
    draw_star, draw_square, draw_circle,
    draw_diamond, draw_text, draw_bisected_square,
    draw_vertical_bisected_diamond, draw_horizontal_bisected_diamond,
    draw_generic)


def rgb255(a, b, c):
    return rgb2hex((a / 255., b / 255., c / 255.))


class Colors(Enum):
    white = rgb255(255, 255, 255)
    blue = rgb255(0, 144, 188)
    green = rgb255(0, 166, 188)
    yellow = rgb255(255, 212, 0)
    light_blue = rgb255(143, 204, 233)
    pink = rgb255(246, 158, 161)
    purple = rgb255(165, 67, 153)
    brown = rgb255(161, 122, 77)
    orange = rgb255(244, 121, 32)
    red = rgb255(237, 28, 36)


def test_color_pallette(ax):
    for i, kv in enumerate(Colors):
        name, color = kv[0], kv[1]
        if name == "?":
            print(kv)
        draw_square(ax, i, 0, color)
        ax.text(i, -0.5, name, ha='center')
    ax.set_xlim(-1, i + 1)
    ax.set_ylim(-1, 1)
