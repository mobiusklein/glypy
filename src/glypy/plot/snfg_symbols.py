from glypy.structure import Stem, Modification
from glypy.utils.enum import Enum
from .cfg_symbols import CFGNomenclature


class SNFGNomenclature(CFGNomenclature):
    class ResidueColor(Enum):
        white = "#ffffff"
        blue = "#0180ff"
        green = "#01A650"
        yellow = "#ffd706"
        light_blue = "#91cce8"
        pink = "#fe87c5"
        purple = "#a21fff"
        brown = "#8e6319"
        orange = "#f37922"
        red = "#fc0200"
        generic = 'white'

    stem_to_color = {
        Stem.glc: ResidueColor.blue,
        Stem.man: ResidueColor.green,
        Stem.gal: ResidueColor.yellow,
        Stem.gul: ResidueColor.orange,
        Stem.alt: ResidueColor.pink,
        Stem.all: ResidueColor.purple,
        Stem.tal: ResidueColor.light_blue,
        Stem.ido: ResidueColor.brown,
        Stem.xyl: ResidueColor.orange,
        Stem.lyx: ResidueColor.yellow,
        Stem.x: ResidueColor.generic,
    }

    def residue_color(self, monosaccharide):
        '''
        Determine which color to use to represent `monosaccharide` under the CFG
        symbol nomenclature.

        Parameters
        ----------
        monosaccharide: |Monosaccharide|
            The residue to be rendered

        Returns
        -------
        ResidueColor.EnumValue

        '''
        if any(mod == Modification.a for p, mod in monosaccharide.modifications.items()):
            return self.resolve_acid_color(monosaccharide)
        if "hex" in [monosaccharide.superclass]:
            if any(mod == Modification.d for p, mod in monosaccharide.modifications.items()) and\
                    monosaccharide.stem == (Stem.gal,):
                return self.ResidueColor.red
        try:
            return self.stem_to_color[monosaccharide.stem[0]]
        except KeyError:
            return self.ResidueColor.generic

    def resolve_acid_color(self, monosaccharide):
        '''
        Resolve the special case in :func:`residue_color` for acidic residues
        '''
        if ('gro' in monosaccharide.stem) and ('gal' in monosaccharide.stem):
            if any(sub.name == 'n_acetyl' for p, sub in monosaccharide.substituents()):
                return self.ResidueColor.purple
            elif any(sub.name == 'n_glycolyl' for p, sub in monosaccharide.substituents()):
                return self.ResidueColor.light_blue
            elif any(sub.name == 'amino' for p, sub in monosaccharide.substituents()):
                return self.ResidueColor.brown
            else:
                return self.ResidueColor.green
        else:
            return self.stem_to_color[monosaccharide.stem[0]]


def test_color_pallette(ax):
    d = SNFGNomenclature()
    for i, kv in enumerate(SNFGNomenclature.ResidueColor):
        name, color = kv[0], kv[1]
        if name == "?":
            print(kv)
        d.draw_square(ax, i, 0, color)
        ax.text(i, -0.5, name, ha='center')
    ax.set_xlim(-1, i + 1)
    ax.set_ylim(-1, 1)
