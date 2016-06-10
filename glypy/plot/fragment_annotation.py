from matplotlib import patches
from matplotlib import transforms as mtransforms
from matplotlib.textpath import TextPath
import glypy.plot.cfg_symbols
reload(glypy.plot.cfg_symbols)
from glypy.plot.cfg_symbols import line_to
from glypy.plot.draw_tree import get_link_pair


def isclose(a, b):
    return abs(a - b) < 1e-4

rot90 = mtransforms.Affine2D().rotate_deg(90)


class BondCleavageArtist(object):

    def __init__(self, draw_tree, fragment, ax):
        self.draw_tree = draw_tree
        self.fragment = fragment
        self.ax = ax
        self.artists = []
        self.patches = []
        self.handle_glycosidic()
        self.handle_crossring()
        self.propagate_transform()

    def handle_glycosidic(self):
        for link_id, frag_type in self.fragment.link_ids.items():
            parent, child = self.draw_tree.get_link_pair(link_id)
            if frag_type[1] in "YZ":
                reducing = True
            else:
                reducing = False
            gba = GlycosidicBondCleavageStroke(self.draw_tree, parent, child, reducing,
                                               ax=self.ax, fragment=self.fragment)
            gba.draw()
            self.artists.append(gba)
            self.patches.extend(gba.patches)

    def handle_crossring(self):
        for crossring_id, frag_type in self.fragment.crossring_cleavages.items():
            node = self.draw_tree.get(crossring_id)
            if frag_type[1] == "X":
                reducing = True
            else:
                reducing = False
            cra = CrossRingBondCleavageStroke(
                self.draw_tree, node, reducing, ax, self.fragment)
            cra.draw()
            self.artists.append(cra)
            self.patches.extend(cra.patches)

    def propagate_transform(self):
        transform = self.draw_tree.transform
        if transform is None:
            return
        transform = transform + self.ax.transData
        for patch in self.patches:
            patch.set_transform(transform)


class CleavageStrokeArtistBase(object):

    def __init__(self, draw_tree, fragment, reducing, ax, **options):
        self.draw_tree = draw_tree
        self.ax = ax
        self.reducing = reducing
        self.fragment = fragment
        self.patches = []
        self.x_coords = []
        self.y_coords = []
        self.patch_dict = draw_tree.data['patches']
        self.options = options
        self.edge_size = None

    def save_patch(self, key, value):
        self.patches.append(value)
        if key in self.patch_dict:
            self.patch_dict[key].append([value])
        else:
            self.patch_dict[key] = [[value]]

    def textpath(self, x, y, text, line_weight=0.5, **kwargs):
        fs = kwargs.get("fontsize", 2) * .5
        t_path = TextPath((x, y), s=text, size=fs)
        patch = patches.PathPatch(t_path, facecolor="black", lw=line_weight / 20., zorder=4)
        a = self.ax.add_patch(patch)
        return a


class CrossRingBondCleavageStroke(CleavageStrokeArtistBase):

    def __init__(self, draw_tree, node, reducing, ax, fragment, **options):
        super(CrossRingBondCleavageStroke, self).__init__(
            draw_tree, fragment, reducing, ax, **options)
        self.node = node
        self.mode = None

    def draw(self):
        self.compute_positions()
        self.draw_main_stroke()
        self.draw_edge_stroke()

    def compute_positions(self):
        px, py = self.node.x, self.node.y
        if self.node.children:
            nx, ny = -1, -1
            for next_node in self.node.children:
                cx, cy = next_node.x, next_node.y
                if cy > ny:
                    nx, ny = cx, cy
            self.mode = 'child'
            cx, cy = nx, ny
            xcenter = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
            ycenter = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
            xlength = (xcenter - min(px, cx)) / 2.
            ylength = (ycenter - min(py, cy)) / 2.
            if py < cy:
                lx2, ly2 = self.node.x + xlength, self.node.y + ylength
                lx1, ly1 = self.node.x - xlength, self.node.y - ylength
            else:
                lx1, ly1 = self.node.x - xlength, self.node.y + ylength
                lx2, ly2 = self.node.x + xlength, self.node.y - ylength
            edge_size = abs(ly1 - ly2)

        elif self.node.parent:
            nx, ny = self.parent
            self.mode = 'parent'

            cx, cy = px, py
            px, py = nx, ny

            xcenter = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
            ycenter = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
            xlength = (xcenter - min(px, cx)) / 2.
            ylength = (ycenter - min(py, cy)) / 2.
            if py < cy:
                lx2, ly2 = xcenter + xlength, ycenter + ylength
                lx1, ly1 = xcenter - xlength, ycenter - ylength
            else:
                lx1, ly1 = xcenter - xlength, ycenter + ylength
                lx2, ly2 = xcenter + xlength, ycenter - ylength
            edge_size = abs(ly1 - ly2)

        else:
            print "Orphaned Node. Cross-ring fragment drawing may be bad"

        self.x_coords = [lx1, lx2]
        self.y_coords = [ly1, ly2]
        self.edge_size = edge_size

    def draw_main_stroke(self, **kwargs):
        graphic_options = self.options.copy()
        graphic_options.update(kwargs)

        color = graphic_options.get('color', 'red')
        lw = graphic_options.get('lw', 2)

        lx1, lx2 = self.x_coords
        ly1, ly2 = self.y_coords

        line = line_to(self.ax, lx1, ly1, lx2, ly2,
                       color=color, zorder=10, lw=lw)
        self.save_patch(self.fragment.name, line)
        line.set_gid(self.draw_tree.uuid + '-' + self.fragment.name)

    def draw_edge_stroke(self, **kwargs):
        graphic_options = self.options.copy()
        graphic_options.update(kwargs)

        color = graphic_options.get('color', 'red')
        lw = graphic_options.get('lw', 2)

        lx1, lx2 = self.x_coords
        ly1, ly2 = self.y_coords

        if self.reducing:
            line = line_to(self.ax, lx1, ly2, lx1, ly2 -
                           self.edge_size, color=color, zorder=3, lw=lw)
            line.set_gid(self.draw_tree.uuid + '-' +
                         self.fragment.name + "-direction")
            self.save_patch(self.fragment.name + "_direction", line)
        else:
            line = line_to(self.ax, lx2, ly2, lx2, ly2 +
                           self.edge_size, color=color, zorder=3, lw=lw)
            line.set_gid(self.draw_tree.uuid + '-' +
                         self.fragment.name + "-direction")
            self.save_patch(self.fragment.name + "_direction", line)


class GlycosidicBondCleavageStroke(CleavageStrokeArtistBase):

    def __init__(self, draw_tree, parent, child, reducing, ax, fragment, **options):
        super(GlycosidicBondCleavageStroke, self).__init__(
            draw_tree, fragment, reducing, ax, **options)

        self.parent = parent
        self.child = child

        self.mode = None

    def draw(self):
        self.compute_positions()
        self.draw_main_stroke()
        self.draw_edge_stroke()

    def compute_positions(self):
        parent, child = self.parent, self.child
        px, py = parent.x, parent.y
        cx, cy = child.x, child.y

        if isclose(py, cy) and not isclose(px, cx):
            mode = 'offsequence'
            center = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
            length = (center - min(px, cx))
            lx1, ly1 = center, py + length
            lx2, ly2 = center, py - length
            edge_size = abs(ly1 - ly2) / 2.

        elif not isclose(py, cy) and isclose(px, cx):
            mode = "sequential"
            center = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
            length = center - min(py, cy)
            lx1, ly1 = px + length, center
            lx2, ly2 = px - length, center
            edge_size = abs(lx1 - lx2) / 2.
        else:
            mode = "branching"
            xcenter = (max(px, cx) - min(px, cx)) / 2. + min(px, cx)
            ycenter = (max(py, cy) - min(py, cy)) / 2. + min(py, cy)
            xlength = (xcenter - min(px, cx)) / 2.
            ylength = (ycenter - min(py, cy)) / 2.
            if py < cy:
                lx2, ly2 = xcenter + xlength, ycenter + ylength
                lx1, ly1 = xcenter - xlength, ycenter - ylength
            else:
                lx1, ly1 = xcenter - xlength, ycenter + ylength
                lx2, ly2 = xcenter + xlength, ycenter - ylength
            edge_size = abs(ly1 - ly2)

        self.mode = mode
        self.edge_size = edge_size
        self.x_coords = [lx1, lx2]
        self.y_coords = [ly1, ly2]

    def draw_main_stroke(self, **kwargs):
        graphic_options = self.options.copy()
        graphic_options.update(kwargs)

        color = graphic_options.get('color', 'red')
        lw = graphic_options.get('lw', 2)

        lx1, lx2 = self.x_coords
        ly1, ly2 = self.y_coords

        line = line_to(self.ax, lx1, ly1, lx2, ly2,
                       color=color, zorder=3, lw=lw)
        self.save_patch(self.fragment.name, line)
        line.set_gid(self.draw_tree.uuid + '-' + self.fragment.name)

    def draw_edge_stroke(self, **kwargs):
        graphic_options = self.options.copy()
        graphic_options.update(kwargs)

        color = graphic_options.get('color', 'red')
        lw = graphic_options.get('lw', 2)

        lx1, lx2 = self.x_coords
        ly1, ly2 = self.y_coords
        # if self.mode == "offsequence":
        #    lx1, lx2, _ = rot90.get_matrix().dot([lx1, lx2, 0])
        #    ly1, ly2, _ = rot90.get_matrix().dot([ly1, ly2, 0])
        if self.reducing:
            line = line_to(self.ax, lx1, ly2, lx1, ly2 -
                           self.edge_size, color=color, zorder=3, lw=lw)
            line.set_gid(self.draw_tree.uuid + '-' +
                         self.fragment.name + "-direction")
            self.save_patch(self.fragment.name + "_direction", line)
        else:
            line = line_to(self.ax, lx2, ly2, lx2, ly2 +
                           self.edge_size, color=color, zorder=3, lw=lw)
            line.set_gid(self.draw_tree.uuid + '-' +
                         self.fragment.name + "-direction")
            self.save_patch(self.fragment.name + "_direction", line)
