from matplotlib import transforms as mtransform


class Affine2DProxy(object):
    def __init__(self, transform=None, text_transform=None, rotate_handler=lambda x: None):
        if transform is None:
            transform = mtransform.Affine2D()
        self.transform = transform
        if text_transform is None:
            text_transform = mtransform.Affine2D()
        self.text_transform = text_transform
        self.rotate_handler = rotate_handler

    def scale(self, x, y=None):
        self.transform.scale(x, y)
        self.text_transform.scale(x, y)

    def translate(self, x, y):
        self.transform.translate(x, y)
        self.text_transform.translate(x, y)

    def rotate_deg(self, deg):
        self.transform.rotate_deg(deg)
        self.text_transform.rotate_deg(deg)
        self.rotate_handler(deg)

    def rotate_deg_around(self, x, y, deg):
        self.transform.rotate_deg_around(x, y, deg)
        self.text_transform.rotate_deg_around(x, y, deg)
        self.rotate_handler(deg)

    def __add__(self, other):
        return Affine2DProxy(self.transform + other, self.rotate_handler)

    def __iadd__(self, other):
        self.transform += other
        return self

    def set(self, artist):
        artist.set_transform(self.transform)


class TransformDispatch(object):
    def __init__(self, transforms=None, rotate_handler=lambda x: None):
        if transforms is None:
            transforms = []
        self.transforms = transforms
        self.recalculated_transforms = []
        self.rotation_counter = 0
        self.rotate_handler = rotate_handler

    def add(self, transform, recalculated=False):
        if recalculated:
            self.recalculated_transforms.append(transform)
        else:
            self.transforms.append(transform)

    def scale(self, x, y=None):
        for transform in self.transforms:
            transform.scale(x, y)

        for transform in self.recalculated_transforms:
            transform.scale(x, y)

    def translate(self, x, y):
        for transform in self.transforms:
            transform.translate(x, y)

        for transform in self.recalculated_transforms:
            transform.translate(x, y)

    def rotate_deg(self, deg):
        self.rotation_counter -= deg
        for transform in self.transforms:
            transform.rotate_deg(deg)

        for transform in self.recalculated_transforms:
            transform.rotate_deg(deg)
        self.recalculated_transforms = self.rotate_handler(self.rotation_counter)

    def rotate_deg_around(self, x, y, deg):
        self.rotation_counter -= deg
        for transform in self.transforms:
            transform.rotate_deg_around(x, y, deg)

        for transform in self.recalculated_transforms:
            transform.rotate_deg_around(x, y, deg)
        self.recalculated_transforms = self.rotate_handler(self.rotation_counter)
