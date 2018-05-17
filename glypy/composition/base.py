
class ChemicalCompositionError(Exception):
    pass


def composition_factory(*args, **kwargs):  # pragma: nocover
    '''
    Late-binding factory function. Used to reconstitute Composition objects
    that have been pickled, where the C-extension version may not be available
    on deserialization.
    '''
    try:
        from .composition import CComposition as Composition
    except ImportError:
        from .composition import PComposition as Composition
    return Composition()


def formula(composition):
    """Convert a :class:`~.Composition` into a :class:`str` form that can
    be parsed.

    Parameters
    ----------
    composition: :class:`~.Composition`
        A chemical composition

    Returns
    -------
    :class:`str`
    """
    return ''.join("%s%d" % (k, v) for k, v in sorted(composition.items()))
