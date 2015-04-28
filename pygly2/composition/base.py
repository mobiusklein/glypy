
class ChemicalCompositionError(Exception):
    pass


def composition_factory(*args, **kwargs):  # pragma: nocover
    '''
    Late-binding factory function. Used to reconstitute Composition objects
    that have been pickled, where the C-extension version may not be available
    on deserialization.
    '''
    try:
        from composition import CComposition as Composition
    except:
        from composition import PComposition as Composition
    return Composition()
