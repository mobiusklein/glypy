'''
Glyconnect
----------

A simple dialect of the Glyconnect/GlycoMod glycan composition notation.
'''
import re

from typing import Dict, Union, Optional, List, Tuple

from urllib.parse import quote
from dataclasses import dataclass, field
from typing import List, Optional, Type, Generic, TypeVar

from glypy.structure.glycan_composition import (
    FrozenGlycanComposition,
    FrozenMonosaccharideResidue,
    SubstituentResidue)
from glypy.structure.glycan import Glycan


try:
    import requests
except ImportError:
    requests = None

#: The set of defined symbols and their mappings.
defined_symbols: Dict[str, Union[SubstituentResidue, FrozenMonosaccharideResidue]] = {
    "Hex": FrozenMonosaccharideResidue.from_iupac_lite("Hex"),
    "HexNAc": FrozenMonosaccharideResidue.from_iupac_lite('HexNAc'),
    "dHex": FrozenMonosaccharideResidue.from_iupac_lite('dHex'),
    "NeuAc": FrozenMonosaccharideResidue.from_iupac_lite("NeuAc"),
    "NeuGc": FrozenMonosaccharideResidue.from_iupac_lite("NeuGc"),
    "S": SubstituentResidue("sulfate"),
    "Su": SubstituentResidue("sulfate"),
    "Sulpho": SubstituentResidue("sulfate"),
    "P": SubstituentResidue("phosphate"),
    "Ph": SubstituentResidue("phosphate"),
    "Phospho": SubstituentResidue("phosphate"),
    "Xyl": FrozenMonosaccharideResidue.from_iupac_lite("Xyl"),
    "HexA": FrozenMonosaccharideResidue.from_iupac_lite("HexA"),
    "Pent": FrozenMonosaccharideResidue.from_iupac_lite("Pen"),
    "Kdn": FrozenMonosaccharideResidue.from_iupac_lite("Kdn"),
}


def _invert_mapping(table: Dict[str, Union[SubstituentResidue, FrozenMonosaccharideResidue]]) -> Dict[Union[SubstituentResidue, FrozenMonosaccharideResidue], str]:
    inverted = {}
    for k, v in table.items():
        if v in inverted:
            if len(k) > len(inverted[v]):
                continue
        inverted[v] = k
    return inverted


monosaccharide_to_symbol = _invert_mapping(defined_symbols)


def _generate_pattern(symbols: List[str]) -> re.Pattern:
    symbols = sorted(symbols, key=len, reverse=True)
    return re.compile(f"({'|'.join(symbols)})(\d+?)")


tokenizer = re.compile(r"([^:\s]+):(\d+)")
undelimited_tokenizer = _generate_pattern(defined_symbols)


def loads(string):
    '''Parse a GlyConnect glycan composition into a :class:`~.FrozenGlycanComposition`

    Parameters
    ----------
    string: str
        The string to parse

    Returns
    -------
    :class:`~.FrozenGlycanComposition`

    Raises
    ------
    :class:`KeyError`: Raised if a key isn't defined by the GlyConnect dialect
    '''
    tokens = tokenizer.findall(string)
    if not tokens:
        tokens = undelimited_tokenizer.findall(string)
    gc = FrozenGlycanComposition()
    for mono, count in tokens:
        mono = defined_symbols[mono]
        count = int(count)
        gc[mono] += count
    return gc


def dumps(composition):
    '''Encode :class:`~.GlycanComposition` or :class:`~.Glycan` into the GlyConnect
    glycan composition text format.

    Parameters
    ----------
    composition: :class:`~.GlycanComposition` or :class:`~.Glycan`
        The structure to format

    Returns
    -------
    :class:`str`

    Raises
    ------
    :class:`KeyError`: Raised if a key isn't defined by the GlyConnect Compozitor dialect
    '''
    if isinstance(composition, Glycan):
        composition = FrozenGlycanComposition.from_glycan(composition)
    tokens = []
    for key, value in composition.items():
        key = monosaccharide_to_symbol[key]
        tokens.append("%s:%d" % (key, value))
    return ' '.join(tokens)


API_SERVER = "https://glyconnect.expasy.org/api"


def from_glytoucan_id(glytoucan_id):
    response = requests.post(
        f"{API_SERVER}/structures/search/glytoucan",
        data={"glytoucan_id": glytoucan_id})
    response.raise_for_status()
    data = response.json()
    return data


@dataclass
class RecordBase:

    @classmethod
    def from_dict(cls, data):
        return cls(**data)


@dataclass
class TaxonomyRecord(RecordBase):
    id: int
    taxonomy_id: str
    common_name: Optional[str] = None
    species: Optional[str] = None


@dataclass
class UniprotProteinAccessionRecord(RecordBase):
    uniprot_acc: str
    uniprot_id: Optional[str] = None
    glygen: Optional[str] = None
    nextprot: Optional[str] = None
    genecards: Optional[str] = None
    glycodomain: Optional[str] = None


@dataclass
class ProteinRecord(RecordBase):
    id: int
    name: str
    taxonomy: TaxonomyRecord
    uniprots: List[UniprotProteinAccessionRecord]

    @classmethod
    def from_dict(cls, data: dict):
        tax = data.get("taxonomy")
        if tax:
            tax = TaxonomyRecord.from_dict(tax)
        uniprots = list(map(UniprotProteinAccessionRecord.from_dict,
                            data.get("uniprots", [])))
        return cls(data['id'], data['name'], tax, uniprots)


@dataclass
class SourceRecord(RecordBase):
    type: str
    name: str
    id: int
    ref: Optional[str] = None
    ontology: Optional[str] = None
    brenda_id: Optional[str] = None


@dataclass
class Source(RecordBase):
    source: List[SourceRecord]
    taxons: List[TaxonomyRecord]

    @classmethod
    def from_dict(cls, data: dict):
        source = [SourceRecord.from_dict(x) for x in data.get("source", [])]
        taxons = [TaxonomyRecord.from_dict(x) for x in data.get("taxons", [])]
        return cls(source, taxons)


@dataclass
class CellLine(RecordBase):
    cellosaurus_id: str
    id: int
    is_problematic: bool
    name: str


@dataclass
class Disease:
    id: int
    name: str
    do_id: Optional[str] = None
    taxons: List[TaxonomyRecord] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: dict):
        taxons = [TaxonomyRecord.from_dict(x) for x in data.get("taxons", [])]
        return cls(data['id'], data['name'], data.get('do_id'), taxons)


@dataclass
class CompositionRecord(RecordBase):
    format_byonic: str
    format_condensed: str
    format_glyconnect: str
    format_numeric: str
    id: int
    mass: float
    mass_monoisotopic: float
    reviewed: bool
    glytoucan_id: Optional[str] = None

    def parse(self):
        return loads(self.format_glyconnect)


@dataclass
class StructureRecord(RecordBase):
    glycan_core: str
    glycan_type: str
    has_image: bool
    id: int
    is_undefined: bool
    reviewed: bool
    glytoucan_id: Optional[str] = None


@dataclass
class CompozitorGlycan(RecordBase):
    composition: CompositionRecord
    structure: StructureRecord
    taxonomy: Optional[TaxonomyRecord]
    protein: Optional[ProteinRecord]

    @classmethod
    def from_dict(cls, data: dict):
        comp = CompositionRecord.from_dict(data['composition'])
        struct = StructureRecord.from_dict(data['structure'])
        protein = ProteinRecord.from_dict(data['protein'])
        taxonomy = TaxonomyRecord.from_dict(data['taxonomy'])
        return cls(comp, struct, protein, taxonomy)


T = TypeVar("T", bound=RecordBase)


@dataclass
class APICollectionProperty(Generic[T]):
    url: str
    record_type: Type[T]

    def __get__(self, obj, objtype=None) -> List[T]:
        if obj is None:
            return self
        result = obj._cache.get(self.url)
        if result is not None:
            return result
        resp = requests.get(self.url)
        resp.raise_for_status()
        data = resp.json()
        result = [self.record_type.from_dict(d) for d in data]
        obj._cache[self.url] = result
        return result

    def __delete__(self, obj):
        del obj._cache[self.url]


@dataclass
class Compozitor:
    _cache: dict = field(default_factory=dict, repr=False)

    proteins = APICollectionProperty(
        f"{API_SERVER}/proteins-all",
        ProteinRecord
    )

    sources = APICollectionProperty(
        f"{API_SERVER}/sources-all",
        Source
    )

    cell_lines = APICollectionProperty(
        f"{API_SERVER}/cell_lines-all",
        CellLine
    )

    diseases = APICollectionProperty(
        f"{API_SERVER}/diseases-all",
        Disease
    )

    def query(self, taxonomy: Optional[str]=None, cell_line: Optional[str]=None,
              protein: Optional[str]=None, disease: Optional[str]=None):
        params = {}
        if taxonomy:
            params['taxonomy'] = (taxonomy)
        if cell_line:
            params['cell_line'] = (cell_line)
        if protein:
            params['protein'] = (protein)
        if disease:
            params['disease'] = (disease)
        resp = requests.get(f"{API_SERVER}/glycosylations", params)
        resp.raise_for_status()
        data = resp.json()
        if isinstance(data, list):
            raise ValueError("Malformed query or invalid response")
        results = []
        if data['results']:
            for res in data['results']:
                results.append(CompozitorGlycan.from_dict(res))
        return results


client = Compozitor()
query = client.query
