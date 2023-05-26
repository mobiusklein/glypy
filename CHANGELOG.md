# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][Keep a Changelog] and this project adheres to [Semantic Versioning][Semantic Versioning].

## [1.0.9] - 2023-05-26

### Added

- `glypy.io.gnome` - An interface to the GNOme subsumption hierarchy.

### Changed

### Deprecated

### Removed

### Fixed
- Fix type conversion error in Py3.10 in `GlycanComposition.serialize`

### Security


## [1.0.8] - 2023-03-18

### Added

### Changed

### Deprecated

### Removed

### Fixed
- `GlycanComposition.mass` no longer ignores options which do not match the first invocation of that
  method on the same object (GH #25)

### Security



## [1.0.7] - 2022-07-22

### Added
- Added two new parameters to `plot` and `draw_tree`, `pad` and `pad_factor` to control
  the padding around a drawn glycan in stand-alone mode. Contributed by @michaelmarty (GH #24)

### Changed

### Deprecated

### Removed

### Fixed

### Security


## [1.0.6] - 2022-01-15

### Added
- Add specialized path for comparing two fully initialized `HashableGlycanComposition` instances
- Updated `GlycoRDF` to use OWL instead of cached turtle

### Changed

### Deprecated

### Removed

### Fixed

### Security


## [1.0.5]

### Changed
1. If a user passes `Pent` to an IUPACLite parser, it will be interpreted as `Pen`.


## [1.0.4] - 2021-10-29

### Added
1. Added support for writing UND sections using `glypy.io.glycoct`. The previous implementation
   is still available explicitly in `glypy.io.glycoct.OrderRespectingGlycoCTWriter`, the new
   behavior is implemented in `UNDOrderRespectingGlycoCTWriter`, which is now the value of
   `glypy.io.glycoct.GlycoCTWriter`.

### Changed
1. Added a C accelerator for `GlycanComposition.serialize`

### Deprecated

### Removed

### Fixed
1. Tolerate empty glycan composition strings in the first pass of the tokenizer again

### Security

## [1.0.3] - 2021-09-02

### Added
1. Added simple validation logic when handling malformed glycan composition string

### Changed

### Deprecated

### Removed

### Fixed

### Security

## [1.0.2] - 2021-06-21

### Added
1. Added C-extensions for core components of the `GlycanComposition` class

### Changed

### Deprecated

### Removed

### Fixed

### Security


## [1.0.1] - 2021-05-04

### Added
1. Added support for more monomer names to `glypy.io.byonic`.

### Changed
1. `FrozenMonosaccharideResidue.from_iupac_lite` now accounts for when the parsed string is different
   from the canonical name (e.g. `"NeuAc"` vs `"Neu5Ac"`). The returned objects are now guaranteed to be
   identical and both cached properly.
2. Improved documentation of extra line formats and IUPAClite.
3. Added support for "O"-prefixed substituents when parsing IUPAC. Does not alter the
   chemical interpretation, but resolves the name without the "O".

### Deprecated

### Removed

### Fixed
1. Fixed several pre-canned queries in `glypy.io.glyspace` which relied on glycan objects having
   an is-a relationship with `glycan:saccharide`. This appears to no longer be the case.

### Security


## [0.12.7] - 2021-02-02

### Added
* Added `glypy.io.byonic` to parse and write Byonic glycan compositions.

### Changed
* `Oct` is no longer automatically considered a `Mannose`
* `glypy.io.iupac` now accepts unknown positions for substituents for `Neu` bases.

### Deprecated

### Removed

### Fixed
* `glyspace` reading commands do not fail if a structure fails to parse due to a missing structure definition or reference

### Security

---

## [Released]

## [0.12.6] - 2020-06-21

### Added

* Added a CHANGELOG file
* `EnumValue` instances now support the XOR `^` operator
* WURCS encoded compositions which have edges that connect everything to everything are now treated as un-localized modifications
* `FrozenGlycanComposition` and its descendants now explicitly have cache field management as part of their API, see `_validate` and `_invalidate`
* Added a draft CFG text parser to `glypy.io.cfg` based heavily on the IUPAC parser.

### Changed

* `GlycanComposition` instances can be compared to any `Mapping`, not just other `GlycanComposition` instances.

### Fixed

* Ketoses have their masses properly calculated [GH18][https://github.com/mobiusklein/glypy/issues/18]
* Added a shim for deprecated `matplotlib.path.get_path_extents` when it is no longer available in recent versions of `matplotlib`
* The composition for phospho-ethanolamine has been corrected
* Fucose and other monosaccharides which are positioned to the side in balanced CFG/SNFG plots are now positioned properly using the topological layout plots

---

<!-- Links -->
[Keep a Changelog]: https://keepachangelog.com/
[Semantic Versioning]: https://semver.org/

<!-- Versions -->
[Unreleased]: https://github.com/mobiusklein/glypy/compare/v1.0.9...HEAD
[Released]: https://github.com/mobiusklein/glypy/releases
[0.12.6]: https://github.com/mobiusklein/glypy/releases/v0.12.6
[0.12.7]: https://github.com/mobiusklein/glypy/releases/v0.12.7
[1.0.0]: https://github.com/mobiusklein/glypy/releases/v1.0.0
[1.0.1]: https://github.com/mobiusklein/glypy/releases/v1.0.1
[1.0.2]: https://github.com/mobiusklein/glypy/releases/v1.0.2
[1.0.3]: https://github.com/mobiusklein/glypy/releases/v1.0.3
[1.0.4]: https://github.com/mobiusklein/glypy/releases/v1.0.4
[1.0.5]: https://github.com/mobiusklein/glypy/releases/v1.0.5
[1.0.6]: https://github.com/mobiusklein/glypy/releases/v1.0.6
[1.0.7]: https://github.com/mobiusklein/glypy/releases/v1.0.7
[1.0.8]: https://github.com/mobiusklein/glypy/releases/v1.0.8
[1.0.9]: https://github.com/mobiusklein/glypy/releases/v1.0.9