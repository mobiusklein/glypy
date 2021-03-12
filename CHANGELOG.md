# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][Keep a Changelog] and this project adheres to [Semantic Versioning][Semantic Versioning].

## [Unreleased]

### Added


### Changed
* `FrozenMonosaccharideResidue.from_iupac_lite` now accounts for when the parsed string is different
  from the canonical name (e.g. `"NeuAc"` vs `"Neu5Ac"`). The returned objects are now guaranteed to be
  identical and both cached properly.

### Deprecated

### Removed

### Fixed

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
[Unreleased]: https://github.com/mobiusklein/glypy/compare/v0.12.7...HEAD
[Released]: https://github.com/mobiusklein/glypy/releases
[0.12.6]: https://github.com/mobiusklein/glypy/releases/v0.12.6
[0.12.7]: https://github.com/mobiusklein/glypy/releases/v0.12.7