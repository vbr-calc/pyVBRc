# v0.3.1

This bugfix release pins unyt < 3.0

## bug fixes
- unyt < 3.0



# v0.3.0

This release adds some isotropic volume averaging via the new `IsotropicMixture` class

## New Features
- added `IsotropicMixture` class to do isotropic volume averaging, see examples.

## Other changes
- `anistotropic` module renamed to `materials`
- switch to pre-commit.ci for linting

# v0.2.0

Main change is adding support for unit-ful arrays via the unyt package.

## New Features

- units support: when loading a VBRc .mat file, arrays will be loaded as unyt arrays if possible (requires VBRc version >= 1.0.0)
- nicely formatted logging with the pyVBRc logger
- anisotropy calculations for aligned inclusions (experimental, consider in beta form)
- documentation updates (scripts are now notebooks in examples/)

## Bug Fixes

None

## Other changes
- code coverage reporting enabled for the repository
- style checks now done by precommit.ci bot
- switch to pyproject.toml build

# v0.1.1

Initial release of pyVBRc, still experimental

## New Features

* python interace for VBRc .mat files
* basic examples
