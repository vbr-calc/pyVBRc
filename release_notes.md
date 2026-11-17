# v0.4.0

This is a maintenance release focused on packaging and infrastructure updates. There are no changes to pyVBRc functionality.

## Breaking changes
- the minimum supported python version is now 3.10 (PR #92, issue #91)

## Infrastructure and maintenance
- add testing for newer python versions (PR #88)
- initialize `uv` for development and update the developer setup; add dependabot config, an AI policy and a standalone `CONTRIBUTING.md`; remove unused templates and dev utilities; remove the deprecated license classifier (PR #92)
- update release workflows to use trusted publishing (PR #96)
