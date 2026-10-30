# Releasing edna2obis

## Where the version lives

The version is defined in exactly three places, which must agree:

| Place | Example |
|---|---|
| `src/_version.py` | `__version__ = "3.3.0"` |
| Git tag on `main` | `v3.3.0` |
| `CHANGELOG.md` section | `## [3.3.0] - 2026-11-15` |

Directory names do not carry a version. GitHub Releases are the place to find the latest version.

## Version numbers

We follow [Semantic Versioning](https://semver.org/): `MAJOR.MINOR.PATCH`.

- **MAJOR**: a user has to change something to keep working. Examples: config keys renamed or removed, input/output paths moved, output columns renamed or removed.
- **MINOR**: new features or options that are off by default or backwards compatible. Examples: a new config option, a new optional output file.
- **PATCH**: bug fixes and documentation only.

Pre-releases use a suffix, e.g. `v4.0.0-rc.1`; they are published as GitHub pre-releases.

## Day to day

Add a line under `## [Unreleased]` in `CHANGELOG.md` in the same PR as the change, in the matching group (`Added`, `Changed`, `Deprecated`, `Removed`, `Fixed`). Mark breaking changes with **Breaking:**.

## Making a release

1. On `develop`, rename `## [Unreleased]` in `CHANGELOG.md` to `## [X.Y.Z] - YYYY-MM-DD`, add a fresh empty `## [Unreleased]` above it, and update the compare links at the bottom.
2. Set `__version__ = "X.Y.Z"` in `src/_version.py`.
3. Open a PR from `develop` to `main`, titled `Release vX.Y.Z`, and merge it.
4. Tag the merge commit on `main` and push the tag:

   ```bash
   git switch main && git pull
   git tag -a vX.Y.Z -m "edna2obis vX.Y.Z"
   git push origin vX.Y.Z
   ```

5. The [release workflow](../.github/workflows/release.yml) checks that the tag matches `src/_version.py`, then creates the GitHub Release using that version's `CHANGELOG.md` section as the notes. If the check fails, fix the files, delete the tag (`git push --delete origin vX.Y.Z`), and tag again.

## Legacy tags

`v1.0.0` is the notebook workflow (`legacy-v1` branch). `v3.0.0-alpha` through `v3.0.2-alpha` were the first script-based versions; there was no 2.x release. Versions up to 3.2.0 used `src-v3/`, `raw-v3/` and `processed-v3/` directories.
