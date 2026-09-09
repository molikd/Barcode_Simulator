# Archival and Preservation Strategy

Barcode_Simulator is a dependency-free POSIX AWK + shell project, so
preservation is straightforward: the git history plus this document is the
whole artifact. Following the ReadMixer precedent reviewed in
[issue #4](https://github.com/molikd/Barcode_Simulator/issues/4):

## Locations

1. **GitHub** (`https://github.com/molikd/Barcode_Simulator`) — primary
   development, issue tracking, and CI. Not itself a long-term archive.
2. **Zenodo (planned)** — mint a DOI for the current state via the
   GitHub–Zenodo integration, then record the DOI in `CITATION.cff`.
   Re-mint per release thereafter.
3. **Software Heritage (planned)** — submit the repository URL for
   universal source-code archival alongside the Zenodo deposit.

## What is preserved

- Full git history (code, tests, docs).
- `CITATION.cff` with the preferred citation to
  [Molik, Pfrender, and Emrich (2020)](https://doi.org/10.3390/mps3010022).
- `Barcode_Simulator_SBOM.json`: CycloneDX SBOM recording that there are
  no runtime library dependencies beyond POSIX AWK, a POSIX shell, and
  (only for `.gz` features) `gzip`.
- `CHANGES.md` release notes and `TESTING.md` for reproducibility.

## Maintenance

- On each release: tag the version, update `CHANGES.md`, verify the
  Zenodo deposit, refresh the DOI in `CITATION.cff`.
- Annually: verify all archival locations resolve; update contact info.
