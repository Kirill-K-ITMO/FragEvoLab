# Third-party components and notices

Inventory date: 2026-10-30. Scope: FragEvoLab 3.1.4 source distribution and the inspected Linux environment. The project's [MIT License](LICENSE) covers its original code; it does not replace third-party licenses or data terms.

This archive contains application source, documentation, example/provenance tables and license-text snapshots. It does not bundle Python, Python package wheels, Vina, Open Babel, ChimeraX, or the RDKit SA-score implementation/model. Dependencies are installed separately. Scientific citations are in [docs/references.md](docs/references.md).

## Declared Python dependencies

Versions below are the pins in `requirements.txt` and `requirements-test.txt`, checked against installed package metadata. Full copied notices are under `licenses/`; [docs/dependency_inventory.json](docs/dependency_inventory.json) records source paths and hashes. Preserve the complete licenses and notices supplied by the actual dependency distribution when redistributing it.

| Component | Version | Role | License / notice |
|---|---|---|---|
| [RDKit](https://www.rdkit.org/) | 2026.3.6 | Chemistry, descriptors, fingerprints and drawings | BSD-3-Clause for core: [core notice](licenses/rdkit/RDKit-core-license.txt). The installed wheel also supplies a distinct MIT packaging notice: [LICENSE.md](licenses/rdkit/LICENSE.md). |
| [Meeko](https://github.com/forlilab/Meeko) | 0.8.0 | PDBQT parameterization | LGPL-2.1, as identified by the package: [license](licenses/meeko/LICENSE). |
| [Gemmi](https://github.com/project-gemmi/gemmi) | 0.7.5 | Dependency of the Meeko stack | MPL-2.0: [license](licenses/gemmi/LICENSE.txt). |
| [NumPy](https://numpy.org/) | 2.3.5 | Numerical dependency | BSD-3-Clause plus bundled-library notices in the installed distribution: [license](licenses/numpy/LICENSE.txt). |
| [SciPy](https://scipy.org/) | 1.17.0 | Numerical dependency | BSD-3-Clause plus bundled-library notices in the installed distribution: [license](licenses/scipy/LICENSE.txt). |
| [Pillow](https://python-pillow.github.io/) | 12.3.0 | Image dependency used by RDKit drawings | MIT-CMU: [license](licenses/Pillow/LICENSE). |
| [Rich](https://github.com/Textualize/rich) | 15.0.0 | Terminal interface and progress | MIT: [license](licenses/rich/LICENSE). |
| [filelock](https://github.com/tox-dev/filelock) | 4.0.12 | Cache/process locks | MIT in this installed version: [license](licenses/filelock/LICENSE). |
| [Jinja2](https://jinja.palletsprojects.com/) | 3.1.6 | HTML templating | BSD-3-Clause: [license](licenses/Jinja2/LICENSE.txt). |
| [pytest](https://docs.pytest.org/) | 9.1.1 | Tests only | MIT: [license](licenses/pytest/LICENSE). |

RDKit's contributed `SA_Score` implementation is imported from the installed RDKit package when available. It carries a separate copyright notice for Novartis Institutes for BioMedical Research Inc. (2013), with BSD three-clause terms. The [complete notice](licenses/rdkit/SA_Score-NOTICE.txt) is retained. Neither `sascorer.py` nor `fpscores.pkl.gz` is copied into this source distribution.

## Observed Python transitive dependencies

These versions describe the inspected environment, not additional pins or a cross-platform lockfile.

| Component | Observed version | License / retained notice |
|---|---|---|
| MarkupSafe | 3.0.4 | [BSD-3-Clause](licenses/MarkupSafe/LICENSE.txt) |
| markdown-it-py | 4.2.0 | [MIT](licenses/markdown-it-py/LICENSE), plus [upstream markdown-it notice](licenses/markdown-it-py/LICENSE.markdown-it) |
| mdurl | 0.1.2 | [MIT](licenses/mdurl/LICENSE) |
| Pygments | 2.21.0 | [BSD-2-Clause](licenses/Pygments/LICENSE) |
| packaging | 26.3 | [Apache-2.0 OR BSD-2-Clause](licenses/packaging/LICENSE); [Apache text](licenses/packaging/LICENSE.APACHE), [BSD text](licenses/packaging/LICENSE.BSD) |
| pluggy | 1.7.0 | [MIT](licenses/pluggy/LICENSE) |
| iniconfig | 2.3.1 | [MIT](licenses/iniconfig/LICENSE) |

On Windows, pytest additionally declares **colorama** (BSD-3-Clause). The [retained notice](licenses/colorama/LICENSE.txt) is from upstream 0.4.6, a reference version rather than an added pin; use the notice supplied by the actual installed version. Python 3.12 is the supported setup, so pytest's conditional exceptiongroup/tomli dependencies for Python below 3.11 are outside this configuration. No claim of a complete Windows binary license inventory is made.

## External programs and runtime

- **AutoDock Vina:** external executable, example-compatible version 1.2.7. Apache-2.0; [license snapshot](licenses/AutoDock-Vina-LICENSE.txt) from the [v1.2.7 source tag](https://github.com/ccsb-scripps/AutoDock-Vina/tree/v1.2.7). If redistributing Vina, include its applicable copyright/license/NOTICE files and any notices required by the distributed binary and its dependencies. This archive contains no Vina executable.
- **Open Babel:** optional external `obabel` fallback, invoked as a separate process. GPL version 2; [COPYING snapshot](licenses/OpenBabel-COPYING.txt) from [3.1.1](https://github.com/openbabel/openbabel/tree/openbabel-3-1-1), provided as a reference rather than an application version pin. Consult the actual installed version's terms and [upstream guidance](https://openbabel.org/docs/Introduction/faq.html) before redistributing Open Babel or its code. This inventory does not relicense it under MIT.
- **Python:** separately installed runtime. Its [license and third-party notices](https://docs.python.org/3/license.html) apply to the selected interpreter distribution. The application uses the standard-library SQLite interface; SQLite's [public-domain statement](https://www.sqlite.org/copyright.html) does not replace Python's own terms.
- **UCSF ChimeraX:** used externally by the dataset preparer for inspection/box placement, not imported or bundled by FragEvoLab. It has a [non-commercial license](https://github.com/RBVI/ChimeraX/blob/develop/License.md), with separate commercial licensing and [citation guidance](https://www.cgl.ucsf.edu/chimerax/docs/credits.html). Its use is acknowledged as preparation provenance, not as a runtime requirement.

## Redistribution scope

License snapshots document the inspected versions; they are not a replacement for the complete notices shipped with a future wheel, executable, container or installer. In particular, NumPy/SciPy binary wheels can contain BLAS/LAPACK/compiler-runtime components, and chemistry/image packages can include additional native libraries. A distributor bundling those binaries must review that exact distribution. No scientific paper citation grants software redistribution rights, and no project citation requirement is added to the MIT license.

No separately vendored third-party implementation was identified in the application source. This is not a proof of the provenance of every line. Calls into RDKit, Meeko and other installed packages are identified above. License text snapshots themselves retain their original terms. No external JavaScript chart library or font bundle is embedded in the HTML-report source.

## Data and assets

The [test-data provenance](docs/test_data.md) distinguishes PDB data, modified receptor files, PubChem structures and project-generated results. The PDB archive's CC0 terms and NCBI/PubChem policies are linked there; the project's MIT license is not asserted over third-party annotations or publications.

`assets/fragevolab-logo.png` is the project logo created for this project with AI assistance. It is not an RDKit, Vina, PubChem or UCSF logo and does not imply affiliation with those projects.

## Copyleft and the MIT source release

See [licensing and distribution](docs/licensing.md) for the source-release rationale and requirements for future bundled distributions. Meeko is LGPL-2.1, Gemmi is MPL-2.0, and the Open Babel fallback uses a separate GPL-2.0 command-line program. NumPy/SciPy Linux notices include **GPL-3.0-or-later WITH GCC-exception-3.1** for libgfortran and LGPL-2.1-or-later for libquadmath. These must not be reduced to a plain GPLv3 label or omitted from binary redistribution notices. The current source architecture does not require changing FragEvoLab's original code from MIT to GPLv3.

## Version-specific package records

These records identify the exact declared or observed Python versions. Metadata is supporting evidence; full license texts and component-specific notices control.

- [rdkit 2026.3.6](https://pypi.org/project/rdkit/2026.3.6/)
- [meeko 0.8.0](https://pypi.org/project/meeko/0.8.0/)
- [gemmi 0.7.5](https://pypi.org/project/gemmi/0.7.5/)
- [numpy 2.3.5](https://pypi.org/project/numpy/2.3.5/)
- [scipy 1.17.0](https://pypi.org/project/scipy/1.17.0/)
- [Pillow 12.3.0](https://pypi.org/project/Pillow/12.3.0/)
- [rich 15.0.0](https://pypi.org/project/rich/15.0.0/)
- [filelock 4.0.12](https://pypi.org/project/filelock/4.0.12/)
- [Jinja2 3.1.6](https://pypi.org/project/Jinja2/3.1.6/)
- [pytest 9.1.1](https://pypi.org/project/pytest/9.1.1/)
- [MarkupSafe 3.0.4](https://pypi.org/project/MarkupSafe/3.0.4/)
- [markdown-it-py 4.2.0](https://pypi.org/project/markdown-it-py/4.2.0/)
- [mdurl 0.1.2](https://pypi.org/project/mdurl/0.1.2/)
- [Pygments 2.21.0](https://pypi.org/project/Pygments/2.21.0/)
- [packaging 26.3](https://pypi.org/project/packaging/26.3/)
- [pluggy 1.7.0](https://pypi.org/project/pluggy/1.7.0/)
- [iniconfig 2.3.1](https://pypi.org/project/iniconfig/2.3.1/)
