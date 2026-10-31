# Licensing and distribution

FragEvoLab's original source code is released under [MIT](../LICENSE). Third-party software and datasets retain their own terms. The current source archive contains no dependency wheels, native libraries, Vina/Open Babel executables, copied Meeko implementation, or ChimeraX software. See [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md) for versions and retained notices.

For this source-only architecture, the dependencies identified below do not require changing the original FragEvoLab source license to GPLv3. This is a project-specific reading of the licenses, not a guarantee about every future installer or legal dispute. License obligations depend on what is copied, linked, modified and distributed, not merely on whether installation is performed by pip.

## Meeko: LGPL, not GPLv3

The installed Meeko 0.8.0 distribution identifies its license as LGPL-2.1 and contains the GNU Lesser General Public License version 2.1. Its broad PyPI classifier mentions LGPLv2+; this project retains the exact supplied license rather than inferring a more specific SPDX expression from the classifier.

FragEvoLab imports the separately installed Meeko Python package through its public API. No Meeko source is incorporated into FragEvoLab. LGPL-2.1 sections 5 and 6 permit a work using the library to have different terms when the applicable requirements are met. The MIT source release permits user modification and does not restrict debugging or reverse engineering. The normal Python import mechanism allows users to install an interface-compatible modified Meeko; the requirements pin identifies a supported version and is not a prohibition on replacement.

Retain the prominent Meeko notice and [LGPL text](../licenses/meeko/LICENSE). If distributing Meeko itself, provide its corresponding source and notices by a method permitted by the license. A frozen executable, statically linked package or bundled environment needs separate consideration of replacement/relinking and source-availability requirements. Do not treat LGPL as permission to remove library notices or relicense Meeko under MIT.

Sources: [Meeko 0.8.0](https://pypi.org/project/meeko/0.8.0/), [upstream license](https://github.com/forlilab/Meeko/blob/v0.8.0/LICENSE), [LGPL-2.1 sections 5-6](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html).

## Open Babel: separate GPL-2.0 command-line program

The optional fallback invokes `obabel` through `subprocess` and exchanges ordinary molecular files. FragEvoLab does not import `openbabel`/`pybel`, link `libopenbabel`, or bundle the executable. Open Babel's own FAQ explicitly distinguishes library linking from using its command-line interface with independent software. On this basis, the current fallback does not require relicensing FragEvoLab under GPL.

The referenced Open Babel 3.1.1 distribution supplies GPL version 2. Do not silently relabel it as GPLv3 or assume that all of its components permit that upgrade. Direct library integration or copying implementation code would require a new licensing assessment. Switching FragEvoLab to GPLv3 is not a universal solution to GPL-2.0-only compatibility issues.

If distributing Open Babel binaries, comply with GPLv2 source-distribution requirements for those binaries. A generic homepage link is not a substitute for an applicable source-distribution method. Users installing it separately obtain it from its distributor.

Sources: [Open Babel FAQ](https://openbabel.org/docs/Introduction/faq.html), [library-use guidance](https://openbabel.org/docs/UseTheLibrary/intro.html), [3.1.1 COPYING](https://github.com/openbabel/openbabel/blob/openbabel-3-1-1/COPYING), [GNU discussion of separate programs](https://www.gnu.org/licenses/old-licenses/gpl-2.0-faq.html#MereAggregation).

## GPLv3 text in NumPy/SciPy notices

The retained Linux wheel notices for NumPy and SciPy include `libgfortran` under `GPL-3.0-or-later WITH GCC-exception-3.1`. The exception is an integral part of the permission: it permits covered runtime libraries to be combined with independent modules under other terms when the compilation conditions are met. Thus the appearance of the GPLv3 text in these files does not, by itself, place FragEvoLab under GPLv3.

Those wheel notices also describe `libquadmath` under LGPL-2.1-or-later. Native-library inventories differ by wheel, operating system and build. NumPy/SciPy's own BSD licenses do not replace bundled-library terms. FragEvoLab distributes neither these wheels nor their native libraries. If a future release bundles them, retain all notices and satisfy source obligations for the actual libraries; the GCC exception does not remove obligations when distributing the runtime library itself.

Sources: [NumPy notice](../licenses/numpy/LICENSE.txt), [SciPy notice](../licenses/scipy/LICENSE.txt), [GCC Runtime Library Exception 3.1](https://www.gnu.org/licenses/gcc-exception-3.1.en.html), [GCC exception FAQ](https://www.gnu.org/licenses/gcc-exception-3.1-faq.html.en).

## Other components

- **Gemmi (MPL-2.0):** copyleft applies at file level. Separate original FragEvoLab files may remain MIT. If distributing Gemmi or modified Gemmi files, satisfy its notice and source-availability requirements. See the [Mozilla FAQ](https://www.mozilla.org/en-US/MPL/2.0/FAQ/).
- **RDKit, Jinja2 and other permissively licensed dependencies:** retain required copyright/license notices if redistributing their code or binaries. RDKit's core BSD notice is distinct from the wheel-packaging MIT notice. Its optional SA-score component has a separate Novartis BSD notice, retained in `licenses/rdkit/`.
- **AutoDock Vina (Apache-2.0):** an external executable installed separately. Apache-2.0 does not require relicensing the original FragEvoLab source. If bundling Vina, preserve the applicable license, attribution/NOTICE material and notices for that exact binary's dependencies. See [Vina's license](https://vina.scripps.edu/LICENSE/).
- **ChimeraX:** used externally for receptor inspection/box placement, not a runtime dependency. Its non-commercial license is not MIT, GPL or a general permission to redistribute ChimeraX. Commercial use requires appropriate UCSF licensing. No ChimeraX software or documentation is bundled. See [UCSF's license text](https://github.com/RBVI/ChimeraX/blob/develop/License.md).

## Source release versus bundled application

For this source release, keep `LICENSE`, `THIRD_PARTY_NOTICES.md`, `licenses/` and the data provenance files. Install dependencies separately and allow compatible library replacement. Scientific references acknowledge methods; they do not replace license compliance.

Before distributing an executable, Docker image, virtual environment or dependency wheelhouse, identify the exact third-party artifacts and their licenses, include required notices and corresponding-source access, and check copyleft replacement/relinking requirements. Preserve upstream rights; do not label the complete bundle as exclusively MIT.

PDB/PubChem inputs have their own [data provenance and terms](test_data.md). The MIT software license does not grant rights over third-party annotations, publications or trademarks. A license review cannot establish the originality of every line solely from its text; contributors must disclose copied or adapted third-party material and retain its notices.
