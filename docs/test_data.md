# JAK example: data provenance and limitations

The [bundled example](../data/README.md) exercises preparation, docking, evolutionary generation, caching and checkpoints. It is not a validated activity, selectivity or redocking benchmark.

## Included files

- [ligands.csv](../data/ligands.csv): 40 unique input molecules with `name,smiles` columns.
- [ligand_sources.csv](../data/ligand_sources.csv): PubChem CIDs, record URLs, identity verification dates and unspecified tetrahedral-center counts.
- [receptor_sources.json](../data/receptor_sources.json): receptor sources, preparation metadata and hashes of the included PDBQT/box files.
- Prepared receptors and box configurations are in `data/proteins/` and `data/configs/`. The matching evolution settings are in the repository root.

## Receptors

| Role | Receptor | Source | Chain | Box center (angstrom) | Box size (angstrom) |
|---|---|---|---|---|---|
| Target | JAK1 | [3EYG](https://www.rcsb.org/structure/3EYG) | A | 10.96, 14.23, -14.78 | 20, 20, 20 |
| Off-target | JAK2 | [6BBV](https://www.rcsb.org/structure/6BBV) | A | -5.49, 1.36, 8.30 | 20, 20, 20 |

Roles are configuration choices, not experimental classifications. Cite the structure entries and primary papers in [references.md](references.md).

Receptors were prepared using Meeko, with ChimeraX used for inspection and box placement, on Windows 10 in a Miniforge environment. JAK1 contains 2,785 atom records and JAK2 contains 2,928, including polar hydrogens. Both contain chain A, with no free ligand or water molecules. Crystal ligands MI1 (3EYG) and D7D (6BBV) supplied spatial references.

Incomplete residues omitted using Meeko's `--delete_bad_res_from_box_radius 4` option:

- JAK1: GLU A946, ASN A950 and GLU A1029.
- JAK2: LYS A1053.

JAK1 retains phosphorylated PTR A1034 and PTR A1035. Omission of incomplete residues does not by itself establish suitability for quantitative docking.

Exact Meeko/ChimeraX versions, original download dates, preparation pH and a complete preparation recipe are unavailable. The JAK1 PRO A912 repair procedure and alternate-location choice at A1058 are not fully documented. Missing-residue reconstruction, receptor energy minimization, protonation-state validation and crystal-pose redocking are not established. The included files are ready for technical execution, but these preparation limitations must be considered in scientific use.

## Ligands

The structures originate from the PubChem Compound records linked in `ligand_sources.csv`. CHEMBL, BDBM and SCHEMBL names are retained input labels; the structural sources here are the linked PubChem records. Verification dates are not acquisition dates. No assay data or experimental selectivity classification is supplied.

Twelve inputs contain unspecified tetrahedral stereocenters. The supplied SMILES are retained without stereoisomer enumeration or assignment of missing configurations. Unspecified stereochemistry is not evidence of a racemate. One prepared 3D structure does not represent an ensemble of stereoisomers. Per-compound counts describe tetrahedral centers, not every possible type of stereoisomerism.

The preparation pipeline adds hydrogens, embeds with ETKDGv3, optimizes using UFF with bounded retries, and produces PDBQT through Meeko. The example uses seed 42. There is no pH-based microstate selection, salt stripping or tautomer enumeration. The 650 Da limit is an example configuration choice, not a universal docking limit. Estimated Ki/pKi conversion uses 298.15 K; this is neither a Vina simulation temperature nor a preparation parameter.

## Attribution and data terms

PDB archive data are provided under CC0 1.0 according to [RCSB policies](https://www.rcsb.org/pages/policies). Cite depositors and associated publications. Prepared PDBQT files are modified derivatives of the deposits.

Retain individual PubChem URLs/CIDs and cite PubChem as described in [references.md](references.md). [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/) distinguish government-created information from third-party contributions. The project's MIT license does not relicense external data. These tables contain structure identifiers and SMILES, not copied assay tables or article text. Experimental activity claims require their own assay sources and conditions.
