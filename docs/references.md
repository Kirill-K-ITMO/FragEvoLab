# Scientific references and software citations

Cite the methods actually used in a reported calculation, together with the software versions, configuration and data sources. A software citation is not a substitute for a license notice; see [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).

## Methods used by the application

| Component | Use in this code | Reference |
|---|---|---|
| RDKit | SMILES parsing, sanitization, descriptors, molecular drawings and fingerprints | [RDKit citation guidance](https://www.rdkit.org/docs/Overview.html#citing-the-rdkit). Cite "RDKit: Open-source cheminformatics", https://www.rdkit.org, with the actual version. The pinned version is 2026.3.6; the [Zenodo collection](https://doi.org/10.5281/zenodo.591637) is a version-discovery link, not a verified DOI for this particular wheel. |
| BRICS | Fragmentation and attachment-type compatibility in `src/evolution/fragments.py`, `operators.py` and `chemistry_rules.py` | Degen J, Wegscheid-Gerlach C, Zaliani A, Rarey M. On the Art of Compiling and Using "Drug-Like" Chemical Fragment Spaces. *ChemMedChem*. 2008;3:1503-1507. [doi:10.1002/cmdc.200800178](https://doi.org/10.1002/cmdc.200800178). See also the [RDKit implementation tutorial](https://greglandrum.github.io/rdkit-blog/posts/2025-08-15-BRICS-tutorial.html). |
| ETKDG | Foundation of the conformer generation method used in ligand preparation | Riniker S, Landrum GA. Better Informed Distance Geometry: Using What We Know To Improve Conformation Generation. *J Chem Inf Model*. 2015;55:2562-2574. [doi:10.1021/acs.jcim.5b00654](https://doi.org/10.1021/acs.jcim.5b00654). |
| ETKDGv3 | `AllChem.ETKDGv3()` in `src/preparation/ligands.py` | Wang S, Witek J, Landrum GA, Riniker S. Improving Conformer Generation for Small Rings and Macrocycles Based on Distance Geometry and Experimental Torsional-Angle Preferences. *J Chem Inf Model*. 2020;60:2044-2058. [doi:10.1021/acs.jcim.0c00025](https://doi.org/10.1021/acs.jcim.0c00025). |
| UFF | Geometry optimization after embedding; bounded convergence retries | Rappe AK, Casewit CJ, Colwell KS, Goddard WA III, Skiff WM. UFF, a full periodic table force field for molecular mechanics and molecular dynamics simulations. *J Am Chem Soc*. 1992;114:10024-10035. [doi:10.1021/ja00051a040](https://doi.org/10.1021/ja00051a040). |
| Meeko | Ligand parameterization and PDBQT serialization | Santos-Martins D et al. Meeko: Molecule Parametrization and Software Interoperability for Docking and Beyond. *J Chem Inf Model*. 2025;65:13045-13050. [doi:10.1021/acs.jcim.5c02271](https://doi.org/10.1021/acs.jcim.5c02271). [Source and documentation](https://github.com/forlilab/Meeko). |
| AutoDock Vina | External docking executable | Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J Comput Chem*. 2010;31:455-461. [doi:10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334). Also cite Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J Chem Inf Model*. 2021;61:3891-3898. [doi:10.1021/acs.jcim.1c00203](https://doi.org/10.1021/acs.jcim.1c00203). |
| Morgan fingerprints | Radius 2, 2048 bits; Tanimoto filtering and diversity summaries in selection, engine and reports | Rogers D, Hahn M. Extended-Connectivity Fingerprints. *J Chem Inf Model*. 2010;50:742-754. [doi:10.1021/ci100050t](https://doi.org/10.1021/ci100050t). Report the actual [RDKit fingerprint settings](https://www.rdkit.org/docs/GettingStartedInPython.html), not just "ECFP"; the code uses RDKit Morgan fingerprints and default chirality handling. |
| Crippen logP | Molecular filtering and descriptor reports | Wildman SA, Crippen GM. Prediction of Physicochemical Parameters by Atomic Contributions. *J Chem Inf Comput Sci*. 1999;39:868-873. [doi:10.1021/ci990307l](https://doi.org/10.1021/ci990307l). |
| TPSA | Descriptor reports via RDKit | Ertl P, Rohde B, Selzer P. Fast Calculation of Molecular Polar Surface Area as a Sum of Fragment-Based Contributions and Its Application to the Prediction of Drug Transport Properties. *J Med Chem*. 2000. [doi:10.1021/jm000942e](https://doi.org/10.1021/jm000942e). |
| SA score | Optional docking-report descriptor loaded from RDKit `Contrib/SA_Score`; not an evolutionary fitness term | Ertl P, Schuffenhauer A. Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. *J Cheminform*. 2009;1:8. [doi:10.1186/1758-2946-1-8](https://doi.org/10.1186/1758-2946-1-8). The RDKit contributed implementation includes modifications to the published method; see its source and [notice](../licenses/rdkit/SA_Score-NOTICE.txt). |
| Open Babel | Optional external `obabel` fallback | O'Boyle NM et al. Open Babel: An open chemical toolbox. *J Cheminform*. 2011;3:33. [doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33). Cite the version only when this backend was actually used. |

The application uses RDKit APIs and its installed BRICS rules. It does not claim to reproduce every implementation detail of the original BRICS publication. Valid valence and compatible attachment types do not demonstrate synthetic accessibility.

Adaptive fragment/operator weighting, bounded mutation policies, target/off-target aggregation, caching and checkpoint handling are project-specific implementation choices. Their source is in `src/evolution/` and `src/core/`; no external algorithm paper is asserted as their exact specification. This release does not use Pareto-front optimization.

The [Ki/pKi conversion note](affinity_estimates.md) documents 298.15 K, units and limitations. Docking scores are not experimental binding free energies, and converted values are not measured inhibition constants.

## Supporting numerical libraries

- Harris CR et al. Array programming with NumPy. *Nature*. 2020;585:357-362. [doi:10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2).
- Virtanen P et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*. 2020;17:261-272. [doi:10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2).

These support the installed chemistry stack; they are not separate docking engines. Gemmi, Pillow, Rich, filelock, Jinja2 and the testing dependencies are identified with versions and licenses in the third-party inventory. A separate scientific paper is not invented for each infrastructure package.

## Test data and preparation tools

- **JAK1, PDB 3EYG:** [structure record](https://www.rcsb.org/structure/3EYG), [dataset DOI](https://doi.org/10.2210/pdb3EYG/pdb). Williams NK et al. Structural study of JAK-specific inhibitor complexes. *J Mol Biol*. 2009;387:219-232. [doi:10.1016/j.jmb.2009.01.041](https://doi.org/10.1016/j.jmb.2009.01.041).
- **JAK2, PDB 6BBV:** [structure record](https://www.rcsb.org/structure/6BBV), [dataset DOI](https://doi.org/10.2210/pdb6BBV/pdb). Vazquez ML et al. Study identifying PF-04965842 as a selective JAK1 clinical candidate. *J Med Chem*. 2018;61:1130-1152. [doi:10.1021/acs.jmedchem.7b01598](https://doi.org/10.1021/acs.jmedchem.7b01598). The paper's JAK1 focus does not change the identity of the JAK2 structure 6BBV.
- **Ligand structures:** Kim S et al. PubChem 2025 update. *Nucleic Acids Res*. 2025;53:D1516-D1525. [doi:10.1093/nar/gkae1059](https://doi.org/10.1093/nar/gkae1059). Individual CIDs and record URLs are in [ligand_sources.csv](../data/ligand_sources.csv). See [PubChem citation guidance](https://pubchem.ncbi.nlm.nih.gov/docs/citation-guidelines). Experimental activity claims would require their own assay/publication sources; none are supplied by this technical dataset.
- **External receptor inspection:** Pettersen EF et al. UCSF ChimeraX: Structure visualization for researchers, educators, and developers. *Protein Sci*. 2021;30:70-82. [doi:10.1002/pro.3943](https://doi.org/10.1002/pro.3943). ChimeraX was used by the dataset preparer for receptor inspection and box placement, not as an application dependency. [Official citation guidance](https://www.cgl.ucsf.edu/chimerax/docs/credits.html).

Details, confirmed preparation steps and unresolved provenance are in [test_data.md](test_data.md). 

## Citing FragEvoLab

Use the root [CITATION.cff](../CITATION.cff) to cite FragEvoLab. Kirill V. Karhanin is the software author. Alexander S. Novikov is acknowledged in [AUTHORS.md](../AUTHORS.md) for scientific supervision, separately from software authorship.
