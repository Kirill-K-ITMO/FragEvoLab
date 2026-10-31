# FragEvoLab 3.1.4 chemical rules

The `src/evolution/chemistry_rules.py` module is used for atom-level mutations and BRICS fragment assembly/replacement. Operators first select the mutation level and operation using the configured probabilities, then select an eligible site and modification. If no eligible modification exists, the attempt is rejected with a recorded reason. Disabled operators are not enabled automatically.

## Atom-level operations

The policy permits neutral C(4), N(3), O(2), F/Cl/Br/I(1), N+(4), and O-(1). This is a conservative set of supported states, not a complete table of elemental valences.

- Adding C/N/O/F replaces one hydrogen. Both implicit hydrogens and bracket-specified hydrogen counts are considered.
- Substitution of a neutral atom with C/N/O requires a compatible bond-order sum; the hydrogen count is recalculated.
- Only eligible terminal neutral atoms can be removed; hydrogens are restored on the neighboring atom.
- Single/double bond changes require sufficient hydrogen capacity at both ends.
- Atom-level operations do not edit aromatic regions. Isotopically labeled atoms, explicit neighboring hydrogen atoms, radicals, and unsupported states are protected.

Every operation requires RDKit sanitization, a connectivity check, and the absence of dummy atoms and radicals. The input molecule is not modified. Existing coordinates are removed; new geometry is generated during ligand preparation.

## BRICS

The module uses compatible RDKit BRICS attachment types. It checks dummy-atom degree and compatibility between the dummy bond order and the new bond order. Bonds to the replaced dummy atoms are removed before the real bond is added. Dummy atoms are allowed only during intermediate assembly; the final molecule must have no unconsumed attachment points.

Labels are assumed to originate from actual BRICS fragmentation. The module does not validate manually fabricated BRICS labels.

## Stereochemistry and scope

The configuration of an edited tetrahedral center is cleared: the algorithm does not infer a new configuration. Stereocenters remote from the edit are preserved. Operations affecting specified E/Z geometry are rejected. Enhanced stereo groups are not supported by this operator set.

Valid valence does not establish molecular stability or synthetic accessibility. Salts are not automatically stripped. Supporting additional states requires dedicated rules and tests.

## Statistics

Mutation events record the operator, site, original and resulting structures, rejection reason, and counts of excluded modifications by reason. Excluded sites are not equivalent to evolutionary attempts. A fragment reward is calculated from the fitness change in an evaluated single replacement. Multistep chains receive a combined reward. Preparation and docking failures are recorded separately without assigning a negative fragment reward.

## Affinity estimates

Chemical rules do not use Ki. Standard reports convert docking scores into estimated Ki/pKi at 298.15 K; these are not experimental constants. See the [conversion methodology](docs/affinity_estimates.md).
