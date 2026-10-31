# Docking-score-derived Ki and pKi

The standard FragEvoLab workflow uses **T = 298.15 K (25 °C)**. This temperature is used for the score conversion only; it is not a Vina simulation-temperature setting.

The existing numerical conversion is preserved:

```text
score_J_per_mol = score_kcal_per_mol * 4184
Ki_M = exp(score_J_per_mol / (8.314463 * T))
pKi = -log10(Ki_M)
```

The numerical Ki is expressed in mol/L relative to the conventional standard concentration of 1 mol/L; pKi is dimensionless. The thermodynamic analogy uses delta_G_standard = R*T*ln(Kd/c_standard). Substitution of a docking score does not make the result a measured dissociation or inhibition constant. The historical field name `ki` is retained for compatibility.

Docking scores approximate aspects of binding and must not be treated as exact Gibbs free energies. These estimates are not experimental Ki values and do not establish biological activity or selectivity. Evolutionary selection uses the configured docking-score fitness, not these derived values.

The default constructor temperature remains 298.15 K. Explicit custom temperatures remain supported by the Python API. Reports generated from imported values cannot infer their original conversion temperature. Their annotation describes the standard workflow and explicitly asks users to verify imported or custom-temperature results.

HTML docking and database reports carry the annotation. CSV docking reports retain their numerical columns and receive an adjacent `.csv.notes.txt` file. The English build uses `name` for the molecule-name column; the legacy Russian input header remains supported. Retain this note when distributing the CSV. Previously generated files are not rewritten. Evolution reports contain docking energies and fitness without a Ki conversion.


Numerical limits: temperature must be finite and greater than zero; the default
remains 298.15 K. pKi is computed directly from the score in logarithmic form.
If Ki exceeds floating-point range, the API returns `None` and the report labels
the numerical limitation. Zero Ki is never converted to pKi=0.
