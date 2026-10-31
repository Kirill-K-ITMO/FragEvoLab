"""affinity - FragEvoLab module."""

import math
from typing import Dict, List, Tuple, Optional

class AffinityCalculator:
    """
    Estimate Ki/pKi from docking scores; the default temperature is 298.15 K.
    A docking score is not a measured free energy or biological activity.
    
    Calculates:
    - Estimated Ki (mol/L)
    - pKi = -log10(Ki)
    - Formatted Ki values
    """
    
    def __init__(self, temperature: float = 298.15):
        """
        Initialize the calculator.
        
        Args:
            temperature: Temperature in kelvin (default: 298.15 K = 25°C)
        """
        self._check_temperature(temperature)
        self.temperature = temperature
        self.R = 8.314463  # Universal gas constant, J/(mol·K)
        self.kcal_to_j = 4184  # Conversion factor from kcal to J
    
    def calculate_ki(self, delta_g: float, temperature: Optional[float] = None) -> Optional[float]:
        """
        Estimate Ki from a docking score (kcal/mol) using a thermodynamic analogy.
        Use the instance temperature by default: 298.15 K.
        This is not an experimental Ki; the standard-state concentration is 1 M.
        
        Formula: Ki = exp(ΔG / RT)
        
        Args:
            delta_g: Docking score in kcal/mol
            temperature: Temperature in K (use self.temperature if None)
            
        Returns:
            Estimated Ki in mol/L (M), or None for floating-point underflow/overflow.
        """
        if temperature is None:
            temperature = self.temperature
        
        exponent = self._exponent(delta_g, temperature)
        try:
            ki = math.exp(exponent)
        except OverflowError:
            return None
        return ki if ki > 0 and math.isfinite(ki) else None

    @staticmethod
    def _check_temperature(temperature):
        if not math.isfinite(temperature) or temperature <= 0:
            raise ValueError("Temperature must be finite and greater than zero kelvin")

    def _exponent(self, delta_g, temperature):
        self._check_temperature(temperature)
        if not math.isfinite(delta_g):
            raise ValueError("Docking score must be finite")
        value = delta_g * (self.kcal_to_j / (self.R * temperature))
        if not math.isfinite(value):
            raise ValueError("Score/temperature exceeds the supported numerical range")
        return value

    def calculate_pki_from_score(self, delta_g, temperature=None):
        return -self._exponent(delta_g, self.temperature if temperature is None else temperature) / math.log(10)

    def calculate_pki(self, ki: float) -> float:
        """
        Calculate pKi from Ki.
        
        Formula: pKi = -log10(Ki)
        
        Args:
            ki: Estimated inhibition constant in mol/L (M)
            
        Returns:
            pKi (dimensionless)
        """
        if ki is None or not math.isfinite(ki) or ki <= 0:
            raise ValueError("Ki must be finite and positive; use calculate_pki_from_score for extreme scores")
        
        return -math.log10(ki)
    
    def format_ki(self, ki: float) -> str:
        """
        Format Ki as a human-readable string with units.
        
        Args:
            ki: Estimated inhibition constant in mol/L (M)
            
        Returns:
            Formatted string (for example: "5.23 nM")
        """
        if ki is None:
            return "Outside floating-point range"
        if ki >= 1:
            return f"{ki:.2e} M"
        elif 1e-3 <= ki < 1:
            return f"{ki * 1e3:.2f} mM"
        elif 1e-6 <= ki < 1e-3:
            return f"{ki * 1e6:.2f} µM"
        elif 1e-9 <= ki < 1e-6:
            return f"{ki * 1e9:.2f} nM"
        elif 1e-12 <= ki < 1e-9:
            return f"{ki * 1e12:.2f} pM"
        else:
            return f"{ki:.2e} M"
    
    def process_docking_results(self, docking_results: Dict) -> Dict:
        """
        Add estimated Ki and pKi to docking results.
        
        Args:
            docking_results: Dictionary {ligand: {protein: {affinity: ...}}}
            
        Returns:
            Enriched dictionary with Ki, pKi, and ki_formatted
        """
        enriched_results = {}
        
        for ligand_name, proteins_data in docking_results.items():
            enriched_results[ligand_name] = {}
            
            for protein_name, data in proteins_data.items():
                # Copy the original data
                result = data.copy()
                
                # Retrieve ΔG
                delta_g = data.get('affinity', data.get('delta_g'))
                
                if delta_g is not None:
                    # Calculate Ki and pKi
                    ki = self.calculate_ki(delta_g)
                    pki = self.calculate_pki_from_score(delta_g)
                    
                    # Add to the result
                    result['delta_g'] = delta_g
                    result['ki'] = ki
                    result['pki'] = pki
                    result['ki_formatted'] = self.format_ki(ki)
                else:
                    result['delta_g'] = None
                    result['ki'] = None
                    result['pki'] = None
                    result['ki_formatted'] = 'N/A'
                
                enriched_results[ligand_name][protein_name] = result
        
        return enriched_results
    
    def calculate_statistics(self, enriched_results: Dict, protein: str) -> Dict:
        """
        Calculate statistics for a receptor.
        
        Args:
            enriched_results: Results enriched with Ki/pKi
            protein: Receptor name
            
        Returns:
            Statistics dictionary (min, max, mean for ΔG, Ki, pKi)
        """
        delta_g_values = []
        ki_values = []
        pki_values = []
        
        for ligand_data in enriched_results.values():
            if protein in ligand_data:
                data = ligand_data[protein]
                
                if data.get('delta_g') is not None:
                    delta_g_values.append(data['delta_g'])
                
                if data.get('ki') is not None:
                    ki_values.append(data['ki'])
                
                if data.get('pki') is not None:
                    pki_values.append(data['pki'])
        
        stats = {
            'delta_g': None,
            'ki': None,
            'pki': None
        }
        
        if delta_g_values:
            stats['delta_g'] = {
                'min': min(delta_g_values),
                'max': max(delta_g_values),
                'mean': sum(delta_g_values) / len(delta_g_values)
            }
        
        if ki_values:
            stats['ki'] = {
                'min': min(ki_values),
                'max': max(ki_values),
                'min_formatted': self.format_ki(min(ki_values)),
                'max_formatted': self.format_ki(max(ki_values))
            }
        
        if pki_values:
            stats['pki'] = {
                'min': min(pki_values),
                'max': max(pki_values),
                'mean': sum(pki_values) / len(pki_values)
            }
        
        return stats
    
    def get_best_ligands_by_pki(self, 
                                enriched_results: Dict, 
                                protein: str, 
                                n: int = 10) -> List[Tuple[str, float, float]]:
        """
        Return the top N ligands for a receptor ranked by pKi.
        
        Args:
            enriched_results: Enriched results
            protein: Receptor name
            n: Number of ligands
            
        Returns:
            List of (ligand_name, delta_g, pKi) tuples sorted by descending pKi
        """
        ligands_data = []
        
        for ligand_name, proteins_data in enriched_results.items():
            if protein in proteins_data:
                data = proteins_data[protein]
                
                if data.get('pki') is not None:
                    ligands_data.append((
                        ligand_name,
                        data['delta_g'],
                        data['pki']
                    ))
        
        # Sort by pKi (higher is better)
        ligands_data.sort(key=lambda x: x[2], reverse=True)
        
        return ligands_data[:n]
