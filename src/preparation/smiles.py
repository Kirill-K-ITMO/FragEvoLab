"""smiles - FragEvoLab module."""

import re
import csv
from collections import defaultdict
from itertools import product
from pathlib import Path
from typing import List, Dict, Optional

class SmilesGenerator:
    """
    Generate SMILES strings from templates with substituents.
    
    Use templates such as 'C{1}C{2}C(=O)O', where {n} marks a substituent position.
    """
    
    def __init__(self, project_root: Optional[Path] = None):
        """
        Initialize the SMILES generator.
        
        Args:
            project_root: Project root directory (default: current directory)
        """
        self.project_root = project_root or Path.cwd()
        self.data_dir = self.project_root / "data"
        self.ligands_csv = self.data_dir / "ligands.csv"
        
        # Dictionaries storing generation data
        self.template = ""
        self.radicals_dict = defaultdict(list)
        self.exclusion_rules = []
        self.generated_smiles = []
        self.combinations_info = []
    
    def set_template(self, template: str) -> bool:
        """
        Set the SMILES generation template.
        
        Args:
            template: String with {n} placeholders, for example 'C{1}C{2}C(=O)O'
            
        Returns:
            True if the template is valid, False otherwise
        """
        # Check for at least one placeholder
        pattern_ids = set(re.findall(r'\{(\d+)\}', template))
        
        if not pattern_ids:
            print("Error: the template contains no {n} substituent positions")
            return False
        
        self.template = template
        self.radicals_dict.clear()
        self.exclusion_rules.clear()
        self.generated_smiles.clear()
        self.combinations_info.clear()
        
        print(f"Template set: {template}")
        print(f"Substituent positions found: {sorted(pattern_ids, key=int)}")
        
        return True
    
    def add_radicals(self, position: str, radicals: List[str]) -> bool:
        """
        Add substituents for a specified template position.
        
        Args:
            position: Position number (for example, '1')
            radicals: List of substituents for this position (empty string = no substituent)
            
        Returns:
            True if substituents were added successfully
        """
        if not self.template:
            print("Error: set a template using set_template() first")
            return False
        
        # Check that the position exists in the template
        pattern_ids = set(re.findall(r'\{(\d+)\}', self.template))
        if position not in pattern_ids:
            print(f"Error: position {{{position}}} does not exist in the template")
            return False
        
        if not radicals:
            print(f"Error: the substituent list for position {{{position}}} is empty")
            return False
        
        self.radicals_dict[position] = radicals
        print(f"Substituents added for position {{{position}}}: {radicals}")
        
        return True
    
    def add_exclusion_rule(self, rule: Dict[str, str]) -> bool:
        """
        Add a rule excluding combinations.
        
        Args:
            rule: Dictionary {position: substituent}, for example {'1': '(C)', '2': '(CC)'}
            
        Returns:
            True if the rule was added successfully
        """
        # Check that all positions in the rule exist
        pattern_ids = set(re.findall(r'\{(\d+)\}', self.template))
        
        for position in rule.keys():
            if position not in pattern_ids:
                print(f"Error: position {position} does not exist in the template")
                return False
        
        self.exclusion_rules.append(rule)
        print(f"Exclusion rule added: {rule}")
        
        return True
    
    def generate_combinations(self) -> List[str]:
        """
        Generate all SMILES combinations from the template and substituents.
        
        Returns:
            List of generated SMILES strings
        """
        if not self.template:
            print("Error: no template has been set")
            return []
        
        if not self.radicals_dict:
            print("Error: no substituents have been added")
            return []
        
        # Create lists of position numbers and their substituents
        keys = sorted(self.radicals_dict.keys(), key=int)
        values = [self.radicals_dict[k] for k in keys]
        
        # Check that substituents are defined for all template positions
        pattern_ids = set(re.findall(r'\{(\d+)\}', self.template))
        missing_positions = pattern_ids - set(keys)
        
        if missing_positions:
            print(f"Warning: no substituents defined for positions: {sorted(missing_positions, key=int)}")
            return []
        
        # Generate the Cartesian product of all substituent lists
        combinations = []
        combinations_info = []
        
        for combo in product(*values):
            # Copy the template
            temp_str = self.template
            combo_dict = {}
            
            # Replace all placeholders with substituents
            for position, radical in zip(keys, combo):
                temp_str = temp_str.replace(f"{{{position}}}", radical)
                combo_dict[position] = radical
            
            # Remove any remaining empty parentheses
            temp_str = temp_str.replace("()", "")
            
            combinations.append(temp_str)
            combinations_info.append(combo_dict)
        
        print(f"\nCombinations generated: {len(combinations)}")
        
        # Apply exclusion rules
        filtered_combinations = []
        filtered_info = []
        
        for smiles, combo_dict in zip(combinations, combinations_info):
            exclude = False
            
            for rule in self.exclusion_rules:
                # Check whether all positions in the rule match
                match = all(
                    position in combo_dict and combo_dict[position] == radical
                    for position, radical in rule.items()
                )
                
                if match:
                    exclude = True
                    break
            
            if not exclude:
                filtered_combinations.append(smiles)
                filtered_info.append(combo_dict)
        
        excluded_count = len(combinations) - len(filtered_combinations)
        if excluded_count > 0:
            print(f"Combinations excluded by rules: {excluded_count}")
        
        print(f"Final combination count: {len(filtered_combinations)}")
        
        self.generated_smiles = filtered_combinations
        self.combinations_info = filtered_info
        
        return filtered_combinations
    
    def save_to_csv(self, 
                     project_name: str = "ligand", 
                     append: bool = False,
                     output_file: Optional[Path] = None) -> bool:
        """
        Save generated SMILES to a CSV file.
        
        Args:
            project_name: Prefix for molecule names
            append: If True, append to the existing file; otherwise overwrite it
            output_file: Output file path (default: data/ligands.csv)
            
        Returns:
            True if saved successfully
        """
        if not self.generated_smiles:
            print("Error: no generated SMILES to save")
            return False
        
        # Determine the output file
        if output_file is None:
            output_file = self.ligands_csv
        
        # Create the data directory if absent
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Remove invalid characters from the project name
        project_name = re.sub(r'[^a-zA-Z0-9_]', '_', project_name.strip())
        if not project_name:
            project_name = "ligand"
        
        from io import StringIO
        from src.preparation.inputs import read_ligand_table
        from src.core.calculation_cache import atomic_write
        fields, old_rows = (["name", "smiles"], [])
        if append and output_file.exists():
            fields, old_rows = read_ligand_table(output_file)
        existing_names = {row["name"].casefold() for row in old_rows}
        csv_data = []
        next_index = 1
        for smiles in self.generated_smiles:
            while f"{project_name}{next_index}".casefold() in existing_names:
                next_index += 1
            name = f"{project_name}{next_index}"
            existing_names.add(name.casefold())
            csv_data.append((name, smiles))
            next_index += 1
        buffer = StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=fields)
        writer.writeheader()
        writer.writerows(old_rows)
        writer.writerows({"name": name, "smiles": smiles} for name, smiles in csv_data)
        atomic_write(output_file, buffer.getvalue().encode("utf-8"))

        action = "appended to" if append else "saved to"
        print(f"\n✓ {len(csv_data)} SMILES {action} file: {output_file}")
        print(f"  Molecule names: {csv_data[0][0]} ... {csv_data[-1][0]}")
        
        return True
    
    def interactive_generation(self) -> List[str]:
        """
        Generate SMILES interactively using user prompts.
        
        Returns:
            List of generated SMILES
        """
        print("TEMPLATE-BASED SMILES LIBRARY GENERATION")
        
        # Step 1: Request a template
        while True:
            template_input = input("\nEnter a template with {n} positions (for example, 'c1ccccc1-{1}-{2}'): ").strip()
            if self.set_template(template_input):
                break
        
        # Step 2: Request substituents for each position
        pattern_ids = sorted(set(re.findall(r'\{(\d+)\}', self.template)), key=int)
        
        for position in pattern_ids:
            while True:
                print(f"\nSubstituents for position {{{position}}}:")
                radicals_input = input(
                    f"  Enter substituents separated by spaces/commas (enter 'space' for no substituent): "
                ).strip()
                
                if not radicals_input:
                    print("  Error: input cannot be empty")
                    continue
                
                # Process input
                radicals = [
                    "" if r.lower() == "space" else r
                    for r in re.split(r'[,\s]+', radicals_input)
                    if r
                ]
                
                if radicals:
                    self.add_radicals(position, radicals)
                    break
                else:
                    print("  Error: no substituents entered")
        
        # Step 3: Generate combinations
        smiles_list = self.generate_combinations()
        
        if not smiles_list:
            print("\n❌ Combination generation failed")
            return []
        
        # Show the first few combinations
        print(f"\nExamples of generated SMILES:")
        for i, smiles in enumerate(smiles_list[:10], 1):
            print(f"  {i}. {smiles}")
        
        if len(smiles_list) > 10:
            print(f"  ... and another {len(smiles_list) - 10} combinations")
        
        # Step 4: Exclude combinations (optional)
        exclude_option = input("\nExclude any combinations? (yes/no): ").strip().lower()
        
        if exclude_option in ['да', 'yes', 'y', 'д']:
            print("\nEnter combinations to exclude in the format: <position><substituent>")
            print("For example, '1(C) 2(CC)' excludes combinations with (C) at position {1} AND (CC) at position {2}")
            print("Enter 'stop' to finish")
            
            while True:
                rule_input = input("\nEnter a combination to exclude (or 'stop'): ").strip()
                
                if rule_input.lower() == 'stop':
                    break
                
                if not rule_input:
                    continue
                
                # Parse the rule
                rule_dict = {}
                valid_rule = True
                
                for part in rule_input.split():
                    match = re.match(r'^(\d+)(\S+)$', part)
                    
                    if not match:
                        print(f"  ❌ Invalid format in '{part}'. Use: <number><substituent>")
                        valid_rule = False
                        break
                    
                    position, radical = match.groups()
                    
                    if position not in pattern_ids:
                        print(f"  ❌ Error: position {position} does not exist in the template")
                        valid_rule = False
                        break
                    
                    rule_dict[position] = radical
                
                if valid_rule and rule_dict:
                    self.add_exclusion_rule(rule_dict)
            
            # Regenerate with exclusion rules
            if self.exclusion_rules:
                print("\nRegenerate with exclusion rules...")
                smiles_list = self.generate_combinations()
        
        # Step 5: Save
        save_option = input("\nSave results to data/ligands.csv? (yes/no): ").strip().lower()
        
        if save_option in ['да', 'yes', 'y', 'д']:
            project_name = input("Enter the project name (prefix for molecule names): ").strip()
            if not project_name:
                project_name = "ligand"
            
            append_option = input("Append to the existing file? (yes/no, default no): ").strip().lower()
            append = append_option in ['да', 'yes', 'y', 'д']
            
            self.save_to_csv(project_name=project_name, append=append)
        
        return smiles_list

def generate_from_parameters(template: str, 
                             radicals: Dict[str, List[str]],
                             exclusions: Optional[List[Dict[str, str]]] = None,
                             project_name: str = "ligand",
                             save: bool = True,
                             project_root: Optional[Path] = None) -> List[str]:
    """
    Generate SMILES from parameters without interactive prompts.
    Suitable for programmatic use.
    
    Args:
        template: SMILES template with {n} positions
        radicals: Dictionary {position: [list of substituents]}
        exclusions: List of exclusion rules
        project_name: Prefix for molecule names
        save: Whether to save results to CSV
        project_root: Project root directory
        
    Returns:
        List of generated SMILES
        
    """
    generator = SmilesGenerator(project_root=project_root)
    
    # Set the template
    if not generator.set_template(template):
        return []
    
    # Add substituents
    for position, radical_list in radicals.items():
        if not generator.add_radicals(position, radical_list):
            return []
    
    # Add exclusion rules
    if exclusions:
        for rule in exclusions:
            generator.add_exclusion_rule(rule)
    
    # Generation
    smiles_list = generator.generate_combinations()
    
    # Save
    if save and smiles_list:
        generator.save_to_csv(project_name=project_name, append=True)
    
    return smiles_list

