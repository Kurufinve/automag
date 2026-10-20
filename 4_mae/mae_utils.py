"""
Utility functions for MAE workflow scripts.

This module provides common functionality shared across MAE workflow scripts
to eliminate code redundancy and improve maintainability.

Contains:
- Cell type determination from input parameters
- Structure file loading with cell type awareness
- Parameter extraction (U, J, kpts, encut)
- Configuration file handling
- Path construction utilities
"""

import glob
import os
from pathlib import Path
from typing import Tuple, Optional, Dict, Any
import numpy as np
from pymatgen.core.structure import Structure


def get_cell_type_from_params(standardize_cell: bool = True,
                              use_primitive_cell: bool = True) -> Tuple[str, str]:
    """
    Determine cell type identifiers from standardization parameters.
    
    Args:
        standardize_cell: Whether standardization is applied
        use_primitive_cell: Whether to use primitive cell (only if standardize_cell=True)
    
    Returns:
        Tuple of (expected_cell_type, cell_type):
            - expected_cell_type: Used for structure file search ('original', 'primitive', 'conventional')
            - cell_type: Used for MAE directory naming ('input_cell', 'primitive_cell', 'conventional_cell')
    
    Examples:
        >>> get_cell_type_from_params(False, False)
        ('original', 'input_cell')
        >>> get_cell_type_from_params(True, True)
        ('primitive', 'primitive_cell')
        >>> get_cell_type_from_params(True, False)
        ('conventional', 'conventional_cell')
    """
    if not standardize_cell:
        expected_cell_type = 'original'
        cell_type = 'input_cell'
    elif use_primitive_cell:
        expected_cell_type = 'primitive'
        cell_type = 'primitive_cell'
    else:
        expected_cell_type = 'conventional'
        cell_type = 'conventional_cell'
    
    return expected_cell_type, cell_type


def load_processed_structure(configuration: str,
                             standardize_cell: bool = True,
                             use_primitive_cell: bool = True,
                             fallback_to_legacy: bool = True,
                             verbose: bool = True) -> Tuple[Optional[str], Optional[Structure], str, str]:
    """
    Load processed structure file with cell type awareness.
    
    This function implements the standard pattern used across all MAE scripts
    for finding and loading the correct processed structure file.
    
    Args:
        configuration: Magnetic configuration name (e.g., 'fm1', 'afm1')
        standardize_cell: Whether standardization was applied
        use_primitive_cell: Whether primitive cell was used
        fallback_to_legacy: If True, try legacy naming patterns if primary search fails
        verbose: If True, print status messages
    
    Returns:
        Tuple of (structure_path, structure, expected_cell_type, cell_type):
            - structure_path: Path to loaded structure file (None if not found)
            - structure: Loaded pymatgen Structure object (None if not found)
            - expected_cell_type: Cell type identifier for file search
            - cell_type: Cell type identifier for directory naming
    
    Examples:
        >>> path, struct, exp_type, folder_type = load_processed_structure('fm1')
        Found processed structure: setting001_fm1_primitive.vasp
          Cell type: primitive
        >>> path, struct, exp_type, folder_type = load_processed_structure('afm1', standardize_cell=False)
        Found processed structure: setting001_afm1_original.vasp
          Cell type: original
    """
    # Determine cell types
    expected_cell_type, cell_type = get_cell_type_from_params(standardize_cell, use_primitive_cell)
    
    # Search for processed file with the specific cell type
    # Pattern: setting*_{configuration}_{cell_type}.vasp
    processed_files = glob.glob(f'setting*_{configuration}_{expected_cell_type}.vasp')
    
    if processed_files:
        structure_path = processed_files[0]
        try:
            structure = Structure.from_file(structure_path)
            if verbose:
                print(f"Found processed structure: {structure_path}")
                print(f"  Cell type: {expected_cell_type}")
            return structure_path, structure, expected_cell_type, cell_type
        except Exception as e:
            if verbose:
                print(f"Warning: Could not load structure from {structure_path}: {e}")
            return None, None, expected_cell_type, cell_type
    
    if fallback_to_legacy:
        # Fallback: try legacy naming patterns
        if verbose:
            print(f"Warning: No processed structure file found matching pattern: setting*_{configuration}_{expected_cell_type}.vasp")
            print(f"Trying legacy naming patterns...")
        
        # Try old "processed" naming (without cell type suffix)
        legacy_processed = glob.glob(f'*{configuration}_processed.vasp')
        legacy_standardized = glob.glob(f'*{configuration}_standardized.vasp')
        
        if legacy_processed:
            structure_path = legacy_processed[0]
            try:
                structure = Structure.from_file(structure_path)
                if verbose:
                    print(f"Using legacy processed file: {structure_path}")
                return structure_path, structure, expected_cell_type, cell_type
            except Exception as e:
                if verbose:
                    print(f"Warning: Could not load structure from {structure_path}: {e}")
                return None, None, expected_cell_type, cell_type
        
        elif legacy_standardized:
            structure_path = legacy_standardized[0]
            try:
                structure = Structure.from_file(structure_path)
                if verbose:
                    print(f"Using legacy standardized file: {structure_path}")
                return structure_path, structure, expected_cell_type, cell_type
            except Exception as e:
                if verbose:
                    print(f"Warning: Could not load structure from {structure_path}: {e}")
                return None, None, expected_cell_type, cell_type
    
    # Not found
    if verbose:
        print("ERROR: Structure file not found!")
        print(f"Expected pattern: setting*_{configuration}_{expected_cell_type}.vasp")
        if fallback_to_legacy:
            print(f"  Or legacy: *{configuration}_processed.vasp")
            print(f"  Or legacy: *{configuration}_standardized.vasp")
    
    return None, None, expected_cell_type, cell_type


def extract_hubbard_uj_from_params(params: Dict[str, Any]) -> Tuple[float, float]:
    """
    Extract Hubbard U and J values from VASP parameters.
    
    Finds the U and J values for the first magnetic atom (ldaul > 0).
    
    Args:
        params: VASP parameters dictionary containing 'ldauu', 'ldauj', 'ldaul'
    
    Returns:
        Tuple of (U, J) values
    
    Examples:
        >>> params = {'ldauu': [5.2, 0, 0], 'ldauj': [0.9, 0, 0], 'ldaul': [2, -1, -1]}
        >>> extract_hubbard_uj_from_params(params)
        (5.2, 0.9)
        >>> params = {'ldauu': [0.0], 'ldauj': [0.0], 'ldaul': []}
        >>> extract_hubbard_uj_from_params(params)
        (0.0, 0.0)
    """
    ldauu_val = params.get('ldauu', [0.0])
    ldauj_val = params.get('ldauj', [0.0])
    ldaul_val = params.get('ldaul', [])
    
    # Find the first magnetic atom (ldaul > 0)
    try:
        magnetic_index = next(i for i, x in enumerate(ldaul_val) if x > 0)
        U = ldauu_val[magnetic_index]
        J = ldauj_val[magnetic_index]
    except (StopIteration, IndexError):
        # No magnetic atoms or empty lists
        U = 0.0
        J = 0.0
    
    return U, J


def extract_convergence_params_from_params(params: Dict[str, Any],
                                          use_first: bool = True) -> Tuple[int, int]:
    """
    Extract k-points and ENCUT from VASP parameters.
    
    Handles both single values and lists. If list, returns first value by default.
    
    Args:
        params: VASP parameters dictionary containing 'kpts', 'encut'
        use_first: If True and params are lists, use first value; if False, raise error
    
    Returns:
        Tuple of (kpts_val, encut_val)
    
    Examples:
        >>> params = {'kpts': 20, 'encut': 830}
        >>> extract_convergence_params_from_params(params)
        (20, 830)
        >>> params = {'kpts': [20, 25, 30], 'encut': [830, 900]}
        >>> extract_convergence_params_from_params(params, use_first=True)
        (20, 830)
    """
    kpts = params.get('kpts')
    encut = params.get('encut')
    
    # Handle kpts
    if isinstance(kpts, list):
        if use_first:
            kpts_val = kpts[0]
        else:
            raise ValueError("kpts is a list; expected single value")
    else:
        kpts_val = kpts
    
    # Handle encut
    if isinstance(encut, list):
        if use_first:
            encut_val = encut[0]
        else:
            raise ValueError("encut is a list; expected single value")
    else:
        encut_val = encut
    
    return int(kpts_val), int(encut_val)


def construct_mae_directory_path(path_to_automag: str,
                                 formula: str,
                                 reduced_formula: str,
                                 calculator: str,
                                 configuration: str,
                                 U: float,
                                 J: float,
                                 kpts_val: int,
                                 encut_val: int,
                                 cell_type: str,
                                 n_atoms: int,
                                 struct_suffix: str = '') -> Path:
    """
    Construct standardized MAE directory path.
    
    Creates path following the pattern:
    {AUTOMAG_PATH}/CalcFold/{reduced_formula}{struct_suffix}/{formula}_{cell_type}/{calculator}/{configuration}/mae_U{U}_J{J}_K{kpts}_EN{encut}_{cell_type}_{n_atoms}atoms
    
    Args:
        path_to_automag: Base AUTOMAG path
        formula: Chemical formula (e.g., 'Fe12O18')
        reduced_formula: Reduced chemical formula (e.g., 'Fe2O3')
        calculator: Calculator type (e.g., 'vasp')
        configuration: Magnetic configuration (e.g., 'fm1')
        U: Hubbard U value
        J: Hubbard J value
        kpts_val: K-points value
        encut_val: ENCUT value
        cell_type: Cell type folder name (e.g., 'primitive_cell', 'input_cell', 'conventional_cell', 'supercell')
        n_atoms: Number of atoms
        struct_suffix: Optional structure suffix
    
    Returns:
        Path object for MAE directory
    
    Examples:
        >>> construct_mae_directory_path('/home/user/automag', 'Fe12O18', 'Fe2O3', 'vasp', 'fm1', 
        ...                              5.2, 0.0, 20, 830, 'primitive_cell', 10)
        PosixPath('/home/user/automag/CalcFold/Fe2O3/Fe12O18_primitive_cell/vasp/fm1/mae_U5.2_J0.0_K20_EN830_primitive_cell_10atoms')
    """
    calcfold_path = Path(path_to_automag) / 'CalcFold'
    mae_base_dir = calcfold_path / f"{reduced_formula}{struct_suffix}" / f"{formula}_{cell_type}" / calculator / configuration
    mae_dir = mae_base_dir / f"mae_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type}_{n_atoms}atoms"
    
    return mae_dir


def construct_config_filename(configuration: str,
                              U: float,
                              J: float,
                              kpts_val: int,
                              encut_val: int,
                              cell_type: str,
                              n_atoms: int) -> str:
    """
    Construct standardized MAE configuration filename.
    
    Creates filename following the pattern:
    {configuration}_mae_config_U{U}_J{J}_K{kpts}_EN{encut}_{cell_type}_{n_atoms}atoms.txt
    
    Args:
        configuration: Magnetic configuration (e.g., 'fm1')
        U: Hubbard U value
        J: Hubbard J value
        kpts_val: K-points value
        encut_val: ENCUT value
        cell_type: Cell type folder name (e.g., 'primitive_cell')
        n_atoms: Number of atoms
    
    Returns:
        Configuration filename string
    
    Examples:
        >>> construct_config_filename('fm1', 5.2, 0.0, 20, 830, 'primitive_cell', 10)
        'fm1_mae_config_U5.2_J0.0_K20_EN830_primitive_cell_10atoms.txt'
    """
    return f'mae_config_{configuration}_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type}_{n_atoms}atoms.txt'


def find_config_file(configuration: str,
                    base_dir: Optional[Path] = None,
                    exact_params: Optional[Dict[str, Any]] = None,
                    verbose: bool = True) -> Optional[str]:
    """
    Find MAE configuration file using pattern matching.
    
    Can search for exact match using parameters or use wildcard pattern.
    
    Args:
        configuration: Magnetic configuration name
        base_dir: Directory to search in (defaults to current directory)
        exact_params: If provided, construct exact filename from params dict with keys:
                     'U', 'J', 'kpts_val', 'encut_val', 'cell_type', 'n_atoms'
        verbose: If True, print status messages
    
    Returns:
        Config filename if found, None otherwise
    
    Examples:
        >>> find_config_file('fm1', exact_params={'U': 5.2, 'J': 0.0, 'kpts_val': 20,
        ...                                        'encut_val': 830, 'cell_type': 'primitive_cell',
        ...                                        'n_atoms': 10})
        'fm1_mae_config_U5.2_J0.0_K20_EN830_primitive_cell_10atoms.txt'
    """
    if base_dir is None:
        base_dir = Path.cwd()
    else:
        base_dir = Path(base_dir)
    
    # Try exact match if parameters provided
    if exact_params:
        config_file = construct_config_filename(
            configuration=configuration,
            U=exact_params['U'],
            J=exact_params['J'],
            kpts_val=exact_params['kpts_val'],
            encut_val=exact_params['encut_val'],
            cell_type=exact_params['cell_type'],
            n_atoms=exact_params['n_atoms']
        )
        
        config_path = base_dir / config_file
        if config_path.exists():
            if verbose:
                print(f"Found config file: {config_file}")
            return str(config_path)
    
    # Try wildcard pattern
    pattern = f'mae_config_{configuration}*.txt'
    matching_configs = list(base_dir.glob(pattern))
    
    if matching_configs:
        config_file = str(matching_configs[0])
        if verbose:
            print(f"Using config file: {Path(config_file).name}")
        return config_file
    
    if verbose:
        print(f"Warning: No config file found matching pattern: {pattern}")
    
    return None


def read_mae_config_parameters(config_file: str,
                               verbose: bool = True) -> Dict[str, Any]:
    """
    Read MAE configuration file and extract parameters.
    
    Parses configuration file and returns dictionary of parameters.
    
    Args:
        config_file: Path to configuration file
        verbose: If True, print status messages
    
    Returns:
        Dictionary with extracted parameters (may contain: 'ncl_magmoms', 'cell_type',
        'processed_formula', 'ldauu_val', 'ldauj_val', 'kpts_val', 'encut_val',
        'easy_axis', 'hard_axis')
    
    Examples:
        >>> params = read_mae_config_parameters('fm1_mae_config_U5.2_J0.0_K20_EN830_primitive_cell_10atoms.txt')
        Reading configuration from: fm1_mae_config_U5.2_J0.0_K20_EN830_primitive_cell_10atoms.txt
          → Loaded 10 NCL magnetic moments
          → Cell type: primitive_cell
          → Processed formula: Fe4O6
    """
    params = {}
    
    if verbose:
        print(f"\nReading configuration from: {config_file}")
    
    try:
        with open(config_file, 'r') as f:
            for line in f:
                if 'NCL magmoms:' in line:
                    # Parse the magmoms - format is a list of tuples
                    magmoms_str = line.split('NCL magmoms:')[1].strip()
                    params['ncl_magmoms'] = eval(magmoms_str)
                    if verbose:
                        print(f"  → Loaded {len(params['ncl_magmoms'])} NCL magnetic moments")
                
                elif 'Cell type:' in line:
                    params['cell_type'] = line.split('Cell type:')[1].strip()
                    if verbose:
                        print(f"  → Cell type: {params['cell_type']}")
                
                elif 'Processed formula:' in line:
                    params['processed_formula'] = line.split('Processed formula:')[1].strip()
                    if verbose:
                        print(f"  → Processed formula: {params['processed_formula']}")
                
                elif 'LDAUU:' in line:
                    ldauu_str = line.split('LDAUU:')[1].strip()
                    params['ldauu_val'] = eval(ldauu_str)
                
                elif 'LDAUJ:' in line:
                    ldauj_str = line.split('LDAUJ:')[1].strip()
                    params['ldauj_val'] = eval(ldauj_str)
                
                elif 'K-points:' in line:
                    params['kpts_val'] = int(line.split('K-points:')[1].strip())
                
                elif 'ENCUT:' in line:
                    params['encut_val'] = int(line.split('ENCUT:')[1].strip())
                
                elif line.startswith('Easy axis:'):
                    # Parse easy axis: format is "Easy axis: [x, y, z]"
                    axis_str = line.split('Easy axis:')[1].strip()
                    axis_str = axis_str.strip('[]')
                    params['easy_axis'] = np.array([float(x.strip()) for x in axis_str.split(',')])
                    if verbose:
                        print(f"  → Loaded easy axis: {params['easy_axis']}")
                
                elif line.startswith('Hard axis:'):
                    # Parse hard axis: format is "Hard axis: [x, y, z]"
                    axis_str = line.split('Hard axis:')[1].strip()
                    axis_str = axis_str.strip('[]')
                    params['hard_axis'] = np.array([float(x.strip()) for x in axis_str.split(',')])
                    if verbose:
                        print(f"  → Loaded hard axis: {params['hard_axis']}")
    
    except Exception as e:
        if verbose:
            print(f"Warning: Could not parse some parameters from config: {e}")
    
    return params


def load_collinear_configuration(path_to_automag: str,
                                 formula: str,
                                 configuration: str,
                                 struct_suffix: str = '',
                                 calculator: str = 'vasp',
                                 verbose: bool = True) -> Tuple[str, list, int]:
    """
    Load a magnetic configuration and its final moments from collinear results.

    Resolves the collinear results directory, finds the setting file that
    contains ``configuration``, and reads the final magnetic moments from the
    collinear single-point output.

    Returns:
        Tuple of (structure_path, final_magmoms, setting_number)
    """
    path_to_coll = Path(path_to_automag) / '2_coll' / f'{formula}{struct_suffix}' / f'{formula}_{calculator}'
    if not path_to_coll.exists():
        path_to_coll = Path(path_to_automag) / '2_coll' / f'{formula}_{calculator}'
        if not path_to_coll.exists():
            path_to_coll = Path(path_to_automag) / '2_coll'

    if verbose:
        print(f'Loading configuration from: {path_to_coll}')

    path_to_trials = path_to_coll / f'trials_{formula}'

    setting = 1
    magmom = None
    while (path_to_trials / f'configurations{setting:03d}.txt').exists():
        with open(path_to_trials / f'configurations{setting:03d}.txt', 'rt') as f:
            for line in f:
                values = line.split()
                if values and values[0] == configuration:
                    magmom = [int(item) for item in values[1:]]
                    break
            if magmom is not None:
                break
        setting += 1

    if magmom is None:
        raise IOError(f'Configuration {configuration} not found in {path_to_trials}')

    structure_path = path_to_trials / f'setting{setting:03d}.vasp'

    # Read the final magnetic moments from the collinear single-point results.
    # The CalcFold file is named after the formula passed in by the caller (the
    # same convention as the rest of automag), not a re-derived one.
    calc_results_file = Path(path_to_automag) / 'CalcFold' / \
        f'{formula}{struct_suffix}_singlepoint_{calculator}.txt'

    if verbose:
        print(f'Reading final magmoms from: {calc_results_file}')

    final_magmoms = None
    with open(calc_results_file, 'rt') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        values = line.split()
        if values and values[0] == configuration:
            magmoms_str = lines[i + 1]
            full_magmom_str = magmoms_str.split('final_magmoms=')[1].strip().strip('[]')
            final_magmoms = [float(num) for num in full_magmom_str.split()]
            break

    if final_magmoms is None:
        raise IOError(f'Final magmoms for {configuration} not found in {calc_results_file}')

    if verbose:
        print(f'Found configuration {configuration} in setting{setting:03d}')

    return str(structure_path), final_magmoms, setting


def _map_magmoms_to_transformed_structure(original: Structure,
                                          transformed: Structure,
                                          original_magmoms: list) -> list:
    """
    Map magnetic moments from the original structure onto a transformed cell.

    Matches by species first, then by minimum-image distance, so that
    primitive/conventional conversions keep the correct moment on each site.
    """
    transformed_magmoms = []

    for i, trans_site in enumerate(transformed.sites):
        target_species = trans_site.species
        target_frac = trans_site.frac_coords

        best_match_idx = None
        min_distance = float('inf')

        for j, orig_site in enumerate(original.sites):
            if orig_site.species != target_species:
                continue

            orig_frac = orig_site.frac_coords
            min_image_dist = float('inf')
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        image_frac = orig_frac + np.array([dx, dy, dz])
                        image_cart = original.lattice.get_cartesian_coords(image_frac)
                        target_cart = original.lattice.get_cartesian_coords(target_frac)
                        min_image_dist = min(
                            min_image_dist, np.linalg.norm(image_cart - target_cart)
                        )

            if min_image_dist < min_distance:
                min_distance = min_image_dist
                best_match_idx = j

        if best_match_idx is None:
            frac_dists = np.linalg.norm(original.frac_coords - target_frac, axis=1)
            best_match_idx = int(np.argmin(frac_dists))
            print(f"  WARNING: No species match for site {i} ({target_species}); "
                  f"using nearest neighbour (site {best_match_idx})")

        transformed_magmoms.append(original_magmoms[best_match_idx])

    return transformed_magmoms


def process_structure(structure_path: str,
                      magmoms: list,
                      output_name: str,
                      symmetrize: bool = True,
                      to_primitive: bool = True,
                      symprec: float = 0.1,
                      verbose: bool = True) -> Tuple[str, list, Structure]:
    """
    Standardize (optionally to a primitive cell) and emit NCL magnetic moments.

    Cell selection follows the shared ``standardize_cell`` / ``use_primitive_cell``
    inputs, and magnetic moments are remapped onto the transformed cell.

    Returns:
        Tuple of (processed_structure_path, ncl_magmoms, processed_structure)
    """
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    pmg_structure = Structure.from_file(structure_path)
    original_natoms = len(pmg_structure)
    pmg_structure.add_site_property("magmom", magmoms)

    if not symmetrize:
        cell_type = "original"
    elif to_primitive:
        cell_type = "primitive"
    else:
        cell_type = "conventional"

    def _cell_filename(name: str, cell: str) -> str:
        if name.endswith('_processed.vasp'):
            return f"{name.replace('_processed.vasp', '')}_{cell}.vasp"
        parts = name.rsplit('.', 1)
        if len(parts) == 2:
            return f"{parts[0]}_{cell}.{parts[1]}"
        return f"{name}_{cell}"

    if verbose:
        print(f"\nStructure Processing:")
        print(f"  Original structure: {original_natoms} atoms")
        print(f"  Symmetrize: {symmetrize}")
        print(f"  Primitive cell: {to_primitive if symmetrize else 'N/A (no standardization)'}")
        print(f"  Cell type: {cell_type}")

    if not symmetrize:
        processed_structure = pmg_structure
    else:
        try:
            sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec)
            if to_primitive:
                processed_structure = sga.get_primitive_standard_structure()
            else:
                processed_structure = sga.get_conventional_standard_structure()

            mapped_magmoms = _map_magmoms_to_transformed_structure(
                pmg_structure, processed_structure, magmoms
            )
            if "magmom" in processed_structure.site_properties:
                processed_structure.remove_site_property("magmom")
            processed_structure.add_site_property("magmom", mapped_magmoms)

            if verbose:
                print(f"  → {original_natoms} → {len(processed_structure)} atoms "
                      f"(space group {sga.get_space_group_symbol()})")
        except Exception as e:
            if verbose:
                print(f"  WARNING: Symmetrization failed: {e}")
                print(f"  → Falling back to original structure")
            processed_structure = pmg_structure
            cell_type = "original"

    dynamic_output_name = _cell_filename(output_name, cell_type)
    processed_structure.to(filename=dynamic_output_name, fmt='POSCAR')

    processed_magmom = list(processed_structure.site_properties['magmom'])
    ncl_magmom = [(0, 0, m) for m in processed_magmom]

    if verbose:
        print(f"  → Saved to: {dynamic_output_name}")
        print(f"  → NCL magmoms: {len(ncl_magmom)} values")

    return dynamic_output_name, ncl_magmom, processed_structure
