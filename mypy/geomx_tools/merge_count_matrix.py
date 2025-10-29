#!/usr/bin/env python3
"""
GeoMx Count Matrix Merger
Merges count matrices from multiple GeoMx XLSX files
"""

import pandas as pd
import numpy as np
import os
import argparse
from pathlib import Path
from functools import reduce
import warnings
warnings.filterwarnings('ignore')


def detect_data_type(df):
    """Detect data type from genome build"""
    if 'GenomeBuild' in df.columns:
        builds = df['GenomeBuild'].astype(str).str.lower()
        if builds.str.contains('grcm', na=False).any():
            return 'MOUSE'
        if builds.str.contains('grch', na=False).any():
            if 'CodeClass' in df.columns and 'TargetName' in df.columns:
                neg = df[df['CodeClass'] == 'Negative']
                if not neg.empty and 'Negative Probe' in neg['TargetName'].unique():
                    return 'CTA'
            return 'WTA'
    return 'UNKNOWN'


def merge_count_matrices(input_folder, output_file=None, verbose=True):
    """
    Merge count matrices from BioProbeCountMatrix sheets
    
    Parameters:
    -----------
    input_folder : str
        Folder containing XLSX files
    output_file : str
        Output CSV file path (optional)
    verbose : bool
        Print progress
    
    Returns:
    --------
    dict : Dictionary of merged dataframes by type
    """
    
    BASE_COLS = ['ProbeName', 'ProbeDisplayName', 'TargetName', 'HUGOSymbol', 
                 'Accessions', 'GenomeBuild', 'GenomicPosition', 'AnalyteType', 
                 'CodeClass', 'ProbePool', 'TargetGroup', 'GeneID']
    
    # Find all XLSX files
    xlsx_files = list(Path(input_folder).glob("*.xlsx"))
    
    if verbose:
        print(f"Found {len(xlsx_files)} XLSX files")
    
    # Group by data type
    data_by_type = {'WTA': [], 'CTA': [], 'MOUSE': []}
    
    for xlsx_file in xlsx_files:
        try:
            df = pd.read_excel(xlsx_file, sheet_name="BioProbeCountMatrix", engine='openpyxl')
            
            data_type = detect_data_type(df)
            if data_type == 'UNKNOWN':
                continue
            
            # Rename AOI columns
            file_prefix = xlsx_file.stem
            aoi_cols = [c for c in df.columns if c not in BASE_COLS]
            rename_dict = {col: f"{file_prefix}_{col}" for col in aoi_cols}
            df.rename(columns=rename_dict, inplace=True)
            
            data_by_type[data_type].append(df)
            
            if verbose:
                print(f"  ✓ {xlsx_file.name} ({data_type}): {len(aoi_cols)} AOIs")
                
        except Exception as e:
            if verbose:
                print(f"  ✗ {xlsx_file.name}: {e}")
    
    # Merge by type
    results = {}
    for dtype, dfs in data_by_type.items():
        if not dfs:
            continue
            
        merged = reduce(lambda l, r: pd.merge(l, r, on='ProbeName', how='outer'), dfs)
        results[dtype] = merged
        
        if verbose:
            print(f"\n{dtype}: {merged.shape}")
        
        # Save if output specified
        if output_file:
            out_path = Path(output_file).parent / f"merged_{dtype}.csv"
            merged.to_csv(out_path, index=False)
            if verbose:
                print(f"  Saved to {out_path}")
    
    return results


def main():
    parser = argparse.ArgumentParser(description='Merge GeoMx count matrices')
    parser.add_argument('input_folder', help='Folder containing XLSX files')
    parser.add_argument('-o', '--output', help='Output folder', default='./output')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Create output folder
    Path(args.output).mkdir(parents=True, exist_ok=True)
    
    # Run merger
    output_file = Path(args.output) / "merged.csv"
    results = merge_count_matrices(args.input_folder, str(output_file), args.verbose)
    
    print(f"\n✅ Complete! Merged {len(results)} data types")
    print(f"   Output: {args.output}")


if __name__ == "__main__":
    main()

