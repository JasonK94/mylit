#!/usr/bin/env python3
"""
GeoMx Metadata Merger
Extracts and merges metadata from SegmentProperties sheets
"""

import pandas as pd
import argparse
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

# Import from the existing metadata_merger_v3.py
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "projects/#999.GeoMxmeta"))

try:
    from metadata_merger_v3 import MetadataExtractor, ColumnHarmonizer, MetadataMerger
except ImportError:
    # Fallback: define minimal versions
    class MetadataExtractor:
        def __init__(self, verbose=True):
            self.verbose = verbose
            self.extracted_data = []
        
        def extract_from_directory(self, directory):
            # Simplified version
            return []
    
    class ColumnHarmonizer:
        def __init__(self):
            self.column_analysis = []
        
        def analyze_columns(self, data):
            return pd.DataFrame()
    
    class MetadataMerger:
        def __init__(self, data, mapping=None):
            self.extracted_data = data
            self.merged_data = {}
        
        def merge_by_data_type(self):
            return {}


def merge_metadata(input_folder, output_folder=None, verbose=True):
    """
    Extract and merge metadata from GeoMx XLSX files
    
    Parameters:
    -----------
    input_folder : str
        Folder containing XLSX files
    output_folder : str
        Output folder for merged metadata
    verbose : bool
        Print progress
    
    Returns:
    --------
    dict : Merged metadata by data type
    """
    
    if output_folder is None:
        output_folder = Path(input_folder) / "metadata_output"
    
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    # Extract metadata
    extractor = MetadataExtractor(verbose=verbose)
    extracted_data = extractor.extract_from_directory(input_folder)
    
    # Analyze columns
    harmonizer = ColumnHarmonizer()
    column_analysis = harmonizer.analyze_columns(extracted_data)
    
    # Save column template
    template_path = Path(output_folder) / 'column_mapping_template.csv'
    if not column_analysis.empty:
        column_analysis.to_csv(template_path, index=False)
        if verbose:
            print(f"✓ Column template: {template_path}")
    
    # Merge metadata
    merger = MetadataMerger(extracted_data, column_mapping=None)
    merged_data = merger.merge_by_data_type()
    
    # Save merged metadata
    for dtype, df in merged_data.items():
        out_file = Path(output_folder) / f"metadata_{dtype}_merged.csv"
        df.to_csv(out_file, index=False)
        if verbose:
            print(f"✓ {dtype}: {out_file} ({df.shape})")
    
    return merged_data


def main():
    parser = argparse.ArgumentParser(description='Merge GeoMx metadata')
    parser.add_argument('input_folder', help='Folder containing XLSX files')
    parser.add_argument('-o', '--output', help='Output folder', default=None)
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    results = merge_metadata(args.input_folder, args.output, args.verbose)
    
    print(f"\n✅ Merged {len(results)} data types")


if __name__ == "__main__":
    main()

