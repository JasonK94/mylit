"""
GeoMx Tools - A Python package for processing GeoMx spatial transcriptomics data

Modules:
--------
- merge_count_matrix: Merge count matrices from multiple GeoMx experiments
- merge_metadata: Extract and merge metadata from SegmentProperties  
- pool_columns: Pool and harmonize metadata columns
- create_pm_dictionary: Create structured dictionary for PM annotation

Usage:
------
    from geomx_tools import merge_count_matrices, merge_metadata
    
    # Merge count matrices
    counts = merge_count_matrices('path/to/xlsx/files')
    
    # Merge metadata
    metadata = merge_metadata('path/to/xlsx/files')
"""

__version__ = '0.1.0'
__author__ = 'Jaecheon Lee'

# Import main functions
from .merge_count_matrix import merge_count_matrices
from .merge_metadata import merge_metadata

__all__ = [
    'merge_count_matrices',
    'merge_metadata',
]

