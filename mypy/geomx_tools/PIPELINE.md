# GeoMx Processing Pipeline

## Quick Start

```bash
# Set your data folder
DATA_FOLDER="/path/to/your/geomx/xlsx/files"
OUTPUT_FOLDER="./geomx_output"

# Run complete pipeline
bash run_pipeline.sh $DATA_FOLDER $OUTPUT_FOLDER
```

## Manual Step-by-Step

### 1. Merge Count Matrices

```bash
python -m geomx_tools.merge_count_matrix \
    $DATA_FOLDER \
    -o $OUTPUT_FOLDER/01_counts \
    -v
```

Or in Python:
```python
from geomx_tools import merge_count_matrices

counts = merge_count_matrices(
    'path/to/xlsx',
    output_file='output/merged.csv',
    verbose=True
)
```

### 2. Merge Metadata

```bash
python -m geomx_tools.merge_metadata \
    $DATA_FOLDER \
    -o $OUTPUT_FOLDER/02_metadata \
    -v
```

Or in Python:
```python
from geomx_tools import merge_metadata

metadata = merge_metadata(
    'path/to/xlsx',
    output_folder='output/metadata',
    verbose=True
)
```

### 3. Filter Duplicates

```python
# In Python/Jupyter
from geomx_tools.filters import remove_duplicate_aois

filtered = remove_duplicate_aois(
    'output/01_counts/merged_WTA.csv',
    output_path='output/03_filtered/merged_WTA_filtered.csv'
)
```

### 4. Match Metadata to Matrix

```python
from geomx_tools.matchers import match_metadata_to_matrix

matched = match_metadata_to_matrix(
    metadata_csv='output/02_metadata/metadata_WTA_merged.csv',
    matrix_csv='output/03_filtered/merged_WTA_filtered.csv',
    output_dir='output/04_matched'
)
```

### 5. Create PM Dictionary

```python
from geomx_tools.pm_dict import create_structured_dictionary

pm_dict = create_structured_dictionary(
    metadata_folder='output/02_metadata',
    output_folder='output/05_pm_dictionary'
)
```

### 6. PM Annotation (Manual)

1. Match PM names to dictionary
2. Distribute to PMs
3. Collect annotations
4. Merge back

### 7. Final Merge

```python
from geomx_tools.pm_dict import apply_pm_annotations

final_metadata = apply_pm_annotations(
    pm_dictionary='output/05_pm_dictionary/STRUCTURED_PM_DICTIONARY.csv',
    annotated_files='output/06_annotated/*.csv',
    output_file='output/07_final/final_metadata_annotated.csv'
)
```

## Pipeline Flow

```
Input XLSX Files
        ↓
[1] Merge Count Matrices → merged_*.csv
        ↓
[2] Merge Metadata → metadata_*_merged.csv
        ↓
[3] Filter Duplicates → merged_*_filtered.csv
        ↓
[4] Match Metadata → metadata_*_MATCHED.csv
        ↓
[5] Create PM Dictionary → STRUCTURED_PM_DICTIONARY.csv
        ↓
[6] PM Annotation (MANUAL)
        ↓
[7] Apply Annotations → final_metadata_annotated.csv
```

## Output Structure

```
output/
├── 01_counts/
│   ├── merged_WTA.csv
│   ├── merged_CTA.csv
│   └── merged_MOUSE.csv
│
├── 02_metadata/
│   ├── metadata_WTA_merged.csv
│   ├── metadata_CTA_merged.csv
│   └── metadata_MOUSE_merged.csv
│
├── 03_filtered/
│   ├── merged_WTA_filtered.csv
│   └── ...
│
├── 04_matched/
│   ├── metadata_WTA_MATCHED.csv
│   └── ...
│
├── 05_pm_dictionary/
│   ├── STRUCTURED_PM_DICTIONARY.csv
│   ├── ID_MAPPING_*.csv
│   └── ...
│
├── 06_annotated/
│   └── (PM annotated files)
│
└── 07_final/
    └── final_metadata_annotated.csv
```

