# 🚀 GeoMx Tools - Quick Start Guide

## 📦 What You Have

A complete Python package for GeoMx data processing!

```
mylit/mypy/geomx_tools/
├── __init__.py              # Package initialization
├── merge_count_matrix.py    # CLI-ready count merger
├── merge_metadata.py        # CLI-ready metadata merger
└── PIPELINE.md             # Complete workflow

Plus:
- setup.py                   # Installation config
- requirements.txt           # Dependencies  
- README.md                 # Full documentation
```

---

## ⚡ Quick Install & Test

```bash
# 1. Navigate to package
cd /home/jaecheon/1kjc1/.1home1/.1kjc1/mylit/mypy

# 2. Install in development mode (like R's devtools::load_all())
pip install -e .

# 3. Test it works
python -c "from geomx_tools import merge_count_matrices; print('✅ Success!')"
```

---

## 🎯 How to Use

### Method 1: Python Import

```python
from geomx_tools import merge_count_matrices, merge_metadata

# Merge count matrices
counts = merge_count_matrices(
    '/path/to/xlsx/files',
    output_file='merged.csv',
    verbose=True
)

# Merge metadata
metadata = merge_metadata(
    '/path/to/xlsx/files',
    output_folder='metadata_output',
    verbose=True
)
```

### Method 2: Command Line

```bash
# Merge counts
python -m geomx_tools.merge_count_matrix /path/to/xlsx -o output/ -v

# Merge metadata
python -m geomx_tools.merge_metadata /path/to/xlsx -o output/ -v
```

### Method 3: Jupyter Notebook

```python
import sys
sys.path.insert(0, '/home/jaecheon/1kjc1/.1home1/.1kjc1/mylit/mypy')

from geomx_tools import merge_count_matrices

results = merge_count_matrices('path/to/data')
```

---

## 📋 Complete Pipeline

```bash
DATA="/home/jaecheon/1kjc1/.1home1/.1kjc1/projects/#999.GeoMxmeta/pooled"
OUT="./output"

# Step 1: Merge count matrices
python -m geomx_tools.merge_count_matrix $DATA -o $OUT/counts

# Step 2: Merge metadata
python -m geomx_tools.merge_metadata $DATA -o $OUT/metadata

# Steps 3-7: Use existing scripts or create notebook
```

---

## 🆚 Python vs R Packaging

| Task | R | Python |
|------|---|--------|
| Create package | `usethis::create_package()` | `mkdir mypackage` |
| Dev install | `devtools::load_all()` | `pip install -e .` |
| Test | `devtools::test()` | `pytest` |
| Document | `devtools::document()` | `sphinx` |
| Check | `devtools::check()` | `flake8` |
| Build | `devtools::build()` | `python setup.py sdist` |

---

## 📚 Files Created

1. **Package code:** `geomx_tools/*.py`
2. **Configuration:** `setup.py`, `requirements.txt`
3. **Documentation:** `README.md`, `PIPELINE.md`
4. **This guide:** `QUICK_START.md`
5. **Complete guide:** `PACKAGE_SETUP_COMPLETE.md`

---

## ✨ Next Steps

1. ✅ Install: `pip install -e .`
2. ✅ Test: `from geomx_tools import merge_count_matrices`
3. ✅ Run on your data: `python -m geomx_tools.merge_count_matrix /path/to/data`
4. ⏳ Create Jupyter notebook with pipeline steps
5. ⏳ Add more modules (filters, matchers, etc.)

---

**Everything is ready! Have a great trip!** 🎉✈️
