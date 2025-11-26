# GeoMx Tools - Python Package

A Python package for processing GeoMx spatial transcriptomics data.

## 📦 Package Structure

```
mylit/mypy/
├── geomx_tools/              # Main package directory
│   ├── __init__.py          # Package initialization
│   ├── merge_count_matrix.py   # Count matrix merger
│   ├── merge_metadata.py       # Metadata merger
│   └── ...                     # Other modules
│
├── setup.py                 # Package installation configuration
├── README.md                # This file
├── requirements.txt         # Package dependencies
│
└── tests/                   # Unit tests (optional)
    ├── __init__.py
    └── test_merge.py
```

## 🚀 Installation

### Development Mode (Editable Install)

This is similar to R's `devtools::load_all()`:

```bash
# Navigate to package directory
cd /home/jaecheon/1kjc1/.1home1/.1kjc1/mylit/mypy

# Install in editable/development mode
pip install -e .

# Now you can import from anywhere!
python -c "from geomx_tools import merge_count_matrices; print('✓ Success')"
```

### Standard Installation

```bash
pip install .
```

### From GitHub (future)

```bash
pip install git+https://github.com/yourusername/geomx_tools.git
```

## 📖 Usage

### As Python Package

```python
from geomx_tools import merge_count_matrices, merge_metadata

# Merge count matrices
counts = merge_count_matrices('path/to/xlsx/files', output_file='merged.csv')

# Merge metadata
metadata = merge_metadata('path/to/xlsx/files', output_folder='metadata_output')
```

## 🔧 Package Functions

### Core Functions

| Function Name | Input | Output | Function |
|---|---|---|---|
| `merge_count_matrices` | input_folder (path to xlsx files), output_file (optional), verbose (bool) | DataFrame with merged counts | Merges count matrices from multiple GeoMx XLSX files into a single matrix |
| `merge_metadata` | input_folder (path to xlsx files), output_folder (optional), verbose (bool) | Merged metadata DataFrame | Extracts and merges metadata from SegmentProperties across multiple files |
| `detect_data_type` | df (DataFrame) | String ('probe' or 'target') | Detects whether data contains probe-based or target-based counts |

### As Command Line Tools

After installation, you get CLI commands:

```bash
# Merge count matrices
geomx-merge-counts /path/to/xlsx/files -o output/ -v

# Merge metadata
geomx-merge-metadata /path/to/xlsx/files -o output/ -v
```

### In Jupyter Notebook

```python
import sys
sys.path.insert(0, '/home/jaecheon/1kjc1/.1home1/.1kjc1/mylit/mypy')

from geomx_tools import merge_count_matrices

# Use the functions
results = merge_count_matrices('path/to/data')
```

## 🏗️ Python Package Development Standards

### 1. Package Structure

#### Minimal Package
```
mypackage/
├── mypackage/
│   ├── __init__.py
│   └── module.py
└── setup.py
```

#### Standard Package
```
mypackage/
├── mypackage/           # Source code
│   ├── __init__.py
│   ├── core.py
│   └── utils.py
├── tests/               # Unit tests
│   ├── __init__.py
│   └── test_core.py
├── docs/                # Documentation
├── setup.py            # Installation config
├── README.md           # Package info
├── requirements.txt    # Dependencies
├── LICENSE             # License
└── .gitignore          # Git ignore
```

#### Full Package (Production)
```
mypackage/
├── src/
│   └── mypackage/
│       ├── __init__.py
│       └── ...
├── tests/
├── docs/
├── examples/
├── setup.py
├── setup.cfg
├── pyproject.toml      # PEP 518 build system
├── MANIFEST.in         # Include/exclude files
├── README.md
├── LICENSE
├── CHANGELOG.md
├── requirements.txt
├── requirements-dev.txt
└── .github/
    └── workflows/      # CI/CD
```

### 2. Key Files Explained

#### `__init__.py`
```python
"""Package initialization and public API"""

__version__ = '0.1.0'
__author__ = 'Your Name'

# Import main functions
from .module import function1, function2

# Define public API
__all__ = ['function1', 'function2']
```

#### `setup.py`
```python
from setuptools import setup, find_packages

setup(
    name="mypackage",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'pandas>=1.3.0',
        'numpy>=1.20.0',
    ],
    entry_points={
        'console_scripts': [
            'mycli=mypackage.cli:main',
        ],
    },
)
```

#### `pyproject.toml` (Modern Python)
```toml
[build-system]
requires = ["setuptools>=45", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "mypackage"
version = "0.1.0"
dependencies = [
    "pandas>=1.3.0",
    "numpy>=1.20.0",
]

[project.scripts]
mycli = "mypackage.cli:main"
```

### 3. Development Workflow

#### Step 1: Create Package Structure
```bash
mkdir -p mypackage/mypackage
touch mypackage/mypackage/__init__.py
touch mypackage/setup.py
```

#### Step 2: Install in Development Mode
```bash
cd mypackage
pip install -e .
```

This is equivalent to R's `devtools::load_all()`:
- Changes to code are immediately available
- No need to reinstall after edits
- Can test interactively

#### Step 3: Make Changes & Test
```bash
# Edit code
vim mypackage/module.py

# Test immediately (no reinstall needed!)
python -c "from mypackage import function; function()"
```

#### Step 4: Add Tests
```bash
# Install pytest
pip install pytest

# Create tests
mkdir tests
touch tests/test_module.py

# Run tests
pytest
```

#### Step 5: Build & Distribute
```bash
# Build distribution
python setup.py sdist bdist_wheel

# Upload to PyPI
pip install twine
twine upload dist/*
```

### 4. Best Practices

#### Code Organization
```python
# mypackage/__init__.py
from .core import MainClass
from .utils import helper_function

__all__ = ['MainClass', 'helper_function']
__version__ = '0.1.0'
```

#### Versioning
- Use semantic versioning: `MAJOR.MINOR.PATCH`
- Store version in `__init__.py` or `version.py`
- Update `CHANGELOG.md`

#### Documentation
```python
def function(param1, param2):
    """
    Brief description.
    
    Parameters
    ----------
    param1 : str
        Description of param1
    param2 : int
        Description of param2
    
    Returns
    -------
    result : DataFrame
        Description of return value
    
    Examples
    --------
    >>> function("test", 5)
    <result>
    """
    pass
```

#### Testing
```python
# tests/test_module.py
import pytest
from mypackage import function

def test_function():
    result = function("test", 5)
    assert result is not None
```

### 5. Comparison: Python vs R

| Feature | Python | R |
|---------|--------|---|
| **Package Dir** | `mypackage/` | `R/` |
| **Init File** | `__init__.py` | `NAMESPACE` |
| **Install Config** | `setup.py` | `DESCRIPTION` |
| **Dev Install** | `pip install -e .` | `devtools::load_all()` |
| **Testing** | `pytest` | `testthat` |
| **Documentation** | Sphinx / docstrings | roxygen2 |
| **Build** | `setup.py sdist` | `R CMD build` |
| **Upload** | PyPI (twine) | CRAN |

## 🔧 Development Tools

### Package Management
```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install package in dev mode
pip install -e .

# Install dev dependencies
pip install -e ".[dev]"
```

### Code Quality
```bash
# Format code
black geomx_tools/

# Lint code
flake8 geomx_tools/

# Type checking
mypy geomx_tools/
```

### Testing
```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=geomx_tools

# Run specific test
pytest tests/test_merge.py::test_function_name
```

### Documentation
```bash
# Build docs (Sphinx)
cd docs
make html

# View docs
open _build/html/index.html
```

## 📚 Resources

### Python Packaging
- [Python Packaging User Guide](https://packaging.python.org/)
- [setuptools documentation](https://setuptools.pypa.io/)
- [PyPA Sample Project](https://github.com/pypa/sampleproject)

### Development Tools
- [pytest](https://docs.pytest.org/) - Testing
- [black](https://black.readthedocs.io/) - Code formatting
- [flake8](https://flake8.pycqa.org/) - Linting
- [Sphinx](https://www.sphinx-doc.org/) - Documentation

### Modern Python Packaging
- [PEP 517](https://www.python.org/dev/peps/pep-0517/) - Build system
- [PEP 518](https://www.python.org/dev/peps/pep-0518/) - pyproject.toml
- [Poetry](https://python-poetry.org/) - Modern dependency management
- [Flit](https://flit.readthedocs.io/) - Simple packaging

## 🎯 Next Steps

1. ✅ Install package in dev mode: `pip install -e .`
2. ✅ Test imports: `from geomx_tools import merge_count_matrices`
3. ✅ Run pipeline notebook: See `ipynb/geomx_pipeline.ipynb`
4. ⏳ Add unit tests
5. ⏳ Add more documentation
6. ⏳ Publish to PyPI (optional)

---

**Version:** 0.1.0  
**Author:** Jaecheon Lee  
**License:** MIT (or your choice)

