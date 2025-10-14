# GELATO Documentation

This directory contains the Sphinx documentation for the GELATO project.

## Building the Documentation

### Install Required Packages

```bash
pip install -r requirements.txt
```

### Build HTML Documentation

```bash
make html
```

The generated HTML documentation will be output to the `_build/html/` directory.

### View in Browser

```bash
# Linux/macOS
open _build/html/index.html

# Windows
start _build/html/index.html
```

## Other Formats

Sphinx supports various output formats:

- `make html` - HTML documentation
- `make latexpdf` - PDF documentation (requires LaTeX)
- `make epub` - ePub format
- `make man` - man pages
- `make text` - Plain text

To see all available targets:

```bash
make help
```

## Documentation Structure

- `conf.py` - Sphinx configuration file
- `Doxyfile` - Doxygen configuration file (for C++ documentation)
- `index.rst` - Top page
- `installation.rst` - Installation instructions
- `tutorial.rst` - Beginner-friendly tutorial
- `usage.rst` - Usage guide
- `examples.rst` - Detailed examples and sample explanations
- `modules.rst` - Module overview
- `api.rst` - Complete Python API reference
- `cpp_api.rst` - Complete C++ API reference
- `api/` - Individual module API documentation

## Automated Documentation Generation

### Python Documentation
This documentation is automatically generated from Python docstrings using `sphinx.ext.autodoc`.
When you add or edit docstrings in the code, rebuild the documentation to see the changes reflected.

### C++ Documentation
C++ documentation is automatically generated using **Doxygen + Breathe**:
1. Doxygen extracts documentation from C++ sources in the `src/` folder in XML format
2. Breathe (Sphinx extension) reads the XML and integrates it into Sphinx documentation
3. Running `make html` automatically executes Doxygen before building with Sphinx

#### C++ Documentation Coverage
- **Coordinate**: ECEF/ECI/NED coordinate transformations, quaternion operations
- **Air**: U.S. Standard Atmosphere model
- **Earth**: Earth model parameters (WGS84)
- **Gravity**: Gravity calculations
- **IIP**: Instantaneous Impact Point calculations

## Documentation Content

### Newly Created Documentation

#### tutorial.rst
Comprehensive beginner's tutorial guide:
- Step-by-step installation instructions
- Running basic examples
- Analyzing output results
- Customizing problem setup
- Understanding optimization behavior and troubleshooting

#### examples.rst
Detailed usage examples and reference:
- Detailed explanation of input files (JSON, CSV formats)
- Rocket configuration parameters
- Launch conditions and terminal conditions setup
- Flight constraint definitions
- Creating user-defined constraint functions
- Batch processing and parameter studies
- Visualization with KML files

### Existing Documentation

- **installation.rst**: Installation instructions (conda/pip, C++ build, troubleshooting)
- **usage.rst**: Basic usage and I/O file descriptions
- **modules.rst**: Module structure overview
- **api.rst**: Complete Python API reference
- **cpp_api.rst**: Complete C++ API reference (Doxygen integration)

## Build Status

✅ Documentation build successful (0 warnings)
✅ Python API documentation complete
✅ C++ API documentation complete (Doxygen + Breathe integration)
✅ Beginner tutorial added
✅ Detailed examples and reference added
✅ All modules documented in English

## GitHub Pages Deployment

This documentation is configured for automatic deployment to GitHub Pages:

### Automatic Deployment
- GitHub Actions workflow automatically builds and deploys documentation on push to `master` or `main` branch
- Documentation will be available at: `https://<username>.github.io/GELATO/`

### Manual Deployment
If you need to deploy manually:

```bash
cd docs
make html
# Push the _build/html/ directory contents to gh-pages branch
```

### Configuration Files
- `.github/workflows/docs.yml` - GitHub Actions workflow for automatic deployment
- `docs/.nojekyll` - Ensures GitHub Pages correctly serves Sphinx static files
- `docs/conf.py` - Contains GitHub Pages configuration
