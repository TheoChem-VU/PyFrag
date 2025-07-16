# PyFrag Cleanup and Simplification Recommendations

## 1. Eliminate Code Duplication in Backend Implementations

### Current State
- Multiple identical copies of `PyFrag.py` across backends (6+ files with 90%+ duplicate code)
- Identical argument parsing logic repeated in each backend
- Multiple copies of `PyFragModules.py` with identical functions
- Shell scripts (`pyfragparce.sh`) with identical parsing functions

### Recommended Solution
Create a unified architecture:

```
src/pyfrag/
├── core/
│   ├── __init__.py
│   ├── argument_parser.py      # Unified argument parsing
│   ├── modules.py              # Consolidated PyFragModules
│   ├── input_parser.py         # Unified input file parsing
│   └── job_runner.py           # Main execution logic
├── backends/
│   ├── __init__.py
│   ├── base.py                 # Base backend class
│   ├── adf_new.py             # ADF-specific logic only
│   ├── adf_old.py             # Legacy ADF logic only
│   ├── gaussian.py            # Gaussian-specific logic only
│   └── orca.py                # ORCA-specific logic only
└── cluster_runner.py          # Already exists and good
```

**Benefits**:
- Reduce codebase size by ~60-70%
- Single point of maintenance for argument parsing
- Easier to add new backends
- Consistent behavior across all backends

## 2. Unify Input File Parsing

### Current State
- Multiple shell scripts with identical `grep -A 200` patterns
- Duplicate functions: `jobsubargue`, `pyfragargue`, `adfargue`
- Complex shell-based parsing that's hard to maintain

### Recommended Solution
Create a Python-based input parser:

```python
# src/pyfrag/core/input_parser.py
class InputParser:
    def __init__(self, input_file):
        self.input_file = input_file
        self.sections = {}

    def parse_section(self, section_name, end_marker=None):
        """Parse a section between markers like JOBSUB...JOBSUB END"""

    def extract_arguments(self):
        """Extract all PyFrag arguments from input file"""

    def generate_job_script(self, backend, template):
        """Generate cluster job submission script"""
```

**Benefits**:
- Eliminate 8+ shell script files
- More robust parsing (handles edge cases better)
- Easier testing and debugging
- Python-based = better error handling

## 3. Consolidate Argument Parsing

### Current State
Each `PyFrag.py` has 40+ identical argument definitions:
```python
parser.add_argument("--ircpath", type=str, action='append', nargs='*', help='IRC coordinate file')
parser.add_argument("--fragment", type=int, action='append', nargs='*', help='atom number for each fragment')
# ... repeated in 6+ files
```

### Recommended Solution
```python
# src/pyfrag/core/argument_parser.py
class PyFragArgumentParser:
    @staticmethod
    def create_parser():
        parser = argparse.ArgumentParser(description='PyFrag Fragment Analysis')

        # Coordinate inputs
        coord_group = parser.add_argument_group('Coordinate inputs')
        coord_group.add_argument("--ircpath", type=str, action='append', nargs='*')
        # ... all arguments defined once

        return parser

    @staticmethod
    def process_arguments(args, backend_type):
        """Process arguments with backend-specific logic"""
```

**Benefits**:
- Single source of truth for all arguments
- Easier to add new analysis options
- Consistent help text and validation
- Backend-specific argument processing

## 4. Simplify Shell Script Architecture

### Current State
- Complex shell functions duplicated across files
- Mix of shell and Python parsing
- Hard to debug and maintain

### Recommended Solution
```bash
# Simplified pyfragparce.sh (one file for all backends)
#!/bin/bash
export HOSTPYFRAG=${HOSTPYFRAG:-$(pwd)}
python ${HOSTPYFRAG}/src/pyfrag/core/input_parser.py "$@"
```

All parsing logic moves to Python for better maintainability.

## 5. Modernize the Package Structure

### Current State (pyproject.toml is good, but structure needs work)
- Mixed old and new Python packaging patterns
- Some legacy scripts in root directories
- Inconsistent import patterns

### Recommended Structure
```
src/pyfrag/
├── __init__.py
├── main.py                    # Entry point (already exists)
├── cluster_runner.py          # Entry point (already exists)
├── core/
│   ├── __init__.py
│   ├── argument_parser.py
│   ├── input_parser.py
│   ├── modules.py             # Consolidated from PyFragModules.py
│   └── job_runner.py
├── backends/
│   ├── __init__.py
│   ├── base.py
│   └── [specific backends]
├── parsers/                   # Move from qmworks/parsers/
│   ├── __init__.py
│   └── [specific parsers]
└── utils/
    ├── __init__.py
    └── [utility functions]
```

## 6. Eliminate Deprecated/Unused Code

### Files to Consider Removing (after user confirmation)
- `host/standalone/adf_newopen/` (seems to be unused)
- `host/standalone/adf_openorb/` (deprecated as mentioned)
- Multiple duplicate PLAMS installations in different backends
- Legacy shell scripts in `host/argueparce/`

## 7. Improve Configuration Management

### Current State
- Hard-coded paths in shell scripts
- Environment variables scattered across files

### Recommended Solution
```python
# src/pyfrag/core/config.py
class PyFragConfig:
    def __init__(self):
        self.pyfrag_home = os.environ.get('PYFRAGHOME', self.detect_home())
        self.host_pyfrag = os.environ.get('HOSTPYFRAG', self.detect_host())

    def get_backend_path(self, backend_name):
        return os.path.join(self.pyfrag_home, 'backends', f'{backend_name}.py')
```

## 8. Create Unified Entry Points

### Update pyproject.toml console scripts
```toml
[project.scripts]
pyfrag = "pyfrag.main:main"                    # Keep existing
pyfrag-cluster = "pyfrag.cluster_runner:main" # Keep existing
pyfrag-parse = "pyfrag.core.input_parser:main" # New unified parser
```

## Implementation Priority

### Phase 1 (High Impact, Low Risk)
1. Create unified argument parser
2. Consolidate PyFragModules.py files
3. Create configuration management

### Phase 2 (Medium Risk)
1. Create unified input parser (Python-based)
2. Refactor shell scripts to use Python parser
3. Create backend base classes

### Phase 3 (Structural Changes)
1. Eliminate duplicate backend files
2. Restructure package layout
3. Remove deprecated code

## Estimated Impact
- **Code Reduction**: ~60-70% fewer lines
- **Maintainability**: Single point of change for common functionality
- **Reliability**: Python parsing is more robust than shell scripts
- **Testability**: Much easier to write unit tests
- **Performance**: Minimal impact, possibly faster startup

## Backward Compatibility
All changes can be implemented while maintaining the existing command-line interface:
- Same `pyfrag job.in` commands work
- Same input file format
- Same output format
- Virtual environment integration preserved
