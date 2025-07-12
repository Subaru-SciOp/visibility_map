# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an astronomical observatory visibility mapping tool for Subaru Telescope. The project calculates observable time windows for celestial objects based on moon phase, moon separation, and altitude constraints. It consists of two main Python scripts:

- `calculate_observable_time.py`: Core calculation engine that computes observable time using HEALPix sky maps
- `create_allskymap.py`: Visualization tool that generates interactive Bokeh plots from calculation results

## Development Environment

The project uses:
- **Python**: Requires Python >=3.11 
- **Package manager**: `uv` (preferred) with `uv.lock` file, fallback to `pip` with `requirements.txt`
- **Virtual environment**: Located at `.venv/` (already configured)

### Setup Commands

```bash
# Install dependencies with uv (preferred)
uv sync

# Or install with pip
pip install -r requirements.txt
```

## Core Workflows

### Running Calculations

Individual calculation:
```bash
python3 calculate_observable_time.py START_DATE END_DATE OUTPUT_FILE.parquet [options]
```

Batch calculations (multiple moon separations):
```bash
./run_calculations.bash
```

### Creating Visualizations

```bash
python3 create_allskymap.py INPUT_FILES OUTPUT_FILE.html START_DATE END_DATE [options]
```

### Configuration

Observatory and calculation parameters are configured via `.env` file:
- `OBSERVATION_SITE`: Observatory location (default: "subaru")
- `UTC_OFFSET`: Time zone offset in hours (default: -10)
- `MIN_ALTITUDE`/`MAX_ALTITUDE`: Altitude constraints in degrees
- `N_JOBS`: Parallel processing jobs count

## Development Tools

Code quality tools are configured in `pyproject.toml`:

```bash
# Code formatting
black .

# Type checking  
mypy .

# Linting
flake8 .

# Testing
pytest
```

## Architecture

### Data Pipeline
1. **Input**: Date ranges and astronomical constraints
2. **Processing**: HEALPix-based sky tessellation with parallel computation
3. **Output**: Parquet files containing observable time maps
4. **Visualization**: Interactive Bokeh plots using uranography library

### Key Dependencies
- **astroplan/astropy**: Astronomical calculations and coordinate systems
- **healpy**: HEALPix sky map operations
- **joblib**: Parallel processing
- **bokeh**: Interactive web-based plotting
- **uranography**: Astronomical visualization framework (LSST project)
- **pandas/pyarrow**: Data processing and Parquet I/O

### File Structure
- `data/`: Pre-calculated parquet files for different semesters and moon separations
- `html/`: Generated visualization outputs
- `uranography_example/`: Sample data and notebook demonstrations
- `run_calculations.bash`: Batch processing script for multiple moon separation values

## Observatory-Specific Details

Configured for **Subaru Telescope** operations:
- Location: Maunakea, Hawaii (HST timezone, UTC-10)
- Typical observing constraints: 30-90° altitude, various moon separation requirements
- Semester-based scheduling (S24A, S24B, S25A, S25B format)
- Standard moon separation values: 30°, 45°, 60°, 90°

## Data Formats

- **Input**: Date strings (YYYY-MM-DD format)
- **Processing**: HEALPix NSIDE=32 resolution (12,288 sky pixels)
- **Output**: Parquet files with columns for pixel indices and observable hours
- **Visualization**: HTML files with interactive sky maps