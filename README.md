# A Python script to calculate observable time by moon phase and moon separation

## Overview

`calculate_observable_time.py` is a script designed to calculate the observable time for celestial objects based on specific criteria (moon phase and moon separation).

`create_allskymap.py` is a script to generate a Bokeh interactive plot to show the calculated results.

## Definition of *Observable*

- Between `sunset + DOMEOPEN_OFFSET` and `sunrise - DOMEOPEN_OFFSET`
- Moon illumination between `min_moon_illumination` and `max_moon_illumination` for each moon phase or moon altitude is below zero.
- Moon separation to the object is larger than `min_moon_separation`
- Object altitude is between `MIN_ALTITUDE` and `MAX_ALTITUDE`

## Features

- Computes observable time windows for celestial objects.
- Takes into account various astronomical constraints.
- Outputs results in a user-friendly format.

## Requirements

- Python >3.10
- Required Python packages: `numpy`, `astropy`, `astroquery`, `healpy`, `joblib`, `pandas`, `pyarrow`, `bokeh`, `python-dotenv`, `colorcet`, and `uranography`

## Usage

To run the script, use the following command:
```console
$ python3 ./calculate_observable_time.py -h
usage: calculate_observable_time.py [-h] [--nside NSIDE] [--min-moon-illumination MIN_MOON_ILLUMINATION MIN_MOON_ILLUMINATION MIN_MOON_ILLUMINATION]
                                    [--max-moon-illumination MAX_MOON_ILLUMINATION MAX_MOON_ILLUMINATION MAX_MOON_ILLUMINATION] [--min-moon-separation MIN_MOON_SEPARATION]
                                    start_date end_date [outfile]

Calculate observable time

positional arguments:
  start_date            Start date of the observing run (e.g., 2024-01-01)
  end_date              End date of the observing run, inclusive (e.g., 2024-12-31)
  outfile               Output file (astropy.table compatible format, parquet recommended)

options:
  -h, --help            show this help message and exit
  --nside NSIDE         Healpix nside
  --min-moon-illumination MIN_MOON_ILLUMINATION MIN_MOON_ILLUMINATION MIN_MOON_ILLUMINATION
                        Min moon illumination (dark, gray, bright), default is [None, 0.25, 0.65]
  --max-moon-illumination MAX_MOON_ILLUMINATION MAX_MOON_ILLUMINATION MAX_MOON_ILLUMINATION
                        Max moon illumination (dark, gray, bright), default is [0.25, 0.65, None]
  --min-moon-separation MIN_MOON_SEPARATION
                        Min moon separation in degrees, default is 30
```

Also, configuration can be made by `.env` file.

```conf
# observation site
OBSERVATION_SITE="subaru"

# Offset between UTC and the local time at the observation site (h)
UTC_OFFSET=-10

# time interval for calculation (min)
DELTA_TIME=10

# time offset before opening and closing the dome after/before sunset and sunrise (min)
DOMEOPEN_OFFSET=30

# minimum and maximum altitude to be observed (degree)
MIN_ALTITUDE=30
MAX_ALTITUDE=90

# number of parallel jobs
N_JOBS=8
```

To create an all skymap based on the reqult file, `create_allskymap.py` will do the job.

```console
$ python3 ./create_allskymap.py -h
usage: create_allskymap.py [-h] [--run-name RUN_NAME] [--vmin VMIN] [--vmax VMAX] input_files [input_files ...] output_file start_date end_date

Create plot

positional arguments:
  input_files          Input files (parquet, ecsv, or ecsv.gz)
  output_file          Output file (html)
  start_date           Start date (YYYY-MM-DD)
  end_date             End date (YYYY-MM-DD)

options:
  -h, --help           show this help message and exit
  --run-name RUN_NAME  Run name (default: s24a)
  --vmin VMIN          Minimum value for colorbar
  --vmax VMAX          Maximum value for colorbar
```

### Example

```sh
python3 ./calculate_observable_time.py 2025-09-01 2025-09-10 tmp.parquet
python3 ./create_allskymap.py tmp.parquet observable_time_tmp.html 2025-09-01 2025-09-10 --run-name="test" --vmin=0 --vmax=100
```
