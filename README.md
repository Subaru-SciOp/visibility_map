# A Python script to calculate observable time for HSC

## Overview
`calculate_observable_time.py` is a script designed to calculate the observable time for celestial objects based on specific criteria. This script is particularly intended for users of the HSC-queue observations to decide observing conditions (moon phase and moon separation).

`create_allskymap.py` is a script to generate a Bokeh interactive plot to show the calculated results.

## Definition of *Observable*
- Between sunset+HSC_DOMEOPEN_OFFSET and sunrise-HSC_DOMEOPEN_OFFSET
- Moon illumination < max_moon_illumination or moon altitude is below zero
- Moon separation to the object > min_moon_separation
- Object altitude is between MIN_ALTITUDE and MAX_ALTITUDE

In the code, the flag is defined as following.

```
flag_illum = moon_illum <= max_moon_illumination
flag_sep = moon_sep >= min_moon_separation * u.deg
flag_moon_alt = moon_alt < 0 * u.deg

flag_obj_alt = np.logical_and(obj_alt >= min_alt, obj_alt <= max_alt)
flag_observable = np.logical_or(flag_illum, flag_moon_alt)
flag_observable = np.logical_and(flag_observable, flag_sep)
flag_observable = np.logical_and(flag_observable, flag_obj_alt)
```

## Features
- Computes observable time windows for celestial objects.
- Takes into account various astronomical constraints.
- Outputs results in a user-friendly format.

## Requirements
- Python >3.10
- Required Python packages: `numpy`, `astropy`, `astroquery`, `healpy`, `joblib`, `pandas`, `bokeh`, `python-dotenv`, `colorcet`, and `uranography`

## Usage
To run the script, use the following command:
```console
$ python3 ./calculate_observable_time.py -h
usage: calculate_observable_time.py [-h] [--nside NSIDE] [--max-moon-illumination MAX_MOON_ILLUMINATION] [--min-moon-separation MIN_MOON_SEPARATION] observing_runs [outfile]

Calculate observable time for HSC

positional arguments:
  observing_runs        Observing runs csv file
  outfile               Output csv file

options:
  -h, --help            show this help message and exit
  --nside NSIDE         Healpix nside
  --max-moon-illumination MAX_MOON_ILLUMINATION
                        Max moon illumination
  --min-moon-separation MIN_MOON_SEPARATION
                        Min moon separation
```

Also, configuration can be made by `.env` file.
```
# observation site
OBSERVATION_SITE="subaru"

# time interval for calculation (min)
DELTA_TIME=10

# time offset before opening and closing the dome (min)
HSC_DOMEOPEN_OFFSET=30

# minimum and maximum altitude to be observed (degree)
MIN_ALTITUDE=30
MAX_ALTITUDE=90

# number of parallel jobs
N_JOBS=24
```


### Example
```sh
python calculate_observable_time.py observing_runs_s24b.csv output/observable_time.csv
```

Here, the input CSV file should look like the following.

```
run_begin,run_end
2024-08-01,2024-08-06
2024-09-24,2024-10-07
2024-11-25,2024-12-05
2024-12-19,2025-01-06
```
