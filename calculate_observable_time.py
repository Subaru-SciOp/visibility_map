#!/usr/bin/env python3

import argparse
import os
import sys
import time
import warnings
from contextlib import redirect_stderr
from datetime import timedelta

import healpy as hp
import joblib
import numpy as np
import pandas as pd
from astroplan import (
    AltitudeConstraint,
    FixedTarget,
    MoonIlluminationConstraint,
    MoonSeparationConstraint,
    Observer,
    time_grid_from_range,
)
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from dotenv import load_dotenv

warnings.filterwarnings("ignore")

# load environment variables
load_dotenv()

# set the default values for the environment variables
OBSERVATION_SITE = os.getenv("OBSERVATION_SITE", "subaru")
DELTA_TIME = float(os.getenv("DELTA_TIME", 5))
HSC_DOMEOPEN_OFFSET = float(os.getenv("HSC_DOMEOPEN_OFFSET", 30))
MIN_ALTITUDE = float(os.getenv("MIN_ALTITUDE", 30))
MAX_ALTITUDE = float(os.getenv("MAX_ALTITUDE", 90))
N_JOBS = float(os.getenv("N_JOBS", 8))

# the location of the observatory
observer = Observer.at_site(OBSERVATION_SITE)

# time interval for calculation
delta_time = DELTA_TIME * u.minute

# offset for dome open time
hsc_domeopen_offset = HSC_DOMEOPEN_OFFSET * u.minute

# min and max altitude for the object
min_alt, max_alt = MIN_ALTITUDE * u.deg, MAX_ALTITUDE * u.deg


# Function to calculate sunset and sunrise
def calculate_sun_times(date, observer):
    time = Time(date)
    sunset = observer.sun_set_time(time)
    sunrise = observer.sun_rise_time(time, which="next")
    return sunset, sunrise


# Function to calculate moon separation
def compute_constraints(
    target,
    times,
    max_moon_illumination: float = None,
    min_moon_separation: float = None,
):
    # compute the altitude constraint (true if the object is between the altitude limits)
    flag_obj_alt = AltitudeConstraint(min=min_alt, max=max_alt).compute_constraint(
        times, observer, target
    )

    # if the object is below or above the altitude limit, return None
    if np.all(~flag_obj_alt):
        return None, None, None

    # compute the moon separation constraint (true if the moon is further than min_moon_separation)
    with redirect_stderr(open(os.devnull, "w")):
        flag_moon_sep = MoonSeparationConstraint(
            min=min_moon_separation * u.deg
        ).compute_constraint(times, observer, target)

    # if the moon is always too close, return None
    if np.all(~flag_moon_sep):
        return None, None, None

    # moon rise and set conditions are considered in the astroplan.MoonIlluminationConstraint class
    # true if the moon illumination is less than max_moon_illumination
    flag_moon_illum = MoonIlluminationConstraint(
        max=max_moon_illumination
    ).compute_constraint(times, observer, target)

    return flag_moon_sep, flag_moon_illum, flag_obj_alt


def calculate_observable_time(
    targets,
    times,
    max_moon_illumination: float = None,
    min_moon_separation: float = None,
):
    # define a function to compute the observable time for a single target with joblib
    def _compute_constraints_to_joblib(target):
        flag_moon_sep, flag_moon_illum, flag_obj_alt = compute_constraints(
            target,
            times,
            max_moon_illumination=max_moon_illumination,
            min_moon_separation=min_moon_separation,
        )
        # when returned early because of no visibility, return 0.
        if flag_moon_sep is None:
            return 0.0

        # combine the constraints
        flag_observable = flag_moon_illum & flag_moon_sep & flag_obj_alt

        # return the total time the object is observable
        return len(times[flag_observable]) * delta_time.to(u.minute).value

    # parallel calculation with joblib
    parallel = joblib.Parallel(n_jobs=N_JOBS, return_as="generator")
    output_observable_time = parallel(
        joblib.delayed(_compute_constraints_to_joblib)(target) for target in targets
    )

    return np.array(list(output_observable_time), dtype=float)


def main_calculate_observable_time(
    df: pd.DataFrame,
    ra,
    dec,
    max_moon_illumination: float = 0.25,
    min_moon_separation: float = 30,
):

    t_start = time.time()

    # define targets from the healpix pixel coordinates
    targets = [
        FixedTarget(coord=SkyCoord(ra=r, dec=d, unit=(u.deg, u.deg), frame="icrs"))
        for r, d in zip(ra, dec)
    ]

    # array to store the observable time for each healpix pixel
    time_observable = np.zeros(ra.size) * u.minute

    # loop over observing runs
    for i in range(df.index.size):
        print(f"Start: run from {df['run_begin'][i]} to {df['run_end'][i]}")

        t_run_start = time.time()

        current_date = df["run_begin"][i]

        while current_date <= df["run_end"][i]:
            print(f"Start: night calculation for {current_date}")
            t_night_start = time.time()

            sunset, sunrise = calculate_sun_times(current_date, observer)

            times = time_grid_from_range(
                [sunset + hsc_domeopen_offset, sunrise - hsc_domeopen_offset],
                time_resolution=delta_time,
            )

            # calculate the observable time for each target (array with the length of the number of targets)
            res = calculate_observable_time(
                targets,
                times,
                max_moon_illumination=max_moon_illumination,
                min_moon_separation=min_moon_separation,
            )

            time_observable += res * u.minute

            t_night_end = time.time()

            print(
                f"Finished: night calculation for {current_date} ({t_night_end - t_night_start:.2f} seconds)"
            )
            current_date += timedelta(days=1)

        t_run_end = time.time()
        print(
            f"Finished: run from {df['run_begin'][i]} to {df['run_end'][i]} ({t_run_end - t_run_start:.2f} seconds)"
        )

    t_end = time.time()

    print(f"Calculation took {t_end - t_start:.2f} seconds")

    return time_observable.to(u.hour).value


def main(
    observing_runs: str,
    outfile: str,
    nside: int,
    max_moon_illumination: float = 0.25,
    min_moon_separation: float = 30,
):
    # Read observing_runs.csv
    df_runs = pd.read_csv(observing_runs, parse_dates=["run_begin", "run_end"])

    nside = args.nside
    npix = hp.nside2npix(nside)

    print(f"{npix=}")
    print(f"{hp.nside2resol(nside, arcmin=True)=}")

    # Generate ra, dec for each healpix pixel
    ra, dec = hp.pix2ang(nside, np.arange(npix), lonlat=True, nest=False)

    observable_times = main_calculate_observable_time(
        df_runs,
        ra,
        dec,
        max_moon_illumination=max_moon_illumination,
        min_moon_separation=min_moon_separation,
    )

    # create a pandas dataframe to store the results
    dfout = pd.DataFrame(
        {
            "ra": ra,
            "dec": dec,
            "time_observable": observable_times,
        }
    )

    if outfile == sys.stdout:
        print(dfout)
    else:
        dfout.to_csv(outfile, index=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculate observable time for HSC")

    parser.add_argument("observing_runs", help="Observing runs csv file")
    parser.add_argument(
        "outfile", help="Output csv file", nargs="?", default=sys.stdout
    )
    parser.add_argument(
        "--nside",
        help="Healpix nside",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--max-moon-illumination",
        help="Max moon illumination",
        type=float,
        default=0.25,
    )
    parser.add_argument(
        "--min-moon-separation",
        help="Min moon separation",
        type=float,
        default=30,
    )

    args = parser.parse_args()

    main(
        args.observing_runs,
        args.outfile,
        args.nside,
        max_moon_illumination=args.max_moon_illumination,
        min_moon_separation=args.min_moon_separation,
    )
