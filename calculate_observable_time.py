#!/usr/bin/env python3

import argparse
import os
import sys
import time
import warnings
from contextlib import redirect_stderr
from datetime import datetime, timedelta

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
from astropy.table import Table
from astropy.time import Time
from dotenv import load_dotenv

warnings.filterwarnings("ignore")

# load environment variables
load_dotenv()

# set the default values for the environment variables
OBSERVATION_SITE = os.getenv("OBSERVATION_SITE", "subaru")
UTC_OFFSET = float(os.getenv("UTC_OFFSET", -10))
DELTA_TIME = float(os.getenv("DELTA_TIME", 5))
DOMEOPEN_OFFSET = float(os.getenv("DOMEOPEN_OFFSET", 30))
MIN_ALTITUDE = float(os.getenv("MIN_ALTITUDE", 30))
MAX_ALTITUDE = float(os.getenv("MAX_ALTITUDE", 90))
N_JOBS = float(os.getenv("N_JOBS", 8))

# the location of the observatory
observer = Observer.at_site(OBSERVATION_SITE)

# UTC offset for the observatory
utc_offset = UTC_OFFSET * u.hour

# time interval for calculation
delta_time = DELTA_TIME * u.minute

# offset for dome open time
domeopen_offset = DOMEOPEN_OFFSET * u.minute

# min and max altitude for the object
min_alt, max_alt = MIN_ALTITUDE * u.deg, MAX_ALTITUDE * u.deg


# Function to calculate sunset and sunrise
def calculate_sun_times(date, observer):
    # assume date is local time and convert to UTC by adding the UTC offset
    time = Time(date) + utc_offset
    sunset = observer.sun_set_time(time)
    sunrise = observer.sun_rise_time(time, which="next")
    return sunset, sunrise


# Function to calculate moon separation
def compute_constraints(
    target,
    times,
    min_moon_illumination=None,
    max_moon_illumination=None,
    min_moon_separation: float = 30.0,
):
    if min_moon_illumination is None:
        min_moon_illumination = [None, 0.25, 0.65]
    if max_moon_illumination is None:
        max_moon_illumination = [0.25, 0.65, None]

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
    flag_moon_illum = []
    for min_illum, max_illum in zip(min_moon_illumination, max_moon_illumination):
        flag_moon_illum.append(
            MoonIlluminationConstraint(min=min_illum, max=max_illum).compute_constraint(
                times, observer, target
            )
        )

    return flag_obj_alt, flag_moon_sep, flag_moon_illum


def calculate_observable_time(
    targets,
    times,
    min_moon_illumination=None,
    max_moon_illumination=None,
    min_moon_separation: float = None,
):
    if min_moon_illumination is None:
        min_moon_illumination = [None, 0.25, 0.65]
    if max_moon_illumination is None:
        max_moon_illumination = [0.25, 0.65, None]

    # define a function to compute the observable time for a single target with joblib
    def _compute_constraints_to_joblib(target):
        flag_obj_alt, flag_moon_sep, flag_moon_illum = compute_constraints(
            target,
            times,
            min_moon_illumination=min_moon_illumination,
            max_moon_illumination=max_moon_illumination,
            min_moon_separation=min_moon_separation,
        )

        # when returned early because of no visibility, return 0.
        if flag_moon_sep is None:
            return [0.0, 0.0, 0.0]

        # combine the constraints
        flag_observable_dark = flag_moon_illum[0] & flag_moon_sep & flag_obj_alt
        flag_observable_gray = flag_moon_illum[1] & flag_moon_sep & flag_obj_alt
        flag_observable_bright = flag_moon_illum[2] & flag_moon_sep & flag_obj_alt

        # output list
        delta_t = delta_time.to(u.minute).value
        res = [
            float(len(times[flag_observable_dark]) * delta_t),
            float(len(times[flag_observable_gray]) * delta_t),
            float(len(times[flag_observable_bright]) * delta_t),
        ]

        # return the total time the object is observable
        return res

    # parallel calculation with joblib, parallelized over the targets, i.e., (ra, dec)
    parallel = joblib.Parallel(n_jobs=N_JOBS, return_as="generator")
    output_joblib_observable_time = parallel(
        joblib.delayed(_compute_constraints_to_joblib)(target) for target in targets
    )
    df_observable_time = pd.DataFrame(
        {
            "dark": np.zeros(len(targets)),
            "gray": np.zeros(len(targets)),
            "bright": np.zeros(len(targets)),
        }
    )

    for i, t in enumerate(list(output_joblib_observable_time)):
        # print(t)
        df_observable_time["dark"][i] += t[0]
        df_observable_time["gray"][i] += t[1]
        df_observable_time["bright"][i] += t[2]

    # print(df_observable_time)

    # print(len(observable_time))

    return df_observable_time


def main_calculate_observable_time(
    dates,
    ra,
    dec,
    min_moon_illumination=None,
    max_moon_illumination=None,
    min_moon_separation: float = 30,
):

    if min_moon_illumination is None:
        min_moon_illumination = [None, 0.25, 0.65]
    if max_moon_illumination is None:
        max_moon_illumination = [0.25, 0.65, None]

    t_start = time.time()

    # define targets from the healpix pixel coordinates
    targets = [
        FixedTarget(coord=SkyCoord(ra=r, dec=d, unit=(u.deg, u.deg), frame="icrs"))
        for r, d in zip(ra, dec)
    ]

    # array to store the observable time for each healpix pixel
    date_out = np.empty((dates.size, ra.size), dtype=type(dates[0]))
    ipix_out = np.zeros((dates.size, ra.size), dtype=int)
    ra_out = np.zeros((dates.size, ra.size))
    dec_out = np.zeros((dates.size, dec.size))
    time_observable_dark = np.zeros((dates.size, ra.size))
    time_observable_gray = np.zeros((dates.size, ra.size))
    time_observable_bright = np.zeros((dates.size, ra.size))

    # loop over observing runs
    for i in range(dates.size):
        t_date_start = time.time()

        date = dates[i]

        date_out[i, :] = [date.strftime("%Y-%m-%d")] * ra.size

        ipix_out[i, :] = np.arange(ra.size)
        ra_out[i, :] = ra
        dec_out[i, :] = dec

        print(f"Start: night calculation for {date.strftime('%Y-%m-%d')}")

        sunset, sunrise = calculate_sun_times(date, observer)

        times = time_grid_from_range(
            [sunset + domeopen_offset, sunrise - domeopen_offset],
            time_resolution=delta_time,
        )

        # calculate the observable time for each target (array with the length of the number of targets (ra, dec))
        df_obs_time_date = calculate_observable_time(
            targets,
            times,
            min_moon_illumination=min_moon_illumination,
            max_moon_illumination=max_moon_illumination,
            min_moon_separation=min_moon_separation,
        )

        time_observable_dark[i, :] = df_obs_time_date["dark"]
        time_observable_gray[i, :] = df_obs_time_date["gray"]
        time_observable_bright[i, :] = df_obs_time_date["bright"]

        t_date_end = time.time()
        print(
            f"Finished: night calculation for {date.strftime('%Y-%m-%d')} ({t_date_end - t_date_start:.2f} seconds)"
        )

    t_end = time.time()

    print(f"Calculation took {t_end - t_start:.2f} seconds")

    tbret = Table(
        {
            "date": date_out.flatten(),
            "ipix": ipix_out.flatten(),
            "ra": ra_out.flatten() * u.deg,
            "dec": dec_out.flatten() * u.deg,
            "dark": (time_observable_dark.flatten() * u.minute).to(u.hour),
            "gray": (time_observable_gray.flatten() * u.minute).to(u.hour),
            "bright": (time_observable_bright.flatten() * u.minute).to(u.hour),
        }
    )
    return tbret


def main(
    start_date: str,
    end_date: str,
    outfile: str,
    nside: int,
    min_moon_illumination: float = None,
    max_moon_illumination: float = 0.25,
    min_moon_separation: float = 30,
):
    # create a list of dates for the observing runs
    dates = pd.date_range(start=start_date, end=end_date, freq="1D", inclusive="both")

    # healpix arrays
    nside = args.nside
    npix = hp.nside2npix(nside)

    print(f"{npix=}")
    print(f"{hp.nside2resol(nside, arcmin=True)=}")

    # Generate ra, dec for each healpix pixel
    ra, dec = hp.pix2ang(nside, np.arange(npix), lonlat=True, nest=False)

    tbout = main_calculate_observable_time(
        dates,
        ra,
        dec,
        min_moon_illumination=min_moon_illumination,
        max_moon_illumination=max_moon_illumination,
        min_moon_separation=min_moon_separation,
    )

    tbout.meta["site"] = OBSERVATION_SITE
    tbout.meta["nside"] = nside
    tbout.meta["date_start"] = start_date
    tbout.meta["date_end"] = end_date
    tbout.meta["delta_time"] = DELTA_TIME
    tbout.meta["min_moon_illumination"] = min_moon_illumination
    tbout.meta["max_moon_illumination"] = max_moon_illumination
    tbout.meta["min_moon_separation"] = min_moon_separation
    tbout.meta["min_altitude"] = min_alt
    tbout.meta["max_altitude"] = max_alt

    if outfile == sys.stdout:
        print(tbout)
    else:
        tbout.write(
            outfile,
            # formats={"date": "%Y-%m-%d"},
            overwrite=True,
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculate observable time")

    # parser.add_argument("observing_runs", help="Observing runs csv file")
    parser.add_argument(
        "start_date", help="Start date of the observing run (e.g., 2024-01-01)"
    )
    parser.add_argument(
        "end_date", help="End date of the observing run, inclusive (e.g., 2024-12-31)"
    )
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
        "--min-moon-illumination",
        nargs=3,
        help="Min moon illumination (dark, gray, bright), default is [None, 0.25, 0.65]",
        type=list,
        default=[None, 0.25, 0.65],
    )
    parser.add_argument(
        "--max-moon-illumination",
        nargs=3,
        help="Max moon illumination (dark, gray, bright), default is [0.25, 0.65, None]",
        type=list,
        default=[0.25, 0.65, None],
    )
    parser.add_argument(
        "--min-moon-separation",
        help="Min moon separation in degrees, default is 30",
        type=float,
        default=30,
    )

    args = parser.parse_args()

    main(
        # args.observing_runs,
        args.start_date,
        args.end_date,
        args.outfile,
        args.nside,
        min_moon_illumination=args.min_moon_illumination,
        max_moon_illumination=args.max_moon_illumination,
        min_moon_separation=args.min_moon_separation,
    )
