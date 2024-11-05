#!/usr/bin/env python3

# %%
import argparse
import datetime
import os
from collections import OrderedDict
from pprint import pprint

import bokeh
import colorcet as cc
import healpy as hp
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io.registry import IORegistryError
from astropy.table import Table
from bokeh.embed import file_html
from bokeh.io import output_file, output_notebook, show
from bokeh.layouts import column, gridplot, row
from bokeh.models import ColorBar, CustomJS, Div, Label, Node, NumericInput, Paragraph
from uranography.api import MollweideMap


# %%
def generate_plot(
    infiles, outfile, start_date, end_date, run_name=None, vmin=0, vmax=650
):

    start_date = datetime.datetime.strptime(start_date, "%Y-%m-%d")
    end_date = datetime.datetime.strptime(end_date, "%Y-%m-%d")

    print(start_date, end_date)

    tables = []

    if np.all(
        [
            infile.endswith(".parquet")
            or infile.endswith(".ecsv.gz")
            or infile.endswith(".ecsv")
            for infile in infiles
        ]
    ):
        for infile in infiles:
            print(f"Reading {infile}")
            try:
                tb = Table.read(infile)
            except IORegistryError:
                tb = Table.read(infile, format="ascii.ecsv")
            tables.append(tb)
    else:
        raise ValueError("Input file must be in parquet, ecsv.gz, or ecsv format.")

    nside = tables[0].meta["nside"]
    area = hp.nside2pixarea(nside, degrees=True) * u.deg**2  # in square degrees

    print(f"{area=}")

    dict = {
        "files": OrderedDict(),
        "moonsep": [],
        "moonphase": ["dark", "gray", "bright"],
        # "vmin": 0,
        # "vmax": 310,
        "vmin": vmin,
        "vmax": vmax,
        # "cmap": cc.CET_CBL3,
        # "cmap": cc.CET_L8,
        "cmap": cc.bmy,
        "title": f"{run_name.upper()} Observable Time by Moon Separation ({start_date.date()} - {end_date.date()})",
        "note": f"pixel size: {area.value:.2f} sq. deg.",
    }
    # print(f"{len(tables)=}")

    # for moonphase in ["dark", "gray", "bright"]:
    moonsep = [int(tb.meta["min_moon_separation"]) for tb in tables]
    for i, moonsep in enumerate(moonsep):
        # print(i, moonsep)
        dict["files"][f"moonsep{moonsep}"] = infiles[i]
        # dict["moonphase"].append(moonphase)
        dict["moonsep"].append(moonsep)

    pprint(dict["files"])

    vmin, vmax = dict["vmin"], dict["vmax"]

    # interactive plot
    plot = bokeh.plotting.figure(
        width=512 * 2,
        height=256 * 2 + 192,
        title=dict["title"],
        match_aspect=True,
        x_axis_label="RA",
        y_axis_label="Dec",
        tools="pan,wheel_zoom,box_zoom,undo,redo,reset",
        active_drag="box_zoom",
        output_backend="webgl",
    )
    sky = MollweideMap(plot=plot, location="Subaru")

    plot.title.text_font_size = "1.5em"

    cmap = bokeh.transform.linear_cmap("value", dict["cmap"], vmin, vmax)
    color_bar = ColorBar(
        color_mapper=cmap["transform"],
        # border_line_color="black",
        bar_line_color="black",
        major_tick_line_color="black",
        title=f"Observable time in dark time with the minimum moon separation {dict['moonsep'][0]}deg [h]",
        title_text_align="center",
    )

    plot.add_layout(color_bar, "below")

    for i in range(len(tables)):

        df = tables[i].to_pandas()
        df["date"] = pd.to_datetime(df["date"])

        df_agg = (
            df.loc[
                (df["date"] >= start_date) & (df["date"] <= end_date),
                ["ipix", "dark", "gray", "bright"],
            ]
            .groupby("ipix")
            .sum()
        )

        # print(df_agg.head(50))

        if i == 0:

            hpix_ds, hp_cmap, hp_glyph = sky.add_healpix(
                df_agg[dict["moonphase"][0]].to_numpy(), nside=nside, cmap=cmap
            )
            hpix_ds.data[f"{dict['moonphase'][1]}_{dict['moonsep'][i]}"] = df_agg[
                dict["moonphase"][1]
            ].to_numpy()
            hpix_ds.data[f"{dict['moonphase'][2]}_{dict['moonsep'][i]}"] = df_agg[
                dict["moonphase"][2]
            ].to_numpy()

            hover_tool = sky.add_hp_hovertool(coordinates=True, value=None)
            hover_tool.tooltips.append(
                (
                    f"T({dict['moonphase'][0]}, {dict['moonsep'][i]}deg) [h]",
                    "@value{0.0}",
                ),
            )
            hover_tool.tooltips.append(
                (
                    f"T({dict['moonphase'][1]}, {dict['moonsep'][i]}deg) [h]",
                    f"@{dict['moonphase'][1]}_{dict['moonsep'][i]}{{0.0}}",
                ),
            )
            hover_tool.tooltips.append(
                (
                    f"T({dict['moonphase'][2]}, {dict['moonsep'][i]}deg) [h]",
                    f"@{dict['moonphase'][2]}_{dict['moonsep'][i]}{{0.0}}",
                ),
            )
        else:
            # pass
            #
            for i_mp in range(len(dict["moonphase"])):
                hpix_ds.data[f"{dict['moonphase'][i_mp]}_{dict['moonsep'][i]}"] = (
                    df_agg[dict["moonphase"][i_mp]].to_numpy()
                )

                hover_tool.tooltips.append(
                    (
                        f"T({dict['moonphase'][i_mp]}, {dict['moonsep'][i]}deg) [h]",
                        f"@{dict['moonphase'][i_mp]}_{dict['moonsep'][i]}{{0.0}}",
                    ),
                )

    sky.add_graticules(
        graticule_kwargs={
            "min_decl": -90,
            "max_decl": 90,
            "decl_space": 30,
            "min_ra": 0,
            "max_ra": 360,
            "ra_space": 30,
        },
        line_kwargs={"color": "white"},
    )

    # plot.axis.visible = True
    # plot.xaxis.ticker = [0, 60, 120, 180, 240, 300]
    # plot.xaxis.major_label_overrides = {
    #     0: "0",
    #     60: "60",
    #     120: "120",
    #     180: "180",
    #     240: "240",
    # }

    frame_left = Node(target="frame", symbol="left", offset=5)
    frame_bottom = Node(target="frame", symbol="bottom", offset=-5)

    note = Label(
        x=frame_left,
        y=frame_bottom,
        anchor="bottom_left",
        text=dict["note"],
        padding=5,
    )

    plot.add_layout(note)

    plot.min_border_left = 48
    plot.min_border_top = 96
    plot.min_border_bottom = 96

    div = Div(
        text="""<font size='5'>
        Definition for the moon phase:
        </font><br />
        <font size='3'>
        Dark: Moon illumination less than 0.25<br />
        Gray: Moon illumination between 0.25 and 0.65<br />
        Bright: Moon illumination more than 0.65<br />
        Note that the moon phase is different from <a href="https://www.naoj.org/Observing/Proposals/howto.html">the definition of Subaru Telescope</a>.
        </font>
""",
        margin=(0, 0, 0, 48),
        # min_border_left=48,
        # min_border_top=96,
        # min_border_bottom=96,
    )

    grid = column(plot, div)

    output_file(outfile, title=dict["title"])

    saved_file = bokeh.plotting.save(grid)


# %%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create plot")

    parser.add_argument(
        "input_files",
        nargs="+",
        type=str,
        help="Input files (parquet, ecsv, or ecsv.gz)",
    )
    parser.add_argument("output_file", type=str, help="Output file (html)")
    parser.add_argument(
        "start_date",
        type=str,
        help="Start date (YYYY-MM-DD)",
    )
    parser.add_argument(
        "end_date",
        type=str,
        help="End date (YYYY-MM-DD)",
    )
    parser.add_argument(
        "--run-name", type=str, default="s24a", help="Run name (default: s24a)"
    )

    parser.add_argument(
        "--vmin", type=int, default=0, help="Minimum value for colorbar"
    )
    parser.add_argument(
        "--vmax", type=int, default=650, help="Maximum value for colorbar"
    )

    args, unknown = parser.parse_known_args()

    print(args)

    print(f"running with args: {args}")
    generate_plot(
        args.input_files,
        args.output_file,
        args.start_date,
        args.end_date,
        run_name=args.run_name,
        vmin=args.vmin,
        vmax=args.vmax,
    )


# %%
