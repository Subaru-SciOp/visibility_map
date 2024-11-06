#!/usr/bin/env python3

# %%
import argparse
import datetime
from collections import OrderedDict
from pprint import pprint

import bokeh
import bokeh.model
import colorcet as cc
import healpy as hp
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io.registry import IORegistryError
from astropy.table import Table
from bokeh.io import output_file
from bokeh.layouts import column, row
from bokeh.models import (
    ColorBar,
    CustomJS,
    Div,
    HoverTool,
    Label,
    Node,
)
from bokeh.palettes import MediumContrast3
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

    dict_observable_time = {
        "files": OrderedDict(),
        "moonsep": [],
        "moonphase": ["dark", "gray", "bright"],
        "vmin": vmin,
        "vmax": vmax,
        "cmap": cc.bmy,
        "title": f"{run_name.upper()} Observable Time by Moon Separation ({start_date.date()} - {end_date.date()})",
        "note": f"pixel size: {area.value:.2f} sq. deg.",
    }

    moonsep = [int(tb.meta["min_moon_separation"]) for tb in tables]
    dict_observable_time["moonsep"] = moonsep
    for i, moonsep in enumerate(moonsep):
        dict_observable_time["files"][f"moonsep{moonsep}"] = infiles[i]

    pprint(dict_observable_time["files"])

    vmin, vmax = dict_observable_time["vmin"], dict_observable_time["vmax"]

    # interactive plot
    plot = bokeh.plotting.figure(
        width=512 * 2,
        height=256 * 2 + 192,
        title=dict_observable_time["title"],
        match_aspect=True,
        x_axis_label="RA",
        y_axis_label="Dec",
        tools="pan,wheel_zoom,box_zoom,tap,undo,redo,reset",
        active_drag="box_zoom",
        output_backend="webgl",
    )
    sky = MollweideMap(plot=plot, location="Subaru")

    plot.title.text_font_size = "1.5em"

    cmap = bokeh.transform.linear_cmap(
        "value", dict_observable_time["cmap"], vmin, vmax
    )
    color_bar = ColorBar(
        color_mapper=cmap["transform"],
        bar_line_color="black",
        major_tick_line_color="black",
        title=f"Observable time in dark time with the minimum moon separation {dict_observable_time['moonsep'][0]}deg [h]",
        title_text_align="center",
    )

    plot.add_layout(color_bar, "below")

    # add a bar chart
    x_bar = [str(x) for x in dict_observable_time["moonsep"]]
    s_bar = bokeh.models.ColumnDataSource(
        data=dict(
            x=x_bar,
            dark=[0] * len(x_bar),
            gray=[0] * len(x_bar),
            bright=[0] * len(x_bar),
        )
    )
    p_bar = bokeh.plotting.figure(
        x_range=bokeh.models.FactorRange(*x_bar),
        width=512,
        height=380,
        x_axis_label="Minimum Moon Separation [deg]",
        y_axis_label="Observable time [h]",
        output_backend="webgl",
    )
    p_bar.vbar_stack(
        dict_observable_time["moonphase"],
        x="x",
        width=0.6,
        color=MediumContrast3,
        source=s_bar,
        legend_label=dict_observable_time["moonphase"],
    )
    p_bar.y_range.start = 0
    p_bar.legend.location = "top_center"
    p_bar.legend.orientation = "horizontal"
    leg_bar = p_bar.legend[0]
    p_bar.add_layout(leg_bar, "below")
    h_bar = HoverTool(
        tooltips=[("Dark", "@dark"), ("Gray", "@gray"), ("Bright", "@bright")]
    )
    p_bar.add_tools(h_bar)

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

        if i == 0:
            hpix_ds, hp_cmap, hp_glyph = sky.add_healpix(
                df_agg[dict_observable_time["moonphase"][0]].to_numpy(),
                nside=nside,
                cmap=cmap,
            )
            hover_tool = sky.add_hp_hovertool(coordinates=True, value=None)

        for i_mp in range(len(dict_observable_time["moonphase"])):
            hpix_ds.data[
                f"{dict_observable_time['moonphase'][i_mp]}_{dict_observable_time['moonsep'][i]}"
            ] = df_agg[dict_observable_time["moonphase"][i_mp]].to_numpy()

            hover_tool.tooltips.append(
                (
                    f"T({dict_observable_time['moonphase'][i_mp]}, {dict_observable_time['moonsep'][i]}deg) [h]",
                    f"@{dict_observable_time['moonphase'][i_mp]}_{dict_observable_time['moonsep'][i]}{{0.0}}",
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

    frame_left = Node(target="frame", symbol="left", offset=5)
    frame_bottom = Node(target="frame", symbol="bottom", offset=-5)

    note = Label(
        x=frame_left,
        y=frame_bottom,
        anchor="bottom_left",
        text=dict_observable_time["note"],
        padding=5,
    )

    plot.add_layout(note)

    plot.min_border_left = 48
    plot.min_border_top = 96
    plot.min_border_bottom = 96

    plot.js_on_event(
        "tap",
        CustomJS(
            args={
                "hpix_ds": hpix_ds,
                "s_bar": s_bar,
                "moonphase": dict_observable_time["moonphase"],
                "n_moonphase": len(dict_observable_time["moonphase"]),
                "moonsep": dict_observable_time["moonsep"],
                "n_moonsep": len(dict_observable_time["moonsep"]),
            },
            code="""
            const data = s_bar.data;
            if (hpix_ds.selected.indices.length == 1) {
                var selected_index = hpix_ds.selected.indices[0];
                s_bar.data = {
                    x: data['x'],
                    dark: [hpix_ds.data['dark_' + data['x'][0]][selected_index],
                            hpix_ds.data['dark_' + data['x'][1]][selected_index],
                            hpix_ds.data['dark_' + data['x'][2]][selected_index]],
                    gray: [hpix_ds.data['gray_' + data['x'][0]][selected_index],
                            hpix_ds.data['gray_' + data['x'][1]][selected_index],
                            hpix_ds.data['gray_' + data['x'][2]][selected_index]],
                    bright: [hpix_ds.data['bright_' + data['x'][0]][selected_index],
                            hpix_ds.data['bright_' + data['x'][1]][selected_index],
                            hpix_ds.data['bright_' + data['x'][2]][selected_index]],
                }
            }
            if (hpix_ds.selected.indices.length == 0) {
                s_bar.data = {
                    x: data['x'],
                    dark: [0, 0, 0],
                    gray: [0, 0, 0],
                    bright: [0, 0, 0],
                }
            }
            """,
        ),
    )

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
        margin=(0, 0, 0, 32),
        width=512,
    )

    grid = column(plot, row(p_bar, div))

    output_file(outfile, title=dict_observable_time["title"])

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
