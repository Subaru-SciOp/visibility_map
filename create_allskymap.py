#!/usr/bin/env python3

# %%
import argparse
import os
from collections import OrderedDict
from pprint import pprint

import bokeh
import colorcet as cc
import healpy as hp
import numpy as np
import pandas as pd
from astropy import units as u
from bokeh.embed import file_html
from bokeh.io import output_file, output_notebook, show
from bokeh.layouts import column, gridplot, row
from bokeh.models import ColorBar, CustomJS, Div, Label, Node, NumericInput, Paragraph
from uranography.api import MollweideMap


# %%
def generate_plot(datadir, run_name=None, nside=32, notebook=True):

    area = hp.nside2pixarea(nside, degrees=True) * u.deg**2  # in square degrees

    print(f"{area=}")

    dict = {
        "files": OrderedDict(),
        "moonsep": [],
        "moonphase": [],
        "vmin": 0,
        "vmax": 310,
        # "cmap": cc.CET_CBL3,
        # "cmap": cc.CET_L8,
        "cmap": cc.bmy,
        "title": f"{run_name.upper()} HSC Observable Time",
        "note": f"pixel size: {area.value:.2f} sq. deg.",
    }

    for moonphase in ["dark", "gray"]:
        for moonsep in [30, 60, 90]:
            dict["files"][f"moonsep{moonsep}_{moonphase}"] = os.path.join(
                datadir,
                f"observable_time_{run_name}_nside{nside}_moonsep{moonsep}_{moonphase}.csv",
            )
            dict["moonphase"].append(moonphase)
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
        title="Observable time (h)",
        title_text_align="center",
    )

    plot.add_layout(color_bar, "below")

    for i, tup in enumerate(dict["files"].items()):
        k, v = tup
        print(i, k, v)

        if i == 0:
            try:
                df_main = pd.read_csv(v)
            except FileNotFoundError:
                print(f"File not found: {v}")
                raise FileNotFoundError(f"File not found {v}")

            hpix_ds, hp_cmap, hp_glyph = sky.add_healpix(
                df_main["time_observable"].to_numpy(), nside=nside, cmap=cmap
            )

            hover_tool = sky.add_hp_hovertool(coordinates=True, value=None)
            hover_tool.tooltips.append(
                (
                    f"T({dict['moonphase'][i]}, moonsep={dict['moonsep'][i]}deg) [h]",
                    "@value{0.0}",
                )
            )
        else:
            try:
                df_comp = pd.read_csv(v)
            except FileNotFoundError:
                print(f"File not found: {v}")
                continue

            hpix_ds.data[f"count_{k}"] = df_comp["time_observable"].to_numpy()

            hover_tool.tooltips.append(
                (
                    f"T({dict['moonphase'][i]}, moonsep={dict['moonsep'][i]}deg) [h]",
                    f"@count_{k}{{0.0}}",
                )
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
        text=dict["note"],
        padding=5,
    )

    plot.add_layout(note)

    plot.min_border_left = 48
    plot.min_border_top = 96
    plot.min_border_bottom = 96

    grid = column(plot)

    if notebook:
        output_notebook()
        show(grid, notebook_handle=True)
    else:
        output_file(
            f"hsc_observable_time_{run_name}_nside{nside}.html", title=dict["title"]
        )

    print(f"Saving to hsc_observable_time_{run_name}_nside{nside}.html")

    saved_file = bokeh.plotting.save(grid)


# %%
if __name__ == "__main__":
    try:
        get_ipython
        print("running in a notebook")

        generate_plot(
            "output2",
            "hsc_observable_time_s24a_nside32.html",
            run_name="s24a",
            nside=32,
            notebook=True,
        )
        generate_plot(
            "output2",
            "hsc_observable_time_s24b_nside32.html",
            run_name="s24b",
            nside=32,
            notebook=True,
        )
    except NameError:

        parser = argparse.ArgumentParser(description="Create plot")

        parser.add_argument(
            "datadir",
            type=str,
            help="Data directory. Filenames must follow the pattern: observable_time_{{run_name}}_nside{{nside}}_moonsep{{moonsep}}_{{moonphase}}.csv",
        )
        # parser.add_argument("outfile", type=str, help="Output file")
        parser.add_argument(
            "--run_name", type=str, default="s24a", help="Run name (default: s24a)"
        )
        parser.add_argument(
            "--nside", type=int, default=32, help="Healpix nside (default: 32)"
        )
        parser.add_argument("--notebook", action="store_true", help="Notebook mode")

        args, unknown = parser.parse_known_args()

        print(f"running with args: {args}")
        generate_plot(
            args.datadir,
            run_name=args.run_name,
            nside=args.nside,
            notebook=args.notebook,
        )


# %%
