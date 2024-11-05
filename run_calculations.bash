#!/bin/bash

NSIDE=32
DATADIR="data"
START_DATE="2024-02-01"
END_DATE="2025-01-31"

for moonsep in 30 60 90 45; do
    time python3 -u ./calculate_observable_time.py $START_DATE $END_DATE ${DATADIR}/observable_time_nside${NSIDE}_moonsep${moonsep}_${START_DATE//-/}-${END_DATE//-/}.parquet \
        --nside=32 --min-moon-separation=${moonsep} 2>&1 | tee log.nside${NSIDE}_moonsep${moonsep}_${START_DATE//-/}-${END_DATE//-/}.parquet
    date
    echo ""
done
