#!/bin/bash
export LD_PRELOAD=libopenblas.so
nice -n 15 Rscript ./R/stage-2-selection-alt.R > logs/out.log 2>logs/message.log &&\
tail logs/out.log
