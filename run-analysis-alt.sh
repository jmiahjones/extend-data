#!/bin/bash
export LD_PRELOAD=libopenblas.so
nice -n 15 Rscript ./R/analysis-alt.R > logs/out.log 2>logs/message-alt.log &&\
tail logs/out-alt.log
