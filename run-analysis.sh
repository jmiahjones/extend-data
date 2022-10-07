#!/bin/bash
export LD_PRELOAD=libopenblas.so
nice -n 15 Rscript ./R/analysis.R > logs/out.log 2>logs/message.log &&\
tail logs/out.log
