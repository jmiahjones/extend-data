#!/bin/bash
export LD_PRELOAD=libopenblas.so
nice -n 15 Rscript ./R/analysis.R > logs/out.log 2>logs/message.log &&\
tail logs/out.log &&\
nice -n 15 Rscript ./R/analysis-alt.R > logs/out.log 2>logs/message-alt.log &&\
tail logs/out-alt.log &&\
nice -n 15 Rscript ./R/no-noise.R > logs/out.log 2>logs/message-nonoise.log &&\
tail logs/out-nonoise.log
