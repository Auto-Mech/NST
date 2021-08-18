#!/usr/bin/env bash
ulimit -c 0
molpro -n 4 run.inp -o run.out --nouse-logfile --no-xml-output >> stdout.log &> stderr.log