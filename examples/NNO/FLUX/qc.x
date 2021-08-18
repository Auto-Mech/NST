#!/usr/bin/env bash
ulimit -c 0
molpro -n 4 qc.in -o qc.out --nouse-logfile --no-xml-output >> stdout.log &> stderr.log