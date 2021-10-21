#!/usr/bin/env bash
ulimit -c 0
g09 run.inp run.out >> stdout.log &> stderr.log