#!/usr/bin/env bash
ulimit -c 0
g09 qc.in qc.out >> stdout.log &> stderr.log