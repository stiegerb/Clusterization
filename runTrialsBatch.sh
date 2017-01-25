#!/bin/bash

NTRIALS=$1
for trial in $(seq 1 ${NTRIALS}); do
    qsub trialTemplate.sh -F ${trial}
done