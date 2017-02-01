#!/bin/bash

NTRIALS=$1
MODULE=$2
for trial in $(seq 1 ${NTRIALS}); do
    qsub trialTemplate.sh -F "${trial} ${MODULE}"
    #sh trialTemplate.sh ${trial}
done