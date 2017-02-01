#!/bin/bash

NTRIALS=$1
MODULE=$2
LOCAL=$3
for trial in $(seq 1 ${NTRIALS}); do
    if [ "$LOCAL" == "local" ]; then
        sh trialTemplate.sh ${trial} ${MODULE}
    else
        qsub trialTemplate.sh -F "${trial} ${MODULE}"
    fi
done