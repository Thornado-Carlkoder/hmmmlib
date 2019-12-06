#!/bin/bash

mkdir running_times
python running_time_A1.py > running_times/running_time_A1.csv
python running_time_A2.py > running_times/running_time_A2.csv

python running_time_C1.py > running_times/running_time_C1.csv
python running_time_C2.py > running_times/running_time_C2.csv
python running_time_C3.py > running_times/running_time_C3.csv
python running_time_C4.py > running_times/running_time_C4.csv
