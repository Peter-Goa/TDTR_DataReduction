# TDTR_DataReduction
TDTR Data process by Aaron Jerome Schmidt's method. This program can calculate an unlimited-layer structure and can fit unlimited number of parameters at the same time.

Usage:
This code has several modes to process TDTR experament datas:

* Fitting Mode - to get values of unknown parameters by fitting the experimental data.

* Theory Curve Mode - to draw the theory curve at a special experiment condition.

* Sensitivity Mode - to calculate sensitivity of parameters special experiment condition.

* Two-Frequency Fitting - to get a pair of parameters by fitting the data from two modulation 
frequencies at the same time.

* Uncertainty Mode - to calculate the uncertainty of parameters at a special experiment condition.


First, you need to modify the configuration file (there are examples with a suffix of `.tdtrcfg` at `./Configration/Example` for each mode seperatelly). And there are two ways to run this programm.

## First way:
1. copy and modify the configuration file to match yourself's experiment condition
2. run TDTRDataProcess.m
3. select configuration file and other paths follow the program prompts

## Second way:
1. copy and modify the configuration file to match yourself's experiment condition.
2. run TDTRDataProcess(ConfigurationFileAddress, ...)
