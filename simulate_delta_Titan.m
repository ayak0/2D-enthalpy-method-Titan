clc
clear

%read the parameter file
alldata=readtable('data/scenario_params.csv');

%choose the scenario code where you want to simulate Titan deltas
scenario='V1S';

%set parameters for the selected scenario
data=alldata(strcmp(alldata.scenario, scenario), :); %import parameters for the selected scenario
period=data.period; %sea-level change period (see eq. 13 for dimensionless variables)
width=data.width; %channel width (100 m for Vid Flumina, 175 m for Saraswati Flumen)
Rab=data.Rab;

%additional parameters
A=0.005; %max sea level (50 m / length scale of 10000 m)
Tmax=1.0; %cycle times when the simulation ends

%simulate deltas
run('enthalpy_method_Titan.m')