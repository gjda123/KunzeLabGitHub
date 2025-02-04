%% Spike Detection for Kunze Lab
% Use this code for Spike Detection for electrophysiology.

clear
clc

addpath('functions');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%File to analyze
Parameters.Filename=[]; 
%^Leave blank (Parameters.Filename=[]) if you want to select a file with UI, 
% otherwise include full path and file ID
    
    
[Parameters,Data] = load_MEA(Parameters);
[Parameters,Data] = filterElectrodes(Parameters,Data);

%%%% SET YOUR FALLING THRESHOLD FOR SPIKE DETECTION HERE %%%%
Parameters.standard_deviation=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Parameters,Data] = SpikeDetection(Parameters,Data);

scatter(Data.SpikeOutput(:,1),Data.SpikeOutput(:,2));