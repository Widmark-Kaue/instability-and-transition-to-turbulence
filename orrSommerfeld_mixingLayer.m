%%
clc; clear all; close all
path(path, 'src')

%% Constans
Re = 10000; 
alpha = 1;
poiseuilliFlow = @(y) 1-y.^2;