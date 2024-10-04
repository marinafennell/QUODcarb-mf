
% QUODcarb example 1
clear all
% measured TC and TA only

% populate opt structure, see readme for default settings
opt.phscale = 1;        % which phscale are we interested in
opt.printcsv = 1;       % print output to CSV? 1 = on, 0 = off
opt.fname = 'example1.csv'; % file name of desired CSV file


% populate obs structure
obs.TC = 2150; % umol/kg-SW
obs.uTC = 5; % ± 5 umol/kg-SW, 1 sigma
obs.TA = 2300; % umol/kg-SW
obs.uTA = 5; % ± 5 umol/kg-SW, 1 sigma
obs.sal = 32.7; % PSU
obs.usal = 0.02; % ± 0.02 PSU, 1 sigma
obs.tp(1).T = 20; % deg Celsius
obs.tp(1).uT = 0.001; % deg Celsius, 1 sigma
obs.tp(1).P = 0; % dbar
obs.tp(1).uP = 0.005; % ± 0.005 dbar, 1 sigma

[est,obs,sys,iflag,opt] = QUODcarb(obs,opt);

% outputs should be:
% est.TC = 2150                 est.uTC = 4.9942
% est.TA = 2300                 est.uTA = 4.9946
% est.tp(1).ph = 7.8395         est.tp(1).uph = 0.0257
% est.tp(1).pco2 = 706.3685     est.tp(1).upco2 = 158.2209
% est.tp(1).co3 = 118.4753      est.tp(1).uco3 = 5.8916


% uncomment if want to compare output to CO2SYS's output:
% compare(obs,est,opt,1,1,1,'compare.csv');


