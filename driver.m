
% driver to go with QUODcarb

load data.mat;
[in] = data;
nD = length(in);

% choose options for opt structure
opt.K1K2 = 10; % option for K1K2 formulation
opt.KSO4 = 1;  % option for KSO4 formulation
opt.KF   = 2;  % option for KF formulation
opt.TB   = 2;  % option for TB formulation
opt.phscale  = 1;  % 1 = tot, 2 = sws, 3 = free, 4 = NBS
opt.printcsv = 0;  % print est to CSV? 1 = on , 0 = off
% opt.fname    = 'QUODcarb_output.csv'; % don't need it if printcsv is off
opt.fname    = 'output_csv/test.csv';
opt.co2press = 1; % 1 = on, 0 = off
opt.Revelle  = 0; % 1 = on, 0 = off 
opt.printmes = 0; % 1 = on, 0 = off


% read in GOMECC data and put into obs structure
for i = 1:nD   
    % measurements that are independent of (T,P)
    obs(i).TC    = in(5,i); % (umol/kg)
    obs(i).eTC   = 2.01;    % TC error ±2.01 umol/kg
    obs(i).TA    = in(6,i);
    obs(i).eTA   = 1.78; % TA error ±2.01 umol/kg
    obs(i).sal   = in(1,i); % PSU
    obs(i).esal  = 0.001; % 1 new as of 1/23 old = 0.002
    % nutrients P and Si also independent of (T,P)
    obs(i).TP    = in(7,i);
    obs(i).eTP   = 0.0019; % 0.30% meas precision & mean(TP) = 0.6408
    obs(i).TSi   = in(8,i);
    obs(i).eTSi  = 0.0238; % 0.31% meas uncertainty & mean(TSi) = 7.68

    % first (T,P)-dependent measurement for pH
    obs(i).tp(1).T    = 25 ; % degC
    obs(i).tp(1).eT   = 0.05 ; % from cruise report
    obs(i).tp(1).P    = 0.0 ; %in(i+ad,1); % NOT in situ
    obs(i).tp(1).eP   = 0.07 ;
    obs(i).tp(1).ph   = in(9,i); % total scale
    obs(i).tp(1).eph  = 0.001 ; % pg 60, cruise report

    % second (T,P)-dependent measurement for pCO2
    obs(i).tp(2).T     = 20 ; % degC
    obs(i).tp(2).eT    = 0.03 ; % from cruise report
    obs(i).tp(2).P     = 0.0 ; % dbar (surface pressure for pco2)
    obs(i).tp(2).eP    = 0.07 ;
    obs(i).tp(2).pco2  = in(10,i); % (µatm)
    obs(i).tp(2).epco2 = 1.1353; % 0.21% relative std error & avg(pco2) = 540.6128

    % third (T,P)-dependent measurement for CO32-T
    obs(i).tp(3).T    = 25 ; % degC
    obs(i).tp(3).eT   = 0.05 ; % from cruise report
    obs(i).tp(3).P    = 0.0 ; % NOT in situ
    obs(i).tp(3).eP   = 0.07 ;
    obs(i).tp(3).co3  = in(11,i); % (µmol/kg)
    obs(i).tp(3).eco3 = in(11,i)*0.02;  % 2% from Jon Sharp NEW 1/25/24

    % fourth (T,P)-dependent measurement IN SITU
    obs(i).tp(4).T  = in(2,i); % deg C, CTD temp
    obs(i).tp(4).eT = 0.02; % ±0.02 degC
    obs(i).tp(4).P  = in(3,i); % dbar
    obs(i).tp(4).eP = 0.63; % (max) ± 0.63 dbar
end

obs_backup = obs;

%% Q5: All five input
% CT AT pH pCO2 CO3 (Q5) (fid5)
[est,obs,sys,iflag,opt] = QUODcarb(obs,opt);
    % note there are now five outputs from QUODcarb




    