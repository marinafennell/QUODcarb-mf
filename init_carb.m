% make initial guess for QUODcarb using CO2SYS
clearvars;
format long

par1type = 2; % DIC
par1 = 2073.5; % DIC_in µmol/kg
epar1 = 10.92; % error DIC
par2type = 3; % pH
par2 = 8.08; % pH(total)
epar2 = 0.018; % error pH
sal = 37.5; % salinity of sample
esal = 0.04; % salinity error
tempin = 15; % temp of sample
etemp = 0.02; % error temp
presin = 0; % pressure of sample (dbars)
tempout = tempin;
presout = presin;
sil = 2.5; % total Si of sample
esi = 0.02; % error Si
po4 = 0.86; % phosphate of sample
epo4 = 0.02; % error PO4
pHscale = 1; % pH(total scale)NBS=4
k1k2c = 4; % Mehrbach refit
epK = [0.002,0.01,0.02,0.01,0.01,0.02,0.02]; % '' = [0.02, 0.0075, 0.015, 0.01, 0.01, 0.02, 0.02]
kso4c = 1; % Dickson 1990
eBt = ''; % error total boron, 0.02

out(:) = CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,...
                       presin,presout,sil,po4,pHscale,k1k2c,kso4c);

Kh = out(66);
K1 = out(67);
K2 = out(68);    
Kb = out(72);    
Kw = out(71);
Ks = out(74);    
Kf = out(73);
K1p = out(75);
K2p = out(76);    
K3p = out(77);   
Ksi = out(78);

CT = out(02)*1e-6; % mol/kgSW
AT = out(01)*1e-6; % mol/kgSW
BT = out(79)*1e-6; % mol/kgSW
co2 = out(23)*1e-6; % mol/kgSW
hco3 = out(21)*1e-6; % mol/kgSW
co3 = out(22)*1e-6; % mol/kgSW
boh4 = out(24)*1e-6; % mol/kgSW
boh3 = BT - boh4; % mol/kgSW
fco2 = out(20)*1e-6; % atm
oh = out(25)*1e-6; % mol/kgSW
pH = out(37);
h = 10^(-pH); % [H+]
TS = out(81)*1e-6; % mol/kgSW
TF = out(80)*1e-6; % mol/kgSW
TP = out(82)*1e-6; % mol/kgSW
TSi = out(83)*1e-6; % mol/kgSW

hf = h / (1+ (TS/Ks));
HF = TF / (1+ (Kf/hf));
f = TF - HF; 
hso4 = TS / (1+ (Ks/hf));
so4 = TS - hso4;   
den = ((h^3)+(K1p*(h^2))+(K1p*K2p*h)+(K1p*K2p*K3p)); % denominator
po4 = TP * K1p * K2p * K3p / (den); 
hpo4 = TP * K1p * K2p * h / (den); 
h2po4 = TP * K1p * (h^2) / (den); 
h3po4 = TP * (h^3) / (den);  
siooh3 = TSi / (1 + (h/Ksi));
sioh4 = TSi - siooh3;

init_c = [CT; AT; BT; co2; hco3; co3; boh3; boh4; fco2; oh; h; ...
          TS; TF; TP; TSi; hf; HF; f; so4; hso4; po4; hpo4; ...
          h2po4; h3po4; siooh3; sioh4];

save 'init_c' init_c.mat