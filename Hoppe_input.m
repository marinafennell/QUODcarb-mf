% take Hoppe data and put it into new QUODcarb

load input.mat;
[in] = input;

alk = zeros(102,1); ealk = zeros(102,1); ralk = zeros(102,1);
dic = zeros(102,1); edic = zeros(102,1); rdic = zeros(102,1);
pHin = zeros(102,1); epH = zeros(102,1); 
Hin = zeros(102,1); eH = zeros(102,1); rH = zeros(102,1);
sal = zeros(102,1); esal = zeros(102,1);
temp = zeros(102,1); etemp = zeros(102,1);
pres = zeros(102,1); 
sil = zeros(102,1); esil = zeros(102,1);
po4 = zeros(102,1); epo4 = zeros(102,1);

for i = 1:102
    alk(i) = in(i,6)*1e-6; % mol/kg
      ealk(i) = in(i,14)*1e-6; % alkalinity absolute error
      ralk(i) = (alk(i) + ealk(i))/alk(i); % convert absolute error to relative errors
    dic(i) = in(i,7)*1e-6; % mol/kg
      edic(i) = in(i,15)*1e-6; % DIC error
      rdic(i) = (dic(i) + edic(i))/dic(i);
    pHin(i) = in(i,8); % total pH scale
    Hin(i) = 10^(-pHin(i));
    epH(i) = in(i,16); % pH error
      eH(i) = 10^(-pHin(i) - epH(i)) - 10^(-pHin(i));
      rH(i) = (Hin(i) + eH(i))/Hin(i);
    sal(i) = in(i,1); % salinity
      esal(i) = in(i,10);
    temp(i) = in(i,2); % deg C
      etemp(i) = in(i,11); 
    pres(i) = in(i,3); % dbar
    sil(i) = in(i,5)*1e-6; % total Si mol/kg
      esil(i) = in(i,13)*1e-6;
    po4(i) = in(i,4)*1e-6; % tota phosphate mol/kg
      epo4(i) = in(i,12)*1e-6;

end

LOG10 = log(10);
p = @(x) -log10(x); % use p as in "pH" to denote the negative log base 10 function
q = @(x) 10.^(-x);  % use q, i.e., a backward p

meas = NaN(102,26);
w = NaN(102,26);

meas(:,1) = dic(:); w(:,1) = rdic(:);
meas(:,2) = alk(:); w(:,2) = ralk(:);
meas(:,11) = Hin(:); w(:,11) = rH(:);


for i = 1:102

TC = temp(i); % deg C
S = sal(i);
P = pres(i); % dbar
yin = meas(i,:); % measured variables IN
win = w(i,:); % weights IN

[yobs,sigobs,sigx] = QUODcarb(yin,win,TC,S,P);




fprintf([' %7.1f +/- %6.1f  |%7.1f +/- %6.1f  |%7.4f +/- %6.4f  ::::   '...
  '%7.1f +/- %6.1f  |%7.1f +/- %6.1f   |%7.4f +/- %6.4f  |%7.1f +/- %6.1f \n'],  ...
        q(yobs(1)),sigobs(1),q(yobs(2)),sigobs(2),yobs(7),sigobs(7), ...
        q(x(1)),sigx(1), q(x(2)),sigx(2),x(7),sigx(7),q(x(6)),sigx(6));   

end






