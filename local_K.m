% QUODcarb_v2
% starting 04/2022


% calculate equilibrium constants-------------------------------------
function [Kh,K1,K2,Kb,Kw,Ks,Kf,K1p,K2p,K3p,Ksi,p2f] = local_K(T,S,P)
    % T = Temp (deg C) input
    % P = pressure (dbar)    

    TK = T + 273.15; % convert to Kelvin
    Rgas = 83.1451; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
    RT = Rgas.*TK;
    Pbar = P./10;
    IonS = 19.924 .* S ./ (1000-1.005 .* S); % from DOE handbook

    % pCO2 to fCO2 conversion
    Pstd = 1.01325; 
    delC = (57.7 - 0.118.*TK); 
    B = -1636.75 + 12.0408.*TK - 0.0327957.*TK.^2 + 3.16528.*0.00001.*TK.^3;
    p2f = exp((B + 2.*delC).*Pstd./(RT));

    % calculate TF (Riley 1965)--------------------------------------
    TF = (0.000067./18.998).*(S./1.80655); % mol/kg-SW

    % calculate TS (Morris & Riley 1966)-----------------------------
    TS = (0.14./96.062).*(S./1.80655); % mol/kg-SW
    
    % calculate Kh (Weiss 1974)--------------------------------------
    TK100 = TK./100;
    lnKh = -60.2409 + 93.4517 ./ TK100 + 23.3585 .* log(TK100) + S .* ...
        (0.023517 - 0.023656 .* TK100 + 0.0047036 .* TK100 .^2);
    Kh = exp(lnKh);

    % calculate Ks (Dickson 1990a)----------------------------------
    lnKs = -4276.1./TK + 141.328 - 23.093 .* log(TK) + ...
        (-13856./TK + 324.57 - 47.986 .* log(TK)) .* sqrt(IonS) + ...
        (35474./TK - 771.54 + 114.723 .* log(TK)) .* IonS + ...
        (-2698./TK) .* sqrt(IonS) .* IonS + (1776./TK) .* IonS.^2;
    Ks = exp(lnKs) .* (1 - 0.001005 .* S); % converted to mol/kg-SW

    % calculate Kf (Dickson 1979)----------------------------------
    lnKf = 1590.2/TK - 12.641 + 1.525 .* IonS.^0.5;
    Kf = exp(lnKf) .* (1 - 0.001005 .* S); % converted to mol/kg-SW
   
    % pH scale conversion factors (not pressure corrected)-----------
    SWS2tot = (1 + TS ./Ks)./(1 + TS ./Ks + TF./Kf);
    FREE2tot = 1 + TS./Ks;
    
    % calculate fH (Takahashi et al 1982)--------------------------
    fH = 1.2948 - 0.002036 .* TK + (0.0004607 - ...
        0.000001475 .* TK) .* S.^2 ;

    % calculate Kb (Dickson 1990)----------------------------------
    lnKbt = -8966.9 - 2890.53 .* sqrt(S) - 77.942 .* S + ...
        1.728 .* sqrt(S) .* S - 0.0996 .* S.^2;
    lnKb = lnKbt./ TK + 148.0248 + 137.1942 .* sqrt(S) + ...
        1.62142 .* S + (-24.4344 - 25.085 .* sqrt(S) - 0.2474 .* ...
        S) .* log(TK) + 0.053105 .* sqrt(S) .* TK;
    Kb = exp(lnKb)./SWS2tot; % SWS pH scale, mol/kg-SW

    % calculate Kw (Millero 1995)--------------------------------
    lnKw = 148.9802 - 13847.26 ./ TK - 23.6521 .* log(TK) + ...
        (-5.977 + 118.67 ./ TK + 1.0495 .* log(TK)) .* ...
        sqrt(S) - 0.01615 .* S;
    Kw = exp(lnKw);

    % calculate K1p, K2p, K3p, Ksi (Yao and Millero 1995)-------
    lnK1p = -4576.752 ./TK + 115.54 - 18.453 .* log(TK) + ...
        (-106.736./TK + 0.69171) .* sqrt(S) + (-0.65643./TK - 0.01844).*S;
    K1p = exp(lnK1p);

    lnK2p = -8814.715./TK + 172.1033 - 27.927.*log(TK) + ...
        (-160.34./TK + 1.3566).*sqrt(S) + (0.37335./TK - 0.05778).*S;
    K2p = exp(lnK2p);

    lnK3p = -3070.75./TK - 18.126 + (17.27039./TK + 2.81197).*sqrt(S) + ...
        (-44.99486./TK - 0.09984).*S;
    K3p = exp(lnK3p);

    lnKsi = -8904.2./TK + 117.4 - 19.334.*log(TK) + (-458.79./TK + ...
        3.5913).*sqrt(IonS) + (188.74/TK - 1.5998).*IonS + ...
        (-12.1652./TK + 0.07871).*IonS.^2;
    Ksi = exp(lnKsi).*(1 - 0.001005.*S); % convert to mol/kg-SW

    % calculate K1 & K2 (Mehrbach refit by Dickson and Millero 1987)---
    pK1 = 3670.7 ./TK - 62.008 + 9.7944 .* log(TK) - 0.0118.*S + ...
        0.000116.*S.^2;
    K1 = 10.^(-pK1); % SWS pH scale in mol/kg-SW

    pK2 = 1394.7./TK + 4.777 - 0.0184.*S + 0.000118.*S.^2;
    K2 = 10.^(-pK2); % SWS pH scale in mol/kg-SW

    % corrections for pressure---------------------------------------
    % sources: Millero 1995, 1992, 1982, 1979; Takahashi et al. 1982;
    %   Culberson & Pytkowicz 1968; Edmond & Gieskes 1970.

    dV = -25.5 + 0.1271.*T;
    Ka = (-3.08 + 0.0877 .* T) ./1000;
    lnK1fac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on K1

    dV = -15.82 - 0.0219 .* T;
    Ka = (1.13 - 0.1475 .*T)./1000;
    lnK2fac = (-dV + 0.5.*Ka .*Pbar).*Pbar./RT; % pressure effect on K2

    dV = -29.48 + 0.1622.*T - 0.002608.*T.^2;
    Ka = -2.84./1000;
    lnKbfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on Kb

    dV = -20.02 + 0.1119.*T - 0.001409 .*T.^2;
    Ka = (-5.13 + 0.0794.*T)./1000;
    lnKwfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on Kw

    dV = -9.78 - 0.009.*T - 0.000942.*T.^2;
    Ka = (-3.91 + 0.054.*T)./1000;
    lnKffac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on Kf

    dV = -18.03 + 0.0466.*T + 0.000316.*T.^2;
    Ka = (-4.53 + 0.09.*T)./1000;
    lnKsfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on Ks

    dV = -14.51 + 0.1211.*T - 0.000321.*T.^2;
    Ka  = (-2.67 + 0.0427.*T)./1000;
    lnK1pfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on K1p

    dV = -23.12 + 0.1758.*T - 0.002647.*T.^2;
    Ka  = (-5.15 + 0.09  .*T)./1000;
    lnK2pfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on K2p

    dV = -26.57 + 0.202 .*T - 0.003042.*T.^2;
    K  = (-4.08 + 0.0714.*T)./1000;
    lnK3pfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on K3p

    dV = -29.48 + 0.1622.*T - 0.002608.*T.^2;
    K  = -2.84./1000;
    lnKsifac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on Ksi

    % correct all Ks for pressure effects
    K1fac  = exp(lnK1fac);  K1  = K1 .*K1fac;
    K2fac  = exp(lnK2fac);  K2  = K2 .*K2fac;
    Kwfac  = exp(lnKwfac);  Kw  = Kw .*Kwfac;
    Kbfac  = exp(lnKbfac);  Kb  = Kb .*Kbfac;
    Kffac  = exp(lnKffac);  Kf  = Kf .*Kffac;
    Ksfac  = exp(lnKsfac);  Ks  = Ks .*Ksfac;
    K1pfac = exp(lnK1pfac); K1p = K1p.*K1pfac;
    K2pfac = exp(lnK2pfac); K2p = K2p.*K2pfac;
    K3pfac = exp(lnK3pfac); K3p = K3p.*K3pfac;
    Ksifac = exp(lnKsifac); Ksi = Ksi.*Ksifac;

    % CorrectpHScaleConversionsForPressure:
    % fH has been assumed to be independent of pressure.
    SWS2tot  = (1 + TS./Ks)./(1 + TS./Ks + TF./Kf);
    FREE2tot =  1 + TS./Ks;

end