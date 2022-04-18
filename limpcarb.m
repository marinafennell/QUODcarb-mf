%
%

function [f,g,H] = limpcarb(z,y,w,TC,S,P)
% [f,g,H] = limpco2(z,y,w,TC,S,P)
%
% Negative log probability for the co2 system  a.k.a. log improbability i.e., limp!
%
% INPUT:
%
%   z  := [pCT;pAT;pBT;...;psioh4;pKh;pK1;...;pKsi;lam1;...;lam26] where pX := -log10(x)
%   y  := measured components, same size and order as z(1:37), non-measured components set to NaN
%   TC := temperature in deg C
%   S  := Salinity in PSU
%   P  := Pressure in dbar
%   w  := measurement precisions, same size and order as y, non-measured components set to 0 or NaN
%
% OUTPUT:
%
%   f := limp
%   g := grad f w.r.t. x
%   h := hessian of f w.r.t. x

    
% utility functions and constants
    LOG10 = log(10);
    
    p = @(x) -log10( x );  % inverse of q
    
    q = @(x)  10.^( -x );  % inverse of p
    
    dqdx = @(x) - LOG10 * 10.^( -x );  
    d2qdx2 = @(x) -LOG10^2 * 10.^( -x );
    
    x   =  z( 1 : 37);  % measureable variables
    lam =  z( 38 : 63);  % Lagrange multipliers

    %
    % "Mass"-conservation equations (x7)
   
    % mc1 = CT -(co2  +hco3  +co3);
    mc1 = [ 1  0  0 -1 -1 -1  0  0  0  0  0 ...
           0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ...
           0  0  0  0  0  0  0  0  0  0  0 ];
    % mc2 = AT -(hco3 +2*co3 +boh4  +oh + hpo4 + 2*po4 + siooh3 - hf - hso4 - HF - h3po4);
    mc2 = [ 0  1  0  0 -1 -2  0 -1  0 -1  0 ...
           0  0  0  0 +1 +1  0  0 +1 -2 -1  0 +1 -1  0 ...
           0  0  0  0  0  0  0  0  0  0  0 ];
    % mc3 = BT  -(boh3 +boh4);
    mc3 = [ 0  0  1  0  0  0 -1 -1  0  0  0 ...
           0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ...
           0  0  0  0  0  0  0  0  0  0  0 ];
    % mc4 = TS  -(hso4 +so4);
    mc4 = [ 0  0  0  0  0  0  0  0  0  0  0 ...
          +1  0  0  0  0  0  0 -1 -1  0  0  0  0  0  0 ...
           0  0  0  0  0  0  0  0  0  0  0 ];
    % mc5 = TF  -(HF   +f);
    mc5 = [ 0  0  0  0  0  0  0  0  0  0  0 ...
           0 +1  0  0  0 -1 -1  0  0  0  0  0  0  0  0 ...
           0  0  0  0  0  0  0  0  0  0  0 ];
    % mc6 = TP  -(h3po4 +h2po4 +hpo4 +po4);
    mc6 = [ 0  0  0  0  0  0  0  0  0  0  0 ...
           0  0 +1  0  0  0  0  0  0 -1 -1 -1 -1  0  0 ...
           0  0  0  0  0  0  0  0  0  0  0 ];
    % mc7 = TSi -(sioh4 +siooh3);
    mc7 = [ 0  0  0  0  0  0  0  0  0  0  0 ...
           0  0  0 +1  0  0  0  0  0  0  0  0  0 -1 -1 ...
           0  0  0  0  0  0  0  0  0  0  0 ];
    M = zeros(7,37);
    M(1,:) = mc1;
    M(2,:) = mc2;
    M(3,:) = mc3;
    M(4,:) = mc4;
    M(5,:) = mc5;
    M(6,:) = mc6;
    M(7,:) = mc7;
    %M = [mc1; mc2; mc3; mc4; mc5; mc6; mc7];
    
    %
    % Equilibrium constants:

    % pKh = (pco2 - pfco2));          % Kh = [CO2]/fCO2 
    ec1 = [ 0   0   0   1   0   0   0   0  -1   0   0   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0];
    % pK1 = (ph + phco3 - pco2);      % K1 = [HCO3][H]/[CO2]
    ec2 = [ 0   0   0  -1   1   0   0   0   0   0   1   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0];
    % pK2 = (ph + pco3 - phco3);      % K2 = [CO3]pH]/[HCO3]
    ec3 =  [ 0   0   0   0  -1   1   0   0   0   0   1   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0];
    % pKb = (ph + pboh4 - pboh3);     % Kb = [BOH4][H]/[BOH3]
    ec4 = [ 0   0   0   0   0   0  -1   1   0   0   1   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0];
    % pKw = (ph + poh);               % Kw = [OH][H]
    ec5 = [ 0   0   0   0   0   0   0   0   0   1   1   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0];
    % pKs = (phf + pso4 - phso4);     % Ks = [hf][SO4]/[HSO4]
    ec6 = [ 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 ...
           +1   0   0  +1  -1   0   0   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0];
    % pKf = (ph + pf - pHF);          % Kf = [H][F]/[HF]
    ec7 =  [ 0   0   0   0   0   0   0   0   0   0   1   0   0   0   0 ...
           0  -1  +1   0   0   0   0   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0];
    % pK1p= (ph + ph2po4 - ph3po4);   % K1p = [H][H2PO4]/[H3PO4]
    ec8 = [ 0   0   0   0   0   0   0   0   0   0  +1   0   0   0   0 ...
           0   0   0   0   0   0   0  +1  -1   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0];
    % pK2p= (ph + phpo4 - ph2po4);    % K2p = [H][HPO4]/[H2PO4]
    ec9 = [ 0   0   0   0   0   0   0   0   0   0  +1   0   0   0   0 ...
           0   0   0   0   0   0  +1  -1   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0];
    % pK3p= (ph + ppo4 - phpo4);      % K3p = [H][PO4]/[HPO4]
    ec10 = [ 0   0   0   0   0   0   0   0   0   0  +1   0   0   0   0 ...
           0   0   0   0   0  +1  -1   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0   0   0];
    % pKSi= (ph + psiooh3 - psioh4);  % KSi = [H][SiOOH3]/[SiOH4]
    ec11 = [ 0   0   0   0   0   0   0   0   0   0  +1   0   0   0   0 ...
           0   0   0   0   0   0   0   0   0  +1  -1 ...
           0   0   0   0   0   0   0   0   0   0   0];   
    K = zeros(11,37); 
    K(1,:) = ec1;
    K(2,:) = ec2; 
    K(3,:) = ec3;
    K(4,:) = ec4;
    K(5,:) = ec5;
    K(5,:) = ec6; 
    K(7,:) = ec7;
    K(8,:) = ec8;    
    K(9,:) = ec9;
    K(10,:) = ec10; 
    K(11,:) = ec11;   
    %K = [ec1; ec2; ec3; ec4; ec5; ec6;...
           %ec7; ec8; ec9; ec10; ec11];


    %[ K0, K1, K2, p2f ] = local_equic(TC,S,P); % from FP's eg.c
    %pk1 = p( [ K0 / p2f; K1; K2 ] ); % from FP's eg.m
    pk = z(27:37); % already calculated from local_K in QUODcarb (I think) -MF

    % Make a vector of measured quantities
    i = find(~isnan(y));
    y = y(i);
    
    % Make a precision matrix
    W = diag(w(i));
    
    % Build a matrix that Picks out the measured components of x
    I = eye(37); % for chain rule
    P = I(i,:); % picking/pick out the measured ones
    
    c = [    M * q( x )  ;  
             pk - ( K * x )  ] ; % constraint equations
    
    e = P*x - y;               % error
    
    f = 0.5 *  e.' * W * e  + lam.' * c ;  % limp, method of lagrange multipliers
    % -(-1/2 sum of squares) + constraint eqns, minimize f => grad(f) = 0

    if ( nargout > 1 ) % compute the gradient
        
        dcdx = [ M * diag( dqdx( x ) ); ...
                     -K ]; % constraint eqns wrt p(concentrations)
        
        g = [ e.' * W * P +  lam.' * dcdx ,  c.' ];
        
    end
    if ( nargout > 2 ) % compute the Hessian
        
        ddq =  diag( d2qdx2( x ) );
        gg1 = lam(1)*diag(M(1,:))*ddq;
        gg2 = lam(2)*diag(M(2,:))*ddq;
        gg3 = lam(3)*diag(M(3,:))*ddq;
        gg4 = lam(4)*diag(M(4,:))*ddq;
        gg5 = lam(5)*diag(M(5,:))*ddq;
        gg6 = lam(6)*diag(M(6,:))*ddq;
        gg7 = lam(7)*diag(M(7,:))*ddq;
        gg8 = lam(8)*diag(M(8,:))*ddq;
        gg9 = lam(9)*diag(M(9,:))*ddq;
        gg10 = lam(10)*diag(M(10,:))*ddq;
        gg11 = lam(11)*diag(M(11,:))*ddq;
        gg12 = lam(12)*diag(M(12,:))*ddq;
        gg13 = lam(13)*diag(M(13,:))*ddq;
        gg14 = lam(14)*diag(M(14,:))*ddq;
        gg15 = lam(15)*diag(M(15,:))*ddq;
        gg16 = lam(16)*diag(M(16,:))*ddq;
        gg17 = lam(17)*diag(M(17,:))*ddq;
        gg18 = lam(18)*diag(M(18,:))*ddq;
        gg19 = lam(19)*diag(M(19,:))*ddq;
        gg20 = lam(20)*diag(M(20,:))*ddq;
        gg21 = lam(21)*diag(M(21,:))*ddq;
        gg22 = lam(22)*diag(M(22,:))*ddq;
        gg23 = lam(23)*diag(M(23,:))*ddq;
        gg24 = lam(24)*diag(M(24,:))*ddq;
        gg25 = lam(25)*diag(M(25,:))*ddq;
        gg26 = lam(26)*diag(M(26,:))*ddq;
        H = [  P.'*W*P + ...
            gg1 + gg2 + gg3 + gg4 + gg5 + gg6 + gg7 + ...
            gg8 + gg9 + gg10 + gg11 + gg12 + gg13 + ...
            gg14 + gg15 + gg16 + gg17 + gg18 + gg19 + ...
            gg20 + gg21 + gg22 + gg23 + gg24 + gg25 + gg26, ...
            dcdx.' ; ...
            dcdx , zeros(18) ];     
    end

end