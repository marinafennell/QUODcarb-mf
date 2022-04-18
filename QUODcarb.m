% QUODcarb_v2
% starting 04/2022

function [yobs,sigobs,sigx] = QUODcarb(meas,Win,TC,S,P)

% z = [ pCT; pAT; pBT; pco2; phco3; pco3; pboh3; pboh4; ...
%       pfco2; poh; ph; pTS; pTF; pTP; pTSi; ...
%       phf; pHF; pf; pso4; phso4; ppo4; phpo4; ...
%       ph2po4; ph3po4; psiooh3; psioh4; ...
%       pKh; pK1; pK2; pKb; pKw; pKs; ...
%       pKf; pK1p; pK2p; pK3p; pKsi; ...
%       lam1; lam2; lam3; lam4; lam5; lam6; lam7; lam8; ...
%       lam9; lam10; lam11; lam12; lam13; lam14; lam15; ...
%       lam16; lam17; lam18; lam19; lam20; lam21; lam22; ...
%       lam23; lam24; lam25; lam26];
% z goes into limpcarb

%----------------------------------------------------------

p = @(x) -log10(x); % use p as in "pH" to denote the negative log base 10 function
q = @(x) 10.^(-x);  % use q, i.e., a backward p, for the inverse of the log10 function

y = p(meas); % -log10
w = p(Win); % -log10

% set lambdas = 0 to enforce equations
lam1 = 0; lam2 = 0; lam3 = 0; lam4 = 0;
lam5 = 0; lam6 = 0; lam7 = 0; lam8 = 0;
lam9 = 0; lam10 = 0; lam11 = 0; lam12 = 0;
lam13 = 0; lam14 = 0; lam15 = 0; lam16 = 0;
lam17 = 0; lam18 = 0; lam19 = 0; lam20 = 0;
lam21 = 0; lam22 = 0; lam23 = 0; lam24 = 0;
lam25 = 0; lam26 = 0;

% create an initial iterate to plug into 'z'
load 'init_c.mat' init_c ;
[Kh,K1,K2,Kb,Kw,Ks,Kf,K1p,K2p,K3p,Ksi,p2f] = local_K(TC,S,P);

pcsys = p(init_c);

Ksys = [Kh;K1;K2;Kb;Kw;Ks;Kf;K1p;K2p;K3p;Ksi];
pKsys(:) = p(Ksys);

z0(:) = [pcsys(:);pKsys(:);...
        lam1;lam2;lam3;lam4;lam5;lam6;lam7;lam8;lam9;...
        lam10;lam11;lam12;lam13;lam14;lam15;lam16;lam17;...
        lam18;lam19;lam20;lam21;lam22;lam23;lam24;lam25;lam26];


yobs = y;
wobs = w;
tol = 1e-9;
gun = @(z) glimpcarb(z,yobs,wobs,TC,S,P);
[z_dah,J,iflag] = newtn(z0,gun,tol);
if (iflag ~=0)
    fprintf('Newton''s method iflag = %i\n',iflag);
end
x = z_dah;

% measurement uncertainty
sigobs = sqrt(1./wobs);
u = q(yobs-sigobs);
l = q(yobs+sigobs);
sigobs(1:37) = 0.5*(u(1:37)-l(1:37));

% posterior uncertainty
C = inv(J);
C = C(1:37,1:37);
x = x(1:37);
sigx = sqrt(diag(C));
u = q(x-sigx);
l = q(x+sigx);
sigx(1:37) = 0.5*(u(1:37)-l(1:37));

end

%-----------------------------------------------------------------

function [g,H] = glimpcarb(z,y,w,TC,S,P)
    [~,g,H] = limpcarb(z,y,w,TC,S,P);
    g = g(:);
end





