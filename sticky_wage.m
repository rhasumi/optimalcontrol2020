%------------------------------------
% 硬直賃金モデル（sticky_wage.m）
%------------------------------------

%cd('/Users/rhasumi/OneDrive - Nezu Foundation of the Musashi Academy/work/note/continuous/search')
%addpath('../gensys');

clear all;

syms c k l mu a cdot kdot ldot mudot adot etac etamu epsilon
syms css kss lss muss ass
syms rho theta eta alpha delta zeta
syms w piw wdot piwdot etapiw nu gamma al

Yt = [c; l; mu; k; a; w; piw];
dotYt = [cdot; ldot; mudot; kdot; adot; wdot; piwdot];
Eta = [etamu; etapiw];
Epsilon = [epsilon];

Eq1 = c^(-theta)-mu;
Eq2 = -w + (1-alpha)*exp(a)*(k/l)^alpha;
Eq3 = (mudot-etamu)+mu*(alpha*exp(a)*(l/k)^(1-alpha)-delta-rho);
Eq4 = -kdot+exp(a)*k^alpha*l^(1-alpha)-delta*k-c-0.5*gamma*piw^2;
Eq5 = -adot-zeta*a+epsilon;
Eq6 = (1.0-nu)*w+nu*al/mu-rho*gamma*piw+theta*(mudot-etamu)/mu*piw+theta*(piwdot-etapiw);
Eq7 = -piw + wdot/w;

System_of_Eq = [Eq1; Eq2; Eq3; Eq4; Eq5; Eq6; Eq7];

neq  = length(System_of_Eq);
CC = zeros(neq,1);
GAM0j = -jacobian(System_of_Eq, dotYt);
GAM1j = jacobian(System_of_Eq, Yt);
PSI0j = jacobian(System_of_Eq, Epsilon);
PPIj = jacobian(System_of_Eq, Eta);

rho = 0.05;
theta = 1.0;
alpha = 0.4;
delta = 0.1;
zeta = 0.1;
eta = 1.0;

nu = 10.0;
gamma = 100.;
%al 

lss = 1.0
kss = ((delta+rho)/alpha)^(-1/(1-alpha));
css = kss^alpha-delta*kss;
%lss =(c_l^(-theta)*(1-alpha)*(k_l)^alpha)^(1/(eta+theta));
%css = c_l*lss;
%kss = k_l*lss;
muss = css^(-theta);
wss = (1-alpha)*kss^alpha;
al = (nu-1.0)/nu*muss*wss;

c = css;
l = lss;
mu = muss;
k = kss;
a = 0.;
cdot = 0.;
ldot = 0.;
mudot =0.;
kdot = 0.;
adot = 0.;
etamu = 0.;
epsilon = 0.;

w = wss;
piw = 0.;
wdot = 0.;
piwdot = 0.;
etapiw = 0.;

check = eval(System_of_Eq)

GAM0 = eval(GAM0j);
GAM1 = eval(GAM1j);
PSI0 = eval(PSI0j);
PPI  = eval(PPIj);

[G1,C,impact,q,aa,b,z,RC] = gensysct(GAM0,GAM1,CC,PSI0,PPI);
%[G1,C,impact,RC,F] = schur_solver(GAM0,GAM1,CC,PSI0,PPI,1,1,1);
RC

dt = 0.1;
nvar = size(Yt, 1);
nirf =1000;

yyirf  = zeros(nvar,nirf);
yyirf(:,1) = impact*[1];
for t = 2:nirf
        yyirf(:,t) = yyirf(:,t-1)+dt.*G1*yyirf(:,t-1);
end

figure(1)
subplot(3, 2, 1)
plot((1:nirf)*dt, yyirf(1,:),'b')
title('c(t)')

subplot(3, 2, 2)
plot((1:nirf)*dt, yyirf(2,:),'b')
title('n(t)')

subplot(3, 2, 3)
plot((1:nirf)*dt, yyirf(6,:),'b')
title('w(t)')

subplot(3, 2, 4)
plot((1:nirf)*dt, yyirf(7,:),'b')
title('\pi_w(t)')

subplot(3, 2, 5)
plot((1:nirf)*dt, yyirf(4,:),'b')
title('k(t)')

subplot(3, 2, 6)
plot((1:nirf)*dt, yyirf(5,:),'b')
title('a(t)')

saveas(1, 'stw.eps')
