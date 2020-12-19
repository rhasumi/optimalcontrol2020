%------------------------------------
% ニューケインジアン・モデル（nk.m）
%------------------------------------

clear all;
%addpath('../gensys');

syms c i ppi m mu a v etac etapi etamu epsilon e
syms cdot idot pidot mdot mudot adot vdot
syms css iss piss mss muss ass vss 
syms rho theta eta eps gamma phipi phiy zeta xi

Yt = [c; i; ppi; m; mu; a; v];
dotYt = [cdot; idot; pidot; mdot; mudot; adot; vdot];
Eta = [etac; etapi; etamu];
Epsilon = [epsilon; e];

Eq1 = c^(-theta)-mu;
Eq2 = -(c/exp(a))^eta+mu*exp(a)*m;
Eq3 = (mudot-etamu)+mu*(i-ppi-rho);
Eq4 = -(i-ppi-(cdot-etac))*ppi+eps/gamma*(m-mss)+(pidot-etapi);
Eq5 = -i + iss+phipi*ppi+phiy*(c-exp(a)*css)+v;
Eq6 = -adot-zeta*a+epsilon;
Eq7 = -vdot-xi*v+e;

System_of_Eq = [Eq1; Eq2; Eq3; Eq4; Eq5; Eq6; Eq7];

neq  = length(System_of_Eq);
CC = zeros(neq,1);
GAM0j = -jacobian(System_of_Eq, dotYt);
GAM1j = jacobian(System_of_Eq, Yt);
PSI0j = jacobian(System_of_Eq, Epsilon);
PPIj = jacobian(System_of_Eq, Eta);

rho = 0.05;
theta = 1.0;
eta = 5.0;
eps = 5.0;
gamma = 2.0e+3;
phipi = 1.5;
phiy = 0.5;
zeta = 0.1;
xi = 0.1;

mss = (eps-1)/eps;
piss = 0;
iss = rho;
css = mss^(1/(eta+theta));
muss = css^(-theta);
ass = 0;
vss = 0;

c = css;
i = iss;
ppi = piss;
m = mss;
mu = muss;
a = ass;
v = vss;
cdot = 0;
idot = 0;
pidot = 0;
mdot = 0;
mudot = 0;
adot = 0;
vdot = 0;
etac = 0;
etapi = 0;
etamu = 0;
epsilon = 0;
e = 0;

check = eval(System_of_Eq)

GAM0 = eval(GAM0j);
GAM1 = eval(GAM1j);
PSI0 = eval(PSI0j);
PPI  = eval(PPIj);

[G1,C,impact,q,aa,b,z,RC] = gensysct(GAM0,GAM1,CC,PSI0,PPI);
%[G1,C,impact,RC,F] = schur_solver(GAM0,GAM1,CC,PSI0,PPI,1,1,1);
RC

dt = 0.1;
nj = 5;
ns = 2;
nvar = nj+ns;
nirf =1000;

% impulse res. to TFP shock

yyirf  = zeros(nvar,nirf);
yyirf(:,1) = impact*[1.; 0];
for t = 2:nirf
        yyirf(:,t) = yyirf(:,t-1)+dt.*G1*yyirf(:,t-1);
end

fs = 14;
figure(1)
subplot(2, 2, 1)
plot((1:nirf)*dt, yyirf(1,:),'b')
title('c(t)','FontSize',fs)

subplot(2, 2, 2)
plot((1:nirf)*dt, yyirf(2,:),'b')
title('i(t)','FontSize',fs)

subplot(2, 2, 3)
plot((1:nirf)*dt, yyirf(3,:),'b')
title('\pi(t)','FontSize',fs)

subplot(2, 2, 4)
plot((1:nirf)*dt, yyirf(6,:),'b')
title('a(t)','FontSize',fs)

saveas(1, 'nktfp.eps')

% impulse res. to MP shock

yyirf  = zeros(nvar,nirf);
yyirf(:,1) = impact*[0; 1];
for t = 2:nirf
        yyirf(:,t) = yyirf(:,t-1)+dt.*G1*yyirf(:,t-1);
end

fs = 14;
figure(2)
subplot(2, 2, 1)
plot((1:nirf)*dt, yyirf(1,:),'b')
title('c(t)','FontSize',fs)

subplot(2, 2, 2)
plot((1:nirf)*dt, yyirf(2,:),'b')
title('i(t)','FontSize',fs)

subplot(2, 2, 3)
plot((1:nirf)*dt, yyirf(3,:),'b')
title('\pi(t)','FontSize',fs)

subplot(2, 2, 4)
plot((1:nirf)*dt, yyirf(7,:),'b')
title('v(t)','FontSize',fs)

saveas(2, 'nkmp.eps')

