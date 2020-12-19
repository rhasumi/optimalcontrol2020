%--------------------
% SIRモデル（sir.m）
%--------------------

N = 1000;
T = 50.;
dt = T/N;
S = zeros(N,1);
I = zeros(N,1);
I(1) = 0.01;
S(1) = 1.0-I(1);
beta = 0.5;
nu = 0.2;

for t = 1:N-1
  S(t+1)=S(t)-beta*S(t)*I(t)*dt;
  I(t+1)=I(t)+beta*S(t)*I(t)*dt-nu*I(t)*dt;
end

plot(0:dt:T-dt, S)
hold on
plot(0:dt:T-dt, I,'--')
legend('S(t)','I(t)','FontSize',18)

saveas(1, 'sir.eps')
