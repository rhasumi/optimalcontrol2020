%------------------------------------
% ラムゼイモデル・相図
% （ramsey_model_pf.m）
%------------------------------------

N = 801;
T = 80;
dt = T/(N-1);
nt = (N-1)/T;
ct = zeros(N,1);
kt = zeros(N,1);

%cdot1 = (alpha*k^(alpha-1)-delta-rho)*c/theta;
%kdot1 = k^alpha-delta*k-c;

rho = 0.05;
theta = 2.0;
alpha = 0.4;
delta = 0.1;
kss = ((delta+rho)/alpha)^(1/(alpha-1));
css = kss^alpha-delta*kss;

k0a = kss*0.5;
k0b = kss*2.0;
ct(N)=css;
niter = 2000;

% SOR 法 (逐次過緩和法, successive overrelaxation method
% 緩和パラメータ

kt(1)=k0a;
kt(2:N) = k0a;
relax = 0.005;

ktold = kt;
kttemp = kt;

for t =  N-1:-1:1
  ct(t) = ct(t+1)-(alpha*kt(t)^(alpha-1)-delta-rho)*ct(t+1)/theta*dt;
end
for t = 1:N-1
  kttemp(t+1)=kttemp(t)+(kttemp(t)^alpha-delta*kttemp(t)-ct(t))*dt;
end
kt(1:N) = ktold(1:N)*(1-relax)+relax*kttemp(1:N);

for iter = 1:niter
  %iter
  try
    for t =  N-1:-1:1
      ct(t) = ct(t+1)-(alpha*kt(t)^(alpha-1)-delta-rho)*ct(t+1)/theta*dt;
    end
    for t = 1:N-1
      kttemp(t+1)=kttemp(t)+(kttemp(t)^alpha-delta*kttemp(t)-ct(t))*dt;
    end
    if sum(abs(imag(ct))+abs(imag(kttemp))) > 0
      throw(1)
    end
    if max(abs(kttemp-kt)) < 0.001
      'converged'
      break
    end
    %'ok'
    ktold = kt;
    kt(1:N) = ktold(1:N)*(1-relax)+relax*kttemp(1:N);
  catch
    'err'
    kt(1:N) = ktold(1:N)*0.5+0.5*kt(1:N);
    %continue
  end
end
iter

kta = kt;
cta = ct;

%--------------------

kt(1) = k0b;
kt(2:N) = kss;
relax = 0.0075;

ktold = kt;
kttemp = kt;

for t =  N-1:-1:1
  ct(t) = ct(t+1)-(alpha*kt(t)^(alpha-1)-delta-rho)*ct(t+1)/theta*dt;
end
for t = 1:N-1
  kttemp(t+1)=kttemp(t)+(kttemp(t)^alpha-delta*kttemp(t)-ct(t))*dt;
end
kt(1:N) = ktold(1:N)*(1-relax)+relax*kttemp(1:N);

for iter = 1:niter
  %iter
  try
    for t =  N-1:-1:1
      ct(t) = ct(t+1)-(alpha*kt(t)^(alpha-1)-delta-rho)*ct(t+1)/theta*dt;
    end
    for t = 1:N-1
      kttemp(t+1)=kttemp(t)+(kttemp(t)^alpha-delta*kttemp(t)-ct(t))*dt;
    end
    if sum(abs(imag(ct))+abs(imag(kttemp))) > 0
      throw(1)
    end
    if max(abs(kttemp-kt)) < 0.001
      'converged'
      break
    end
    %'ok'
    ktold = kt;
    kt(1:N) = ktold(1:N)*(1-relax)+relax*kttemp(1:N);
  catch
    'err'
    kt(1:N) = ktold(1:N)*0.5+0.5*kt(1:N);
    %continue
  end
end
iter

ktb = kt;
ctb = ct;

%------------------------

cdot = @(k, c, dt) (alpha*k^(alpha-1)-delta-rho)*c/theta*dt;
kdot = @(k, c, dt) (k^alpha-delta*k-c)*dt;

kset = 2:0.01:12.0;
plot(kset, kset.^alpha-delta*kset,'black')
hold on
xline(kss,'black')
for i = 1:T-1
  p = quiver(kta((i-1)*nt+1),cta((i-1)*nt+1),kdot(kta((i-1)*nt+1),cta((i-1)*nt+1),1.0),cdot(kta((i-1)*nt+1),cta((i-1)*nt+1),1.0));
  p.Color = 'blue'; p.MaxHeadSize = 0.75;
  p = quiver(ktb((i-1)*nt+1),ctb((i-1)*nt+1),kdot(ktb((i-1)*nt+1),ctb((i-1)*nt+1),1.0),cdot(ktb((i-1)*nt+1),ctb((i-1)*nt+1),1.0));
  p.Color = 'blue'; p.MaxHeadSize = 0.75;
end

kset2 = 2.0:1.0:11.5;
cset2 = 0.8:0.1:2.2;

for i = 1:size(kset2, 2)
  for j = 1:size(cset2, 2) 
    p = quiver(kset2(i),cset2(j),kdot(kset2(i),cset2(j),1.0),cdot(kset2(i),cset2(j),1.0));
    p.Color = 'red'; p.MaxHeadSize = 0.75;
  end
end
xlim([1.75, 11.25])
ylim([0.75, 2.25])
xlabel('k(t)'); ylabel('c(t)'); 
h = gca; h.XAxis.FontSize = 16; h.YAxis.FontSize = 16;

saveas(1, 'ramsey_phase.eps','epsc')
