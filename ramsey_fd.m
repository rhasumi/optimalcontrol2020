%------------------------------------
% ラムゼイモデル・HJB方程式
% （ramsey_model_fd.m）
%------------------------------------

rho = 0.05;
theta = 2.0;
alpha = 0.4;
delta = 0.1;
  
IA = 200;
amin = 0.1;
amax = 20.0;
a = linspace(amin,amax,IA)';
da = (amax-amin)/(IA-1);

maxit= 100;
crit = 10^(-6);
Delta = 100000;

v0 = zeros(IA);
Aa = zeros(IA, IA);
b = zeros(IA, 1);
B0 = (1/Delta + rho)*spdiags(ones(IA),0,IA,IA);
c = zeros(IA, 1);
ss = zeros(IA, 1);
v0 = (2.0*a).^(1-theta)./(1-theta)./rho;

V = v0;
V_old = v0;

dVf = zeros(IA,1);
dVb = zeros(IA,1);
dVc = zeros(IA,1);

for ii = 1:maxit
  disp(ii)
  V_old = V;
  
  dVf(1:IA-1) = (V(2:IA)-V(1:IA-1))/da;
  dVf(IA) = (amax^alpha -delta*amax)^(-theta);
  dVb(2:IA) = (V(2:IA)-V(1:IA-1))/da;
  dVb(1) = (amin^alpha -delta*amin)^(-theta);
  dVc(2:IA-1) = 0.5*(V(3:IA)-V(1:IA-2))/da;
  dVc(1) = dVf(1);
  dVc(IA) = dVb(IA);
  
  cf = dVf.^(-1.0/theta);
  ssf = a.^alpha -delta*a - cf;
  cb = dVb.^(-1.0/theta);
  ssb = a.^alpha -delta*a - cb;
  c0 = a.^alpha -delta*a;
  dV0 = c0.^(-theta);
  cc = dVc.^(-1.0/theta);
  ssc = a.^alpha -delta*a - cc;
  
  If = ssf > 0;
  Ib = ssb < 0;
  Iboth = (ssf > 0).* (ssb < 0);
  Ineither = (ssf <= 0).* (ssb >= 0);
  dV_Upwind = dVf.*(If-Iboth) + dVb.*Ib + dV0.*Ineither;
  c = dV_Upwind.^(-1.0/theta);
  u = (c.^(1.0-theta)-1.0)/(1.0-theta);
  
  Xa = -(Ib).*ssb./da;
  Ya = -(If-Iboth).*ssf./da+(Ib).*ssb./da;
  Za = (If-Iboth).*ssf./da;
  ss = ssf.*(If-Iboth) + ssb.*Ib;
  
  % A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
  Aa = spdiags(Ya(:,1),0,IA,IA)+spdiags(Xa(2:IA,1),-1,IA,IA)+spdiags([0;Za(1:IA-1,1)],1,IA,IA);
  b = u + V./Delta;
  V = (B0-Aa)\b;
  
  dist = max(max(abs(V_old-V)));
  if dist < 0.0000001
    disp('converged')
    break
  end
end

fs = 18;
plot(a, V); title('V(k)','FontSize',fs*1.2); xlabel('k(t)','FontSize',fs)
h = gca; h.XAxis.FontSize = fs*0.9; h.YAxis.FontSize = fs*0.9;
saveas(1, 'ramsey_v.eps')

plot(a, c); title('c(k)','FontSize',fs*1.2); xlabel('k(t)','FontSize',fs)
h = gca; h.XAxis.FontSize = fs*0.9; h.YAxis.FontSize = fs*0.9;
saveas(1, 'ramsey_c.eps')

plot(a, ss); title('s(k)','FontSize',fs*1.2); xlabel('k(t)','FontSize',fs)
h = gca; h.XAxis.FontSize = fs*0.9; h.YAxis.FontSize = fs*0.9;
saveas(1, 'ramsey_s.eps')

