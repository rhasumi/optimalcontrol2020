%--------------------
% 相図（phase.m）
%--------------------

A = diag([-1 1]);

nt = 40;
a0 = [-1; -1];
avec = zeros(nt, 2);
dt = 0.2;
nj = 20;

a0 = [0; 0.01];
for i = 1:nt
  avec(i,:) = expm((i-1)*dt*A)*a0;
end
i = 1;
quiver(avec(i,1),avec(i, 2),avec(i+1,1)-avec(i,1),avec(i+1,2)-avec(i,2))
xlim([-1 1]); ylim([-1 1])
hold on;

for i = 1:nt-1
    p = quiver(avec(i,1),avec(i, 2),avec(i+1,1)-avec(i,1),avec(i+1,2)-avec(i,2));
    p.Color = 'black'; p.MaxHeadSize = 0.5;
    p = quiver(avec(i,1),-avec(i, 2),avec(i+1,1)-avec(i,1),-avec(i+1,2)+avec(i,2));
    p.Color = 'black'; p.MaxHeadSize = 0.5;
end

for j = 0:nj
  y = j/nj*2-1;
  a0 = [1; y];
  
  for i = 1:nt
    avec(i,:) = expm((i-1)*dt*A)*a0;
  end
  
  for i = 1:nt-1
    p = quiver(avec(i,1),avec(i, 2),avec(i+1,1)-avec(i,1),avec(i+1,2)-avec(i,2));
    p.Color = 'black'; p.MaxHeadSize = 0.5;
    p = quiver(-avec(i,1),avec(i, 2),-avec(i+1,1)+avec(i,1),avec(i+1,2)-avec(i,2));
    p.Color = 'black'; p.MaxHeadSize = 0.5;
  end
end

saveas(1, 'phase1.eps')

%--------------------

A = diag([-1 -0.5]);

nt = 80;
a0 = [-1; -1];
avec = zeros(2, nt);
dt = 0.2;
nj = 32;

figure(2)
for j = 1:nj
  x = 2*pi*j/nj;
  a0 = 3*[sin(x); cos(x)];
  
  for i = 1:nt
    avec(:,i) = expm((i-1)*dt*A)*a0;
  end
  
  for i = 1:nt-1
    p = quiver(avec(1,i),avec(2,i),avec(1,i+1)-avec(1,i),avec(2,i+1)-avec(2,i));
    hold on;
    p.Color = 'black'; p.MaxHeadSize = 0.5;
  end
  xlim([-1 1]); ylim([-1 1])
end

saveas(2, 'phase2.eps')
