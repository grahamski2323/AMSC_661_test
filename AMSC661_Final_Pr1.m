%% AMSC 661, Final Problem 1
% Graham Antoszewski
% 17 May 2016
% Objective: Set up BVP problem of stationary heat distribution
% plot results, making sure you have the correct BCs

%% set up mesh
nx = 300;
ny = 150;
xaux = linspace(0,1,nx);
yaux = linspace(0,1,ny);
[x,y] = meshgrid(xaux,yaux);

%% set up A matrix
hx = 1/nx;
hy = 1/ny;
% dimension of A
nxy = nx*ny;
e = ones(nxy,1);
T = spdiags([e -4*e e],[-1 0 1],ny,ny);
T(1,2) = 2; % for neumann BCs
I = speye(nx);
A1 = kron(I,T);
S = spdiags([e,e],[-ny,ny],nxy,nxy);
A = (A1 + S);
% taking into account the periodic boundary conditions
A1 = spdiags(e,-nxy + ny,nxy,nxy);
A2 = spdiags(e,nxy - ny,nxy,nxy);
A = A + A1 + A2;
A = (1/hx)*(1/hy)*A;
% creating right-hand side due to sun
rhs = zeros(nxy,1);
rhs(1:floor(nxy/2),1) = -10; % from the sun, negative since it is (u_xx+u_yy) + f/cp = 0

%% solving the system
uax = A\rhs;
u = zeros(ny,nx);
u(1:ny,1:nx) = reshape(uax,ny,nx);

%% contour plot
 ma = max(max(u));
 mi = min(min(u));

figure;
% the y-indexes of my matrix are top to bottom,
% but my y-indexes of the graph are bottom to top,
% so I have to flip my contour plot to accurately
% show what u is and have boundary conditions on the
% correct sides
contourf(u,linspace(mi,ma,10));
set(gca,'ydir','reverse')
colorbar;
set(gca,'DataAspectRatio',[1,1,1],'Fontsize',14);
xlabel('x','Fontsize',14);
ylabel('y','Fontsize',14);
title('Heat Distribution, Sun Shining on Left Half')

