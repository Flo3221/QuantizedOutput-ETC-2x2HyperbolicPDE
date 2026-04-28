% parameters
a = @(x) ones(size(x));
b = @(x) 2*ones(size(x));
c1 = @(x) 2*ones(size(x));
c2 = @(x) 2*ones(size(x));
p = 0.6;
r = 0.9;

xg = 128; % grid for x (minus one end point determined by boundary conditions)
tmp = linspace(0,1,xg+1); % auxiliary grid
x0 = tmp(1:xg); x1 = tmp(2:xg+1); % grids for v and u

% control and output injection gains
N = 13; % approximation order of power series
[K1, K2, ec] = kernelsc(a(1),b(1),c1(1),c2(1),p,r,x1,x0,N);
[P1, P2, eo] = kernelso(a(1),b(1),c1(1),c2(1),p,x1,x0,N);

% boundary conditions u(0)=pv(0) and v(1)=ru(1)+Ud
% -> 0 is a ghost point for u and 1 is a ghost point for v
% D is backward difference, -D' is forward difference
D = eye(xg) - diag(ones(1,xg-1), -1);

% initializations for finite-difference approximation
A = zeros(2*xg); % system operator 
B = zeros(2*xg, 1); % control operator
C = B'; % obseration operator

% control operator v(1) = U and observation y = u(1)
B(2*xg) = b(1)*xg;
C(xg) = 1;

% fill A
Ik = 1:xg; Il = xg+1:2*xg; % index sets
A(Ik,Ik) = -xg*a(x1).*D; % lambda tranport term
A(Il,Il) = -xg*b(x0).*D'; % mu transport term
A(Ik,Il) = diag(c1(x0)); % c1 term
A(Il,Ik) = diag(c2(x1)); % c2 term
A(Ik(1),Il(1)) = A(Ik(1),Il(1)) + xg*a(0)*p; % u(0)=pv(0)
% A for observer (same as A but v(1)=0 instead of v(1)=rv(0))
Ao = A - [P1; P2]*C; % add observer output injection
A(Il(xg),Ik(xg)) = A(Il(xg),Ik(xg)) + xg*b(1)*r; % v(1)=ru(1)
% observer control operator
Bo = [[P1; P2] + r*B, B];

% initial conditions
v0 = @(x) 10*(1-x);
z0 = [p*v0(x1), v0(x0)]';
z0hat = 1.5*z0;
 mu0 = 8;
 Omega=.6;

% run simulink model
% out = sim('q2x2sim.slx');