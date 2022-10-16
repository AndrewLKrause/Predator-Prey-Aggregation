
% Number of gridpoints and domain size.
N=1000; L=1;

% Parameters from the 2020 paper.
delta1=0.002; delta2 = 0.002;
mu1 = 0.02; mu2=0.02;
eps1=0.0002; eps2=0.0002; 
a11=1; a22=1; a12=1; a21=-1;A1=1; A2=0;
b11=1; b22=1; b12=1; b21=2;

% Kinetic Fecundities (will be multiplied as uf, vg in model)
f = @(x,t,u,v)A1-a11*u-a12*v;g = @(x,t,u,v)A2-a21*u-a22*v;
uss = 1/2; vss = 1/2;
d = @(x)1+0*x;
sigmau = eps1; sigmav = eps2;

% Note that du here is really du(u,v)*Du(x,t) in the paper
du = @(x,t,u,v)delta1*d(x);
dv = @(x,t,u,v)delta2*d(x);

% Generalized Fecundity functions - gradients of these will appear.
F = @(x,t,u,v)mu1*(-b11*u-b12*v);G = @(x,t,u,v)mu2*(-b21*u-b22*v);

% Time span to simulate - from 0 to 200, with 1e4 resolved steps (NOT
% timesteps!)
T = linspace(0,200,1e4);

% Run the simulation with these parameters
[u, v,T] = run1DSim(N,L,f,g,F,G,d,du,dv,sigmau,sigmav,uss,vss,T);

