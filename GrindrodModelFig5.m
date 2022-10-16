
N=1000; L=1;

delta1=0.002; delta2 = 0.002;
mu1 = 0.02; mu2=0.02;
eps1=0.0002; eps2=0.0002; 
a11=1; a22=1; a12=1; a21=-1;A1=1; A2=0;
b11=1; b22=1; b12=1; b21=2;

f = @(x,t,u,v)A1-a11*u-a12*v;g = @(x,t,u,v)A2-a21*u-a22*v;
uss = 1/2; vss = 1/2;
d = @(x)1+0*x;
sigmau = eps1; sigmav = eps2;
du = @(x,t,u,v)delta1*d(x);
dv = @(x,t,u,v)delta2*d(x);
gradF = @(x,t,u,v)mu1*(-b11*u-b12*v);gradG = @(x,t,u,v)mu2*(-b21*u-b22*v);

T = linspace(0,200,1e4);

[u, v,T] = run1DSim(N,L,f,g,gradF,gradG,d,du,dv,sigmau,sigmav,uss,vss,T);

