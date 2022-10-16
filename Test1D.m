
N=200; L=10;

f = @(x,t,u,v)1-u-0.5*v;g = @(x,t,u,v)1-v-0.5*u;
uss = 2/3; vss = 2/3;
d = @(x)1+0*x;
sigmau = 0; sigmav = 0;
du = @(x,t,u,v)d(x);
dv = @(x,t,u,v)10*d(x);
gradF = f;gradG = g;

T = linspace(0,100,1e3);

[u, v,T] = run1DSim(N,L,f,g,gradF,gradG,d,du,dv,sigmau,sigmav,uss,vss,T);

