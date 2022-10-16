% This code simulates the system defined by the arguments. The time range
% *for plotting purposes only* is defined below and is given by T.
function [u, v,T] = run1DSim(N,L,f,g,F,G,d,du,dv,sigmau,sigmav,uss,vss,T)

% Spatial step and diffusion prefactor.
dx = L/(N-1);Dx = (1/dx^2);

% Evaluate a predefined grid x.
x = linspace(0,L,N)';

% Randomly perturb the given steady state.
rng(1);
uss=uss.*abs(1+1e-1*randn(N,1));
vss=uss.*abs(1+1e-1*randn(N,1));
uinit = [uss;vss];

ui = 1:N; vi=N+1:2*N;

% Form the Laplacian
%e = ones(N,1); % Vector of ones to use across the diagonals
%Lap= spdiags([e -2*e e], -1:1, N, N); % Diagonal Laplacian
LapH= spdiags([d(x) -2*d(x) d(x)], -1:1, N, N); % Diagonal weighted Laplacian
LapH(1,1) = -d(1); LapH(N,N) = -d(end); % Neumann boundary conditions
LapH = -LapH./dx^2; % Scale the finite-difference operator
I = speye(N);

RHSu = @(t,u,v)u.*f(x,t,u,v)+Dx*Lap(du(x,t,u,v),u)-Dx*Div(((sigmau*LapH+I)\Grad(F(x,t,u,v))).*d(x).*u);
RHSv = @(t,u,v)v.*g(x,t,u,v)+Dx*Lap(du(x,t,u,v),v)-Dx*Div(((sigmav*LapH+I)\Grad(G(x,t,u,v))).*d(x).*v);

RHS = @(t,U)[RHSu(t,U(ui),U(vi)); RHSv(t,U(ui),U(vi))];

% Unfortunately, as written the system is full and not sparse. This will
% make computation slow...
odeparams = odeset('RelTol',1e-11,'AbsTol',1e-11,'MaxStep',0.02);

[T,U] = ode15s(RHS,T,uinit,odeparams);

u = U(:,ui); v = U(:,vi);

    function Div = Div(H)
        Div = [0; H(3:end)-H(1:end-2);0]/2;
    end
    function Grad = Grad(H)
        Grad = [0; H(3:end)-H(1:end-2);0]/2;
    end
% Version of (u_(i+1)-u_i)*(D_(i+1)+D_i)-(u_i-u_(i-1))*(D_i+D_(i-1))/2
    function Lap = Lap(D,u)
        
        Lap = [(u(2)-u(1)).*(D(2)+D(1));...,
            (u(3:end)-u(2:end-1)).*(D(3:end)+D(2:end-1))-(u(2:end-1)-u(1:end-2)).*(D(2:end-1)+D(1:end-2));...,
            -(u(end)-u(end-1)).*(D(end)+D(end-1))]/2;
        % The factor of 1/2 comes from averaging the D_i in half-steps
    end

end