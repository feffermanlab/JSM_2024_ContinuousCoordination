function Sol = CoordSolve(u0, h, T, tau,Kernel, rhoprime )
%COORDSOLVE solves the Initial value problem
%u_t=int_\OmegaK(x,y)\rho'(u(x,t)-u(y,t))dy in one dimension 
%   u0:     Initial Data
%   h:      Size of the Spatial Mesh
%   T:      Max time
%   tau:    Size of the Time Mesh
%   Kernel: Familiarity Kernel, function on OmegaXOmega
%   rhoprime: Derivative of Recognition function on RR

%%This function solves the IVP with a first order forward difference
%%method and a simple righthand quadrature on the nonlocality. The output
%%is the solution which is a matrix the size of the space time mesh.  
    Omega = h*(1:length(u0));
    EndTime = round(T/tau);
    Sol = zeros(length(u0),EndTime+1);
    Sol(:,1) = u0;
    for t =1:EndTime
        dt = zeros(length(Omega),1);
        for i = 1:length(Omega)
            x = Omega(i);
            K = Kernel(Omega-x);
            rho = rhoprime(Sol(:,t),i);
            integrand = K.*rho';
            crhr = sum(integrand(1:end-1).*h);
            dt(i)=crhr;
        end
        Sol(:,t+1)= Sol(:,t)+dt*tau;
    end
end