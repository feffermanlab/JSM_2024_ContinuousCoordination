function Sol = VarCoordSolve(u0, h, T, tau,Kernel, rhoprime, r)
    Omega = h*(1:length(u0));
    EndTime = round(T/tau);
    Sol = zeros(length(u0),EndTime+1);
    Sol(:,1) = u0;
    for t =1:EndTime
        dt = zeros(length(Omega),1);
        for i = 1:length(Omega)
            x = Omega(i);
            K = Kernel(Omega-x);
            rho = rhoprime(Sol(:,t),i,r);
            integrand = K.*rho' ;
            c=trapz(integrand);
            crhr = sum(integrand(1:end-1).*h);
            dt(i)=crhr;
        end
        Sol(:,t+1)= Sol(:,t)+dt*tau;
    end
end