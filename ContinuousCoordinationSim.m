%Length of Domain:
L=1;

%Spatial Mesh Size
N = 1000;
h=1/(N+1); %There are N gridpoints (excluding endpoints) which are each h apart

%Max time
T = 10;

%Tempotal Mesh Size
TN=100;%FIX THIS
tau = 1/TN; %There are N Temporal gridpoints not including t=0 which are each tau apart

%initial condition
Omega = cat(2,0,h.*(1:N),L);
u = 1./(1+exp(-40.*Omega+20));

%familiarity Parameter
s=0.005;
%Recognition Parameter
r=0.2;


figure 
title('animation')

G= plot(nan,nan);
axis([min(Omega) max(Omega) min(u) max(u)])
for j = 1:200
    %compute derivative wrt time
    dt = zeros(1,N+2);
    for i = 2:length(Omega)-1
        x=Omega(i);
        c=conv(Kernal(Omega-x,s),DifIntSurface(u,x,h,r).*h,'same');
        dt(i)=c(i);
    end
    u = u+(tau.*dt);
    G.XData = Omega;
    G.YData = u;
    drawnow
end

%x = Omega(10);
%plot(Omega,conv(Kernal(Omega-x,0.1),DifIntSurface(u,x,h).*h,'same'))
%hold on
%plot([x,x],[-1,1])
%hold off


function K = Kernal(x,s)
    K = 1/(s*sqrt(2*pi))*exp((x).^2./(-2*s^2));
end

function S = InteligabilitySurface(u,x,h,r)
    t = u(int32(x/h));
    d = t-u;
    S = exp((-1/r).*d.^2);
end

function R = DifIntSurface(u,x,h,r)
    t = u(int32(x/h));
    d = t-u;
    R = (-2/r).*d.*exp((-1/r).*d.^2);
end
