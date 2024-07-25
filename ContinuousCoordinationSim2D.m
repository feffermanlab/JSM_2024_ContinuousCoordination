%Domain Dimensions
X1=1;
X2=1;

%Spatial Mesh Size
N = 50;
h = 1/(N-1); %There are N gridpoints including the endpoints which are each 
             %h apart
%Max time  
T=10;

%Temporal Mesh Size
TN = 20;
tau = 1/TN; %There are TN Temporal Gridpoints not including t=0 which are 
            %each tau apart

%Initial conditions
Omegax1 = X1*h.*(0:N-1);
Omegax2 = X2*h.*(0:N-1);
[OmegaX1,OmegaX2]=meshgrid(Omegax1,Omegax2);
%Omega=meshgrid(Omegax1,Omegax2)

%u0 = @(x1,x2) 1./(1+exp(-80.*(x1-0.5)));
%u0=@(x1,x2) cos(2*pi.*x1).*sin(2*pi.*x2);
u0=@(x1,x2) cos(2*pi.*x1).*(x2.^2);

U0 = u0(OmegaX1,OmegaX2);
figure
mesh(OmegaX1,OmegaX2,U0)

%Recognition parameter
r=1;
%Familiarity parameter
s = 0.1;

%Plot limits:
x1lims = [0,X1*h*(N-1)];
x2lims = [0,X2*h*(N-1)];
zlims = [min(U0,[],'all'),max(U0,[],'all')];

figure
title('animation')
u=U0;
for i = 1:1000
    %compute derivative wrt time
    dt = zeros(size(U0));
    for j = (1:length(OmegaX1))
        for k = (1:length(OmegaX2))
            x1 = OmegaX1(k,j);
            j;
            x2 = OmegaX2(k,j);
            k;
            u;
            K=Kernel(sqrt((x1-OmegaX1).^2+(x2-OmegaX2).^2),s);
            D = DifIntSurfaceCPT(u,j,k,r);
            %c = sum(K.*D.*h^2,'all');
            %intfirst = trapz(K.*D);
            %intsecond = trapz(intfirst)*h^2;
            c = trapz(trapz(K.*D))*h^2;
            dt(k,j)=c;
        end
    end
    dt;
    u = u+(tau.*dt);
    mesh(OmegaX1,OmegaX2,u)
    xlabel('x1');
    ylabel('x2');
    zlabel('u(x1,x2)');
    xlim(x1lims);
    ylim(x2lims);
    zlim(zlims);
    drawnow
end


%function K = Kernel(x)
%    idx = (x==0);
%    K = (1-idx).*(1./x.^4);
%    K(isnan(K))=2^16;
%end

function K = Kernel(x,s)
    K = 1/(s*sqrt(2*pi))*exp((x).^2./(-2*s^2));
end

function R = DifIntSurface(u,j,k,r)
    d = u(k,j)-u;
    R = -2/r.*d.*exp(-1/r.*d.^2);
end


function R = DifIntSurfaceCPT(u,j,k,r)
    d = u(k,j)-u;
    idx = abs(d)<0.25;
    R = (idx).*-2/r.*d.*exp(-1/r.*d.^2);
end
