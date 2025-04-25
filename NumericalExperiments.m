clear
%Length of Domain:
L=1;
%Spatial Mesh Size
N = 200;
h=1/(N-1); %There are N*L gridpoints which are each h apart
tau = 0.01;
%Line width for plots
lw =1.5;

%%CODE FOR NUMERICAL EXPERIMENT 1
Omega = linspace(-1/2,1/2,N);
u1=1/2*Omega;
u2=3/2*Omega;
u3=7/2*Omega;

Sol1=CoordSolve(u1,h,20,tau,@ker,@cptrhoprime);
Sol2=CoordSolve(u2,h,20,tau,@ker,@cptrhoprime);
Sol3=CoordSolve(u3,h,20,tau,@ker,@cptrhoprime);

lmin=0
lmax=4
lres=100

%SolNE1 = NumericalExperiment1(lmin,lmax,lres,@ker,@cptrhoprime,L,h,20,tau);
%writematrix(SolNE1,'DataFiles\NumericalExp1.csv');

tiledlayout(2,3)
nexttile([1,3])
SolNE1 = readmatrix("DataFiles\NumericalExp1.csv");
length(SolNE1(1,:))
labeling = NaN*ones(length(SolNE1(1,:)),1);
labeling(round ((1/2-lmin)/(lmax-lmin)*lres))=1/2;
labeling(round((3/2-lmin)/(lmax-lmin)*lres))=3/2;
labeling(round((7/2-lmin)/(lmax-lmin)*lres))=7/2;
length(labeling)
logbinNE1=log2(SolNE1/(L/h));
figure(1)
h1=heatmap(logbinNE1,'GridVisible','off');
h1.XDisplayLabels = labeling
h1.YDisplayLabels = NaN*ones(length(SolNE1(:,1)),1);
h1.XLabel = 'Initial Slope';
h1.YLabel = 'Range of u(\Omega,T)';
h1.Title = 'Coordination Bifurcations';
h1.Colormap = hot;
hs1 = struct(h1);
ylabel(hs1.Colorbar, "log_2 density of u(\Omega,T)");



%figure('Position',[10,10,300,300])
nexttile
h2=plot(Omega,Sol1(:,1),':',Omega,Sol1(:,end),'-')
ylim([-3.5/2,3.5/2]);
xlabel('x')
ylabel('u(x,t)')
title('u(x,0)=1/2x')
set(h2(1),'Color',"#0000a4","LineWidth",lw)
set(h2(2),'Color','#bc272d',"LineWidth",lw)
legend('u(\Omega,0)','u(\Omega,20)','location','southeast')

%figure('Position',[10,10,300,300])
nexttile
h3=plot(Omega,Sol2(:,1),':',Omega,Sol2(:,end),'-')
ylim([-3.5/2,3.5/2]);
xlabel('x')
ylabel('u(x,t)')
title('u(x,0)=3/2x')
set(h3(1),'Color',"#0000a4","LineWidth",lw)
set(h3(2),'Color','#bc272d',"LineWidth",lw)
legend('u(\Omega,0)','u(\Omega,20)','location','southeast')

%figure('Position',[10,10,300,300])
nexttile
h4=plot(Omega,Sol3(:,1),':',Omega,Sol3(:,end),'-')
ylim([-3.5/2,3.5/2]);
xlabel('x')
ylabel('u(x,t)')
title('u(x,0)=7/2x')
set(h4(1),'Color',"#0000a4","LineWidth",lw)
set(h4(2),'Color','#bc272d',"LineWidth",lw)
legend('u(\Omega,0)','u(\Omega,20)','location','southeast')


%%CODE FOR NUMERICALEXPERIMENT 2
Omega = linspace(-1/2,1/2,N);
u1 = 1./(1+exp(-2.*(Omega)));
u2 = 1./(1+exp(-4.*(Omega)));
u3 = 1./(1+exp(-14.*(Omega)));

Sol1=CoordSolve(u1,h,20,tau,@ker,@cptrhoprime);
Sol2=CoordSolve(u2,h,20,tau,@ker,@cptrhoprime);
Sol3=CoordSolve(u3,h,20,tau,@ker,@cptrhoprime);

lmin=0
lmax=15
lres=100

%SolNE2 = NumericalExperiment2(lmin,lmax,lres,@ker,@cptrhoprime,L,h,20,tau);
%writematrix(SolNE2,'DataFiles\NumericalExp2.csv');

tiledlayout(2,3)
nexttile([1,3])
SolNE2 = readmatrix("DataFiles\NumericalExp2.csv");
length(SolNE2(1,:))
labeling = NaN*ones(length(SolNE2(1,:)),1);
labeling(round ((2-lmin)/(lmax-lmin)*lres))=2;
labeling(round((4-lmin)/(lmax-lmin)*lres))=4;
labeling(round((14-lmin)/(lmax-lmin)*lres))=14;
length(labeling)
logbinNE2=log2(SolNE2/(L/h));
figure(1)
h1=heatmap(logbinNE2,'GridVisible','off');
h1.XDisplayLabels = labeling
h1.YDisplayLabels = NaN*ones(length(SolNE2(:,1)),1);
h1.XLabel = 'Sigmoid Parameter';
h1.YLabel = 'Range of u(\Omega,T)';
h1.Title = 'Coordination Bifurcations';
h1.Colormap = hot;
hs1 = struct(h1);
ylabel(hs1.Colorbar, "log_2 density of u(\Omega,T)");



%figure('Position',[10,10,300,300])
nexttile
h2=plot(Omega,Sol1(:,1),':',Omega,Sol1(:,end),'-')
ylim([0,1]);
xlabel('x')
ylabel('u(x,t)')
title('u(x,0)=(1+Exp(-2x))^{-1}')
set(h2(1),'Color',"#0000a4","LineWidth",lw)
set(h2(2),'Color','#bc272d',"LineWidth",lw)
legend('u(\Omega,0)','u(\Omega,20)','location','southeast')

%figure('Position',[10,10,300,300])
nexttile
h3=plot(Omega,Sol2(:,1),':',Omega,Sol2(:,end),'-')
ylim([0,1]);
xlabel('x')
ylabel('u(x,t)')
title('u(x,0)=(1+Exp(-4x))^{-1}')
set(h3(1),'Color',"#0000a4","LineWidth",lw)
set(h3(2),'Color','#bc272d',"LineWidth",lw)
legend('u(\Omega,0)','u(\Omega,20)','location','southeast')

%figure('Position',[10,10,300,300])
nexttile
h4=plot(Omega,Sol3(:,1),':',Omega,Sol3(:,end),'-')
ylim([0,1]);
xlabel('x')
ylabel('u(x,t)')
title('u(x,0)=(1+Exp(-14x))^{-1}')
set(h4(1),'Color',"#0000a4","LineWidth",lw)
set(h4(2),'Color','#bc272d',"LineWidth",lw)
legend('u(\Omega,0)','u(\Omega,20)','location','southeast')

%%CODE FOR NUMERICAL EXPERIMENT 3
rmin=0
rmax=0.7
rres=200

Omega = linspace(-1/2,1/2,N);
u = Omega;

Sol4=VarCoordSolve(u,h,20,tau,@ker,@VarCptrhoprime,0.04);
Sol5=VarCoordSolve(u,h,20,tau,@ker,@VarCptrhoprime,0.15);
Sol6=VarCoordSolve(u,h,20,tau,@ker,@VarCptrhoprime,0.4);

lw =1.5;

%SolNE3 = NumericalExperiment3(rmin,rmax,rres,@ker,@VarCptrhoprime,L,h,20,tau);
%writematrix(SolNE3,'DataFiles\NumericalExp3.csv');

tiledlayout(2,3)
nexttile([1,3])
SolNE3 = readmatrix("DataFiles\NumericalExp3.csv");
length(SolNE3(1,:))
labeling = NaN*ones(length(SolNE3(1,:)),1);
labeling(round ((0.04-rmin)/(rmax-rmin)*rres))=0.04;
labeling(round((0.15-rmin)/(rmax-rmin)*rres))=0.15;
labeling(round((0.4-rmin)/(rmax-rmin)*rres))=0.4;
length(labeling)
logbinNE3=log2(SolNE3/(L/h));
figure(1)
h1=heatmap(logbinNE3,'GridVisible','off');
h1.XDisplayLabels = labeling
h1.YDisplayLabels = NaN*ones(length(SolNE3(:,1)),1);
h1.XLabel = "Radius of Supp(\rho')";
h1.YLabel = 'Range of u(\Omega,T)';
h1.Title = 'Coordination Bifurcations';
h1.Colormap = hot;
hs1 = struct(h1);
ylabel(hs1.Colorbar, "log_2 density of u(\Omega,T)");



%figure('Position',[10,10,300,300])
nexttile
h2=plot(Omega,Sol4(:,1),':',Omega,Sol4(:,end),'-')
ylim([-1/2,1/2]);
xlabel('x')
ylabel('u(x,t)')
title('r=0.04')
set(h2(1),'Color',"#0000a4","LineWidth",lw)
set(h2(2),'Color','#bc272d',"LineWidth",lw)
legend('u(\Omega,0)','u(\Omega,20)','location','southeast')

%figure('Position',[10,10,300,300])
nexttile
h3=plot(Omega,Sol5(:,1),':',Omega,Sol5(:,end),'-')
ylim([-1/2,1/2]);
xlabel('x')
ylabel('u(x,t)')
title('r=0.15')
set(h3(1),'Color',"#0000a4","LineWidth",lw)
set(h3(2),'Color','#bc272d',"LineWidth",lw)
legend('u(\Omega,0)','u(\Omega,20)','location','southeast')

%figure('Position',[10,10,300,300])
nexttile
h4=plot(Omega,Sol6(:,1),':',Omega,Sol6(:,end),'-')
ylim([-1/2,1/2]);
xlabel('x')
ylabel('u(x,t)')
title('r=0.4')
set(h4(1),'Color',"#0000a4","LineWidth",lw)
set(h4(2),'Color','#bc272d',"LineWidth",lw)
legend('u(\Omega,0)','u(\Omega,20)','location','southeast')

%Numerical Experiment 4
Omega = linspace(-1/2,1/2,N);
u4=7/2*Omega;

Sol4=CoordSolve(u4,h,30,tau,@ker,@cptrhoprime);

tiledlayout(1,4)
nexttile
h2=plot(Omega,Sol4(:,1));
ylim([-1.75,1.75]);
xlabel('x')
ylabel('u(x,t)')
title('t=0')
set(h2(1),'Color',"#0000a4","LineWidth",lw)
legend('u(\Omega,0)','location','southeast')
nexttile
h2=plot(Omega,Sol4(:,round(length(Sol4(1,:))/3)));
ylim([-1.75,1.75]);
xlabel('x')
title('t=10')
set(h2(1),'Color',"#0000a4","LineWidth",lw)
legend('u(\Omega,10)','location','southeast')
nexttile
h2=plot(Omega,Sol4(:,round(2*length(Sol4(1,:))/3)));
ylim([-1.75,1.75]);
xlabel('x')
title('t=20')
set(h2(1),'Color',"#0000a4","LineWidth",lw)
legend('u(\Omega,20)','location','southeast')
nexttile
h2=plot(Omega,Sol4(:,end));
ylim([-1.75,1.75]);
xlabel('x')
title('t=30')
set(h2(1),'Color',"#0000a4","LineWidth",lw)
legend('u(\Omega,30)','location','southeast')



function res = NumericalExperiment1(lmin,lmax,lres,ker,rhoprime,L,h,T,tau)
%NumericalExperiment1 This function finds a solution with many linear
%initial data with the same kernel and recognition function
% lmin:     The smallest slope of the initial data
% lmax:     The largest slope of the initial data
% lres:     The number of different initial data to test
% ker:      The Kernel to be used in the IVP
% rhoprime: The recognition function to be used in the IVP
% L:        The length of the spatial domain
% h:        The size of the spatial mesh
% T:        The max time
% tau:      The size of the temporal mesh

%%This function produces a solution for each initial condition and reports
%%only the solution at the final time T in a matrix the size of lres and
%%the spatial grid. 
    lspace = linspace(lmin,lmax,lres);
    Omega = linspace(-1/2,1/2,round(L/h+1));
    Range = linspace(-lmax/2,lmax/2,round(L/h+1));
    result = zeros(round(L/h+1),lres);
    for i = 1:lres
        l = lspace(i);
        u=Omega*l;
        Sol = CoordSolve(u,h,T,tau,ker,rhoprime);
        final = Sol(:,end);
        for j = 1:length(Range)-1
            result(j,i) = sum(final>= Range(j)& final<Range(j+1));
        end
    end
    res = result
end

function res = NumericalExperiment2(lmin,lmax,lres,ker,rhoprime,L,h,T,tau)
%NumericalExperiment2 This function finds a solution with many sigmoid
%initial data with the same kernel and recognition function
% lmin:     The smallest sigmoid parameter of the initial data
% lmax:     The largest sigmoid parameter of the initial data
% lres:     The number of different initial data to test
% ker:      The Kernel to be used in the IVP
% rhoprime: The recognition function to be used in the IVP
% L:        The length of the spatial domain
% h:        The size of the spatial mesh
% T:        The max time
% tau:      The size of the temporal mesh

%%This function produces a solution for each initial condition and reports
%%only the solution at the final time T in a matrix the size of lres and
%%the spatial grid. 
    lspace = linspace(lmin,lmax,lres);
    Omega = linspace(-1/2,1/2,round(L/h+1));
    Range = linspace(0,1,round(L/h+1));
    result = zeros(round(L/h+1),lres);
    for i = 1:lres
        l = lspace(i);
        u = 1./(1+exp(-l.*(Omega)));
        Sol = CoordSolve(u,h,T,tau,ker,rhoprime);
        final = Sol(:,end);
        for j = 1:length(Range)-1
            result(length(Range)+1-j,i) = sum(final>= Range(j)& final<Range(j+1));
        end
    end
    res = result
end

function res = NumericalExperiment3(rmin,rmax,rres,ker,rhoprime,L,h,T,tau)
%NumericalExperiment3 This function finds a solution with a single linear
%intital condition and kernel and with many different recognition functions
% rmin:     The smallest radius of support for the recognition function 
% lmax:     The largest radius of support for the recognition function
% lres:     The number of different recognition functions to test
% ker:      The Kernel to be used in the IVP
% rhoprime: The recognition function to be used in the IVP (must contain
%               parameter r)
% L:        The length of the spatial domain
% h:        The size of the spatial mesh
% T:        The max time
% tau:      The size of the temporal mesh

%%This function produces a solution for each recognition function and reports
%%only the solution at the final time T in a matrix the size of rres and
%%the spatial grid. 
    rspace = linspace(rmin,rmax,rres);
    Omega = linspace(-1/2,1/2,round(L/h+1));
    u=Omega;
    result = zeros(round(L/h+1),rres);
    for i = 1:rres
        r=rspace(i);
        Sol = VarCoordSolve(u,h,T,tau,ker,rhoprime,r);
        final = Sol(:,end);
        for j= 1:length(Omega)-1
            result(length(Omega)+1-j,i) = sum(final>= Omega(j)& final<Omega(j+1));
        end
    end
    res=result
end

function K = ker(x)
%%KER the kernel function used throughout this numerical experiement
    s=0.5;
    K = 1/(s*sqrt(2*pi))*exp((x).^2./(-2*s^2));
end

function R = rhoprime(u,i)
%a recognition function with only 1 zero on the whole real line
    r=0.5;
    t = u(i);
    d = t-u;
    R = (-2/r).*d.*exp((-1/r).*d.^2);
end

function R= cptrhoprime(u,i)
%a recognition function with compact support
    r=0.2;
    t =u(i);
    d=t-u;
    bool = abs(d)<r;
    rhoprime = -2*(d/r)./((1-(d/r).^2).^2).*exp(-1./(1-(d/r).^2));
    R=bool.*rhoprime;
    R(isnan(R))=0;
end

function R= VarCptrhoprime(u,i,r)
% a recognition function with compact support determined by the parameter r
    t =u(i);
    d=t-u;
    bool = abs(d)<r;
    rhoprime = -2*(d/r)./((1-(d/r).^2).^2).*exp(-1./(1-(d/r).^2));
    R=bool.*rhoprime;
    R(isnan(R))=0;
end