function [t,y]=nh2cl(ph)
%Matlab code to verify NH2CL bacth model results of EPANET-MSX
%input argument: pH value 
%output: t = time in hour 
%        y = matrix to store the time series of the species concentration
%        y[:,j] is the time series vector of spcies with index j  
% species index in the DAE systemn
%HOCL   1       Hypochlorous acid
%NH3    2       Ammonia
%NH2CL  3       Monochloramine
%NHCL2  4       Dichloramine 
%I      5       Unidentified intermediate compound 
%H      6      
%ALK    7       Alkalinity     
%OCL    8       Hypochlorite Ion
%NH4    9       Ammonium Ion
%CO3    10      
%H2CO3  11
%HCO3   12
%OH     13

%M is the Mass matrix for DAE, first 7 equations are differential equations and the
%other 6 are algebraic equations
M = zeros(13,13);
for i = 1:7
    M(i,i)= 1.0;
end
% Initial condition
y0 =zeros(13,1);
y0(3) = 0.05e-3;        % Initial [NH2CL]
y0(6) = 10^(-ph);       % Initial [H]
y0(7) = 0.004;          % Initial [ALK]   
tspan = [0 168];        % Simulation Time
tol = 1e-8*ones(1,13);  % Absolute tolerance 

options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',tol, ...
                 'Vectorized','off', 'BDF', 'on');

[t,y] = ode15s(@f,tspan,y0,options);

%--------------------------------------------------------------------------


function out = f(t,y)

k1 = 1.5e10;
k2 = 7.6e-2;
k3 = 1.0e6;
k4 = 2.3e-3;
k6 = 2.2e8;
k7 = 4.0e5;
k8 = 1.0e8;
k9 = 3.0e7;
k10 = 55.0;

k5 = 2.5e7*y(6)+4.0e4*y(11)+800*y(12);
a1 = k1*y(1)*y(2);
a2 = k2*y(3);
a3 = k3*y(1)*y(3);
a4 = k4*y(4);
a5 = k5*y(3)*y(3);
a6 = k6*y(4)*y(2)*y(6);
a7 = k7*y(4)*y(13);
a8 = k8*y(5)*y(4);
a9 = k9*y(5)*y(3);
a10 = k10*y(3)*y(4);

%Evalutaion of DAEs 
out(1,1) = -a1+a2-a3+a4+a8;        %d[HOCL]/dt      
out(2,1) = -a1+a2+a5-a6;           %d[NH3]/dt
out(3,1) = a1-a2-a3+a4-a5+a6-a9-a10;    %d[NH2CL]/dt
out(4,1) = a3-a4+a5-a6-a7-a8-a10;       %d[NHCL2]/dt
out(5,1) = a7-a8-a9;                    %d[I]/dt
out(6,1) = 0;                           %constant pH
out(7,1) = 0;                           %constant alkalinity
%The following are equilibrium equations
out(8,1) = y(6)*y(8) - 3.16e-8*y(1);
out(9,1) = y(6)*y(2) - 5.01e-10*y(9);
out(10,1) = y(6)*y(10)- 5.01e-11*y(12);
out(11,1) = y(6)*y(12) - 5.01e-7*y(11);
out(12,1) = y(7) - y(12) - 2*y(10) - y(13) + y(6);
out(13,1) = y(6)*y(13) - 1.0e-14;
