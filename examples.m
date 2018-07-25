% This MATLAB program checks the feasibility of LMIs from Theorems 1-3 of the paper 
% A. Selivanov and E. Fridman, "Sampled-data H-infinity filtering of a 2D heat equation under pointwise measurements," in 57th Conference on Decision and Control, 2018. 
D=[1 0; 0 .8]; a=2*pi^2;    % system parameters 
N=36;                       % number of sensors
epsilon=.05;                % parameter of (4) 
l=1/(2*sqrt(N))+epsilon/2;  % subdomain size (13) 
L=5;                        % observer gain 
alpha=.01;                  % decay rate 
%% Continuous-time observer 
if LMI_CDC18_th1(D,a,L,l,alpha)
    disp('Theorem 1: feasible'); 
else
    disp('Theorem 1: not feasible'); 
end
%% H-infty filtering under continuous measurements 
gamma=2.4; % L2-gain 
if LMI_CDC18_th2(D,a,L,l,alpha,gamma)
    disp('Theorem 2: feasible'); 
else
    disp('Theorem 2: not feasible'); 
end
%% H-infty filtering under sampled in time measurements 
gamma=4.6;          % L2-gain 
cmax=1/epsilon^2;   % = max||c_i||_\infty
h=0.001;            % sampling period 
if LMI_CDC18_th3(D,a,L,N,cmax,l,h,alpha,gamma)
    disp('Theorem 3: feasible'); 
else
    disp('Theorem 3: not feasible'); 
end