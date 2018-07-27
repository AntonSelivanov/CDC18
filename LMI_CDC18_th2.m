function flag=LMI_CDC18_th2(D,a,L,l,alpha,gamma)
% This MATLAB program checks the feasibility of LMIs from Theorem 2 of the paper 
% A. Selivanov and E. Fridman, "Sampled-data H-infinity filtering of a 2D heat equation under pointwise measurements," in 57th Conference on Decision and Control, 2018. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% D, a  - parameters of (20)
% L     - observer gain from (5)
% l     - subdomain size (13)
% alpha - decay rate 
% gamma - L2-gain

% Output: 
% flag =1 if feasible, =0 otherwise
%% Decision variables and notations 
sdpvar p1 p2 p3 gamma1 eta lambda1 lambda2 lambda3 lambda4 lambda5 lambda6
gamma2=gamma^2*gamma1; 
P=[p1 p2; p2 p3]; 
pbar=[p1; 2*p2; p3]; 
dbar=[D(1,1); 2*D(1,2); D(2,2)]; 
%% LMIs 
Psi=blkvar; 
Psi(1,1)=2*(a-L+alpha)-(lambda5+lambda6)*pi^2+gamma1; 
Psi(1,3)=1; 
Psi(1,4)=1; 
Psi(1,5)=1; 
Psi(2,2)=-pbar*dbar'-dbar*pbar'+[0 0 lambda4; 0 lambda3*(2*l/pi)^4-2*lambda4 0; lambda4 0 0]; 
Psi(2,3)=-pbar; 
Psi(2,4)=-pbar; 
Psi(2,5)=-pbar; 
Psi(3,3)=-eta/L^2; 
Psi(4,4)=-gamma2; 
Psi(5,5)=-gamma2; 
Psi=sdpvar(Psi); 

PhiNabla=-2*D+2*(a-L+alpha)*P+(2*l/pi)^2*diag([lambda1,lambda2])+diag([lambda5,lambda6]); 
%% Solution of LMIs
LMIs=[Psi<=0, PhiNabla<=0, P>=0, gamma1>=0, lambda5>=0, lambda6>=0, diag([lambda1,lambda2,lambda3])>=eta*ones(3,3)]; 
options=sdpsettings('solver','lmilab','verbose',0);
sol=optimize(LMIs,[],options); 

flag=0; 
if sol.problem == 0
    [primal,~]=check(LMIs); 
    flag=min(primal)>=0; 
else
    yalmiperror(sol.problem) 
end