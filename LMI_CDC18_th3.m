function flag=LMI_CDC18_th3(D,a,L,N,cmax,l,h,alpha,gamma)
% This MATLAB program checks the feasibility of LMIs from Theorem 3 of the paper 
% A. Selivanov and E. Fridman, "Sampled-data H-infinity filtering of a 2D heat equation under pointwise measurements," in 57th Conference on Decision and Control, 2018. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% D, a  - parameters of (20) 
% L     - observer gain from (5) 
% N     - number of sensors 
% cmax  = max||c_i||_\infty with c_i from (21) 
% l     - subdomain size (13)
% h     - maximum sampling period 
% alpha - decay rate 
% gamma - L2-gain 

% Output: 
% flag =1 if feasible, =0 otherwise
%% Decision variables and notations 
sdpvar p1 p2 p3 gamma1 eta nu lambda1 lambda2 lambda3 lambda4 lambda5 lambda6
gamma2=gamma^2*gamma1; 
P=[p1 p2; p2 p3]; 
pbar=[p1; 2*p2; p3]; 
dbar=[D(1,1); 2*D(1,2); D(2,2)]; 
%% LMIs 
Upsilon=blkvar; 
Upsilon(1,1)=2*(a-L+alpha)-(lambda5+lambda6)*pi^2+gamma1; 
Upsilon(1,3)=1; 
Upsilon(1,4)=1; 
Upsilon(1,5)=1; 
Upsilon(1,6)=1; 
Upsilon(1,7)=nu*h*(a-L); 
Upsilon(2,2)=-pbar*dbar'-dbar*pbar'+[0 0 lambda4; 0 lambda3*(2*l/pi)^4-2*lambda4 0; lambda4 0 0]; 
Upsilon(2,3)=-pbar; 
Upsilon(2,4)=-pbar; 
Upsilon(2,5)=-pbar; 
Upsilon(2,6)=-pbar; 
Upsilon(2,7)=nu*h*dbar; 
Upsilon(3,3)=-eta/L^2; 
Upsilon(3,7)=nu*h; 
Upsilon(4,4)=-gamma2; 
Upsilon(4,7)=nu*h; 
Upsilon(5,5)=-gamma2; 
Upsilon(5,7)=nu*h; 
Upsilon(6,6)=-nu; 
Upsilon(6,7)=-nu*h; 
Upsilon(7,7)=-(pi^2*N*nu*exp(-2*alpha*h))/(4*L^2*cmax); 
Upsilon=sdpvar(Upsilon); 

PhiNabla=-2*D+2*(a-L+alpha)*P+(2*l/pi)^2*diag([lambda1,lambda2])+diag([lambda5,lambda6]); 
%% Solution of LMIs
LMIs=[Upsilon<=0, PhiNabla<=0, P>=0, gamma1>=0, lambda5>=0, lambda6>=0, diag([lambda1,lambda2,lambda3])>=eta*ones(3,3)]; 
options=sdpsettings('solver','lmilab','verbose',0);
sol=optimize(LMIs,[],options); 

flag=0; 
if sol.problem == 0
    [primal,~]=check(LMIs); 
    flag=min(primal)>=0; 
else
    yalmiperror(sol.problem) 
end