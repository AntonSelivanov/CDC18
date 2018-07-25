function flag=LMI_CDC18_th3(D,a,L,N,cmax,l,h,alpha,gamma)
% This MATLAB program checks the feasibility of LMIs from Theorem 3 of the paper 
% A. Selivanov and E. Fridman, "Sampled-data H-infinity filtering of a 2D heat equation under pointwise measurements," in 57th Conference on Decision and Control, 2018. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% D, a  - parameters of (1)
% L     - observer gain from (5)
% N     - number of sensors 
% cmax  = max||c_i||_\infty
% l     - subdomain size (13)
% h     - sampling period 
% alpha - decay rate 
% gamma - L2-gain

% Output: 
% flag =1 if feasible, =0 otherwise
%% Decision variables and notations 
sdpvar p1 p2 p3 q1 eta nu lambda1 lambda2 lambda3 lambda4 lambda5 lambda6
if gamma==0
    sdpvar q2
else
    q2=gamma^2*q1; 
end
P=[p1 p2; p2 p3]; 
d1=D(1,1); d2=(D(1,2)+D(2,1))/2; d3=D(2,2); 
%% LMIs 
Upsilon=blkvar; 
Upsilon(1,1)=2*(a-L+alpha)-(lambda5+lambda6)*pi^2+q1; 
Upsilon(1,5)=-L; 
Upsilon(1,6)=1; 
Upsilon(1,7)=-1; 
Upsilon(1,8)=-L; 
Upsilon(1,9)=nu*h*exp(alpha*h)*(a-L); 
Upsilon(2,2)=-2*p1*d1; 
Upsilon(2,3)=-2*(p1*d2+p2*d1); 
Upsilon(2,4)=-p1*d3-p3*d1+lambda4; 
Upsilon(2,5)=L*p1; 
Upsilon(2,6)=-p1; 
Upsilon(2,7)=p1; 
Upsilon(2,8)=L*p1; 
Upsilon(2,9)=nu*h*exp(alpha*h)*d1; 
Upsilon(3,3)=-8*p2*d2+lambda3*(2*l/pi)^4-2*lambda4; 
Upsilon(3,4)=-2*(p2*d3+p3*d2); 
Upsilon(3,5)=2*L*p2; 
Upsilon(3,6)=-2*p2; 
Upsilon(3,7)=2*p2; 
Upsilon(3,8)=2*L*p2; 
Upsilon(3,9)=nu*h*exp(alpha*h)*2*d2; 
Upsilon(4,4)=-2*p3*d3; 
Upsilon(4,5)=L*p3; 
Upsilon(4,6)=-p3; 
Upsilon(4,7)=p3; 
Upsilon(4,8)=L*p3; 
Upsilon(4,9)=nu*h*exp(alpha*h)*d3; 
Upsilon(5,5)=-eta; 
Upsilon(5,9)=-nu*h*exp(alpha*h)*L; 
Upsilon(6,6)=-q2; 
Upsilon(6,9)=nu*h*exp(alpha*h); 
Upsilon(7,7)=-q2; 
Upsilon(7,9)=-nu*h*exp(alpha*h); 
Upsilon(8,8)=-nu*pi^2/4; 
Upsilon(8,9)=-nu*h*exp(alpha*h)*L; 
Upsilon(9,9)=-N*nu/cmax; 
Upsilon=sdpvar(Upsilon); 

Psi=-2*D+2*(a-L+alpha)*P+(2*l/pi)^2*diag([lambda1,lambda2])+diag([lambda5,lambda6]); 
%% Solution of LMIs
LMIs=[Upsilon<=0, Psi<=0, P>=0, q1>=0, lambda5>=0, lambda6>=0, diag([lambda1,lambda2,lambda3])>=eta*ones(3,3)]; 
options=sdpsettings('solver','lmilab','verbose',0);
sol=optimize(LMIs,[],options); 

flag=0; 
if sol.problem == 0
    [primal,~]=check(LMIs); 
    flag=min(primal)>=0; 
else
    yalmiperror(sol.problem) 
end