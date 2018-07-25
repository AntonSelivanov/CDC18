function flag=LMI_CDC18_th2(D,a,L,l,alpha,gamma)
% This MATLAB program checks the feasibility of LMIs from Theorem 2 of the paper 
% A. Selivanov and E. Fridman, "Sampled-data H-infinity filtering of a 2D heat equation under pointwise measurements," in 57th Conference on Decision and Control, 2018. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% D, a  - parameters of (1)
% L     - observer gain from (5)
% l     - subdomain size (13)
% alpha - decay rate 
% gamma - L2-gain

% Output: 
% flag =1 if feasible, =0 otherwise
%% Decision variables and notations 
sdpvar p1 p2 p3 q1 eta lambda1 lambda2 lambda3 lambda4 lambda5 lambda6
if gamma==0
    sdpvar q2
else
    q2=gamma^2*q1; 
end
P=[p1 p2; p2 p3]; 
d1=D(1,1); d2=(D(1,2)+D(2,1))/2; d3=D(2,2); 
%% LMIs 
Xi=blkvar; 
Xi(1,1)=2*(a-L+alpha)-(lambda5+lambda6)*pi^2+q1; 
Xi(1,5)=-L; 
Xi(1,6)=1; 
Xi(1,7)=-1; 
Xi(2,2)=-2*p1*d1; 
Xi(2,3)=-2*(p1*d2+p2*d1); 
Xi(2,4)=-p1*d3-p3*d1+lambda4; 
Xi(2,5)=L*p1; 
Xi(2,6)=-p1; 
Xi(2,7)=p1; 
Xi(3,3)=-8*p2*d2+lambda3*(2*l/pi)^4-2*lambda4; 
Xi(3,4)=-2*(p2*d3+p3*d2); 
Xi(3,5)=2*L*p2; 
Xi(3,6)=-2*p2; 
Xi(3,7)=2*p2; 
Xi(4,4)=-2*p3*d3; 
Xi(4,5)=L*p3; 
Xi(4,6)=-p3; 
Xi(4,7)=p3; 
Xi(5,5)=-eta; 
Xi(6,6)=-q2; 
Xi(7,7)=-q2; 
Xi=sdpvar(Xi); 

Psi=-2*D+2*(a-L+alpha)*P+(2*l/pi)^2*diag([lambda1,lambda2])+diag([lambda5,lambda6]); 
%% Solution of LMIs
LMIs=[Xi<=0, Psi<=0, P>=0, q1>=0, lambda5>=0, lambda6>=0, diag([lambda1,lambda2,lambda3])>=eta*ones(3,3)]; 
options=sdpsettings('solver','lmilab','verbose',0);
sol=optimize(LMIs,[],options); 

flag=0; 
if sol.problem == 0
    [primal,~]=check(LMIs); 
    flag=min(primal)>=0; 
else
    yalmiperror(sol.problem) 
end