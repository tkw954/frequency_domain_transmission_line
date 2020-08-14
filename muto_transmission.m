function T = muto_transmission(omega,rA,rB,l,c,nu,rho)
%T= muto_transmission(omega,rA,rB,l,c,nu,rho)
% Calculates transmission matrix:
% [qA;pA]=T*[qB;pB]
% for inputs of:
% omega: frequency
% rA,rB: pipe radius at A and B end
% l: length
% c: wave speed
% nu: kinematic viscosity
% rho: density
% Follows Muto, Takayoshi, Yoshifumi Kinoshita, and Ryuichi Yoneda. 
% "Dynamic Response of Tapered Fluid Lines: 1st Report, Transfer Matrix 
% and Frequency Response." Bulletin of JSME 24.191 (1981): 809-815.
% DOI:10.1299/jsme1958.24.809

R0=rA;
theta=atan((rB-rA)./l);


s=1j*omega;
S=l./c*s;

t0=l./c;
D=nu*t0./R0^2;
chi0=1j*sqrt(S/D);
a0=pi*rA.^2;
gamma0=sqrt(-besselj(0,chi0,1)./besselj(2,chi0,1));%eq 13
lambda0=1/2*(1-(chi0.*(gamma0.^2-1)./(2*gamma0)).^2);%eq 17
Omegal=s.*gamma0.*l./c.*(1+theta*l./R0.*(1-lambda0));%eq 22



%eq 21
T_11=cosh(Omegal)+c.*lambda0.*theta./(R0.*s.*gamma0).*sinh(Omegal);
T_12=a0./(gamma0.*rho.*c).*sinh(Omegal);
T_21=(gamma0.*rho.*c)./a0.*sinh(Omegal);
T_22=cosh(Omegal)-c.*lambda0.*theta./(R0.*s.*gamma0).*sinh(Omegal);

kq=1./(1+lambda0.*theta.*l./R0);
kp=1./(1-lambda0.*theta.*l./R0);


T=[reshape(T_11.*kq,1,1,[]) reshape(T_12.*kp,1,1,[]);
    reshape(T_21.*kq,1,1,[]) reshape(T_22.*kp,1,1,[])];



end

