%% Liquid properties

nu = 100e-6; %(m^2/s) kinematic viscosity
rho = 870; %(kg/m^3) density
K = 1.5e9; %(Pa) bulk modulus
c=sqrt(K/rho); % (m/s) wave speed

%% Pipeline Dimensions

L=1000;%(m) pipe length
rA=4*25.4e-3/2;%(m) pipe inner radius
rB=0.8*rA;

%% calc Transmission Matrix
N_cycles=50;
N_omega=100e3;
omega_cycle=c/(2*L)*2*pi;%(rad/s) full cycle frequency

omega=omega_cycle/N_cycles*(0:N_omega-1);

T = muto_transmission(omega,rA,rB,L,c,nu,rho);

%PA/QA with load resistance R_L
ZcA=rho*c/(pi*rA^2);%characteristic impedance
ZcB=rho*c/(pi*rB^2);%characteristic impedance
R_L=1/20*ZcB;%(Pa/(m^3/s)) load resistance
PA_QA_R=(squeeze(T(2,1,:))+R_L*squeeze(T(2,2,:)))./(squeeze(T(1,1,:))+R_L*squeeze(T(1,2,:)));

%% calculate impulse response
G=PA_QA_R;
G(1)=0;%fix DC gain
G_fft=[G; flipud(conj(G(2:end)))];%populate redundent part

U=ones(size(G_fft));%input (ones for unit impulse)
%Calculate frequency domain, given input from G(s)=X(s)/U(s)
X=G_fft.*U;

x_imp=ifft(X);%impulse response
x_imp=x_imp-mean(x_imp(round(N_omega*2*0.9)+(1:10)));%fix DC gain so that impulse response trails to zero near end

dt=2*pi/(omega(end)*2);%(s) sample time
t=((0:(numel(x_imp)-1))*dt)';%time

x_step=cumsum(x_imp);

%% plot
figure(1)
loglog(omega/(pi*c/L),abs(PA_QA_R)/ZcA)
ylabel('|P_A/Q_A|/Z_{cA}')
xlabel('\omega/(\pic/L)')
grid on

figure(2)
plot(t/(2*L/c),x_imp/ZcA)
grid on
ylabel('P_A/Z_{cA}/(1 m^3/s)')
xlabel('tc/(2L)')
title('Unit impulse response')

figure(3)
plot(t/(2*L/c),x_step/ZcA)
grid on
ylabel('P_A/Z_{cA}/(1 m^3/s)')
xlabel('tc/(2L)')
title('Unit step response')

