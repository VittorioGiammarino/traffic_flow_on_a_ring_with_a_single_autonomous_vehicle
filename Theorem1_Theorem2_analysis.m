%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Source Code for Paper: Traffic flow on a ring with a single autonomous
%%%                        vehicle: an interconnected perspective
%%% Journal: Transactions on Intelligent transportation Systems, 2019
%%% Authors: Vittorio Giammarino, Simone Baldi, Paolo Frasca and Maria
%%%          Laura Delle Monache
%%%
%%% Copyright Vittorio Giammarino, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
length_standard = 4.5; %vehicle length
N = 22; % Number of vehicles in the platoon
Ring_length = (260*N)/22;
v_max = 9.751; % Max speed in m/s
h_eq = (Ring_length/N); % equilibrium headway
d_s = 6; % safe distance
V = @(h,l_v,d_s) v_max*((tanh(h-l_v-d_s) + tanh(l_v + d_s))/(1 + tanh(l_v + d_s))); %Velocity non-linear function. h = headway, l_v = Veh length, d_s = safe distance
dV = @(h,l_v,d_s) (v_max/(1+tanh(l_v + d_s)))*(1-(tanh(h - l_v -d_s))^2); % dV/dh derivative of V in the headway
v_eq = V(h_eq, length_standard, d_s); % equilibrium velocity of the ring
k_bar = dV(h_eq, length_standard, d_s); % HV param
b_bar = 0.5; % HV param
a_bar = 20/(h_eq)^2; % HV param
%% Theorem 1: Necessary and sufficient condition for stability
alpha_1 = b_bar*k_bar;
alpha_2 = a_bar+b_bar;
alpha_3 = a_bar;

delta_lambda_plus = zeros(N,1);
delta_lambda_minus = zeros(N,1);

for i =1:N
    gamma = alpha_2 -alpha_3*cos(2*pi*(i-1)/N);
    phi = -alpha_3*sin(2*pi*(i-1)/N);
    eta = 4*alpha_1*(sin(2*pi*(i-1)/N));
    r = gamma^2 - phi^2 - 4*alpha_1*(1 - cos(2*pi*(i-1)/N));

    if (eta-2*gamma*phi) == 0
        delta_lambda_plus(i,1) = 0.5*(-gamma + sqrt(r));
        delta_lambda_minus(i,1) = 0.5*(-gamma - sqrt(r));
    else
        delta_lambda_plus(i,1) = 0.5*(-gamma + sqrt(0.5*(sqrt(r^2 + (eta +2*gamma*phi)^2) + r)));
        delta_lambda_minus(i,1) = 0.5*(-gamma - sqrt(0.5*(sqrt(r^2 + (eta +2*gamma*phi)^2) + r)));
    end
end

delta_lambda = [delta_lambda_plus; delta_lambda_minus]; % real part of the eigenvalues of the system

%% Transfer functions for P_N and Gamma_i (See Fig.2 and Fig.3)
a=20; %HV param
b=0.5; %HV param
a_bar = a/(h_eq*h_eq); %HV param linearized
b_bar = b; %HV param linearized
% transfer functions of P_N and Gamma_i
Pd_num = [1 0];
Pd_den = [1 (a_bar+b_bar) b_bar*k_bar];
P_num = [a_bar b_bar*k_bar];
P_den = [1 (a_bar+b_bar) b_bar*k_bar];
PP = tf(P_num,P_den); % Gamma
PP_d = tf(Pd_num,Pd_den); % P_N

[p_gamma, z_gamma] = pzmap(PP);
[p_HD_ring,z_HD_ring] = pzmap((PP_d)/(1-PP^N));
% Bode diagram for HVs
figure()
bode(PP)
grid on
title("Bode diagram for HVs")
% Pole and zeros distribution of Gamma 
figure()
plot(p_gamma, 'rx')
hold on
plot(z_gamma,0 , 'bo')
hold on
grid on
title("Pole and zeros distribution of Gamma")
%% Poles distribution of HVs on the ring for different sizes (different number of vehicles N)

[p_HD_ring_5,z_HD_ring_5] = pzmap((PP_d)/(1-PP^5));
[p_HD_ring_10,z_HD_ring_10] = pzmap((PP_d)/(1-PP^10));
[p_HD_ring_15,z_HD_ring_15] = pzmap((PP_d)/(1-PP^15));
[p_HD_ring_20,z_HD_ring_20] = pzmap((PP_d)/(1-PP^20));
[p_HD_ring_25,z_HD_ring_25] = pzmap((PP_d)/(1-PP^25));

figure()
plot3(imag(p_HD_ring_5),real(p_HD_ring_5),5*ones(length(p_HD_ring_5),1),'rd')
hold on
plot3(imag(p_HD_ring_10),real(p_HD_ring_10),10*ones(length(p_HD_ring_10),1),'mo')
hold on
plot3(imag(p_HD_ring_15),real(p_HD_ring_15),15*ones(length(p_HD_ring_15),1),'c*')
hold on
plot3(imag(p_HD_ring_20),real(p_HD_ring_20),20*ones(length(p_HD_ring_20),1),'k.')
grid on
plot3(imag(p_HD_ring_25),real(p_HD_ring_25),25*ones(length(p_HD_ring_25),1),'bx')
grid on
xlabel('Imm axis')
ylabel('Real axis')
zlabel('Number of Vehicles')
ylim([-0.8 0.3])
xlim([-1.5 1.5])
zlim([0 30])
legend('N=5', 'N=10', 'N=15', 'N=22', 'N=25')
title("Poles distribution for different number of HVs N on the ring")

figure()
plot(real(p_HD_ring_5),imag(p_HD_ring_5),'rd')
hold on
plot(real(p_HD_ring_10),imag(p_HD_ring_10),'mo')
hold on
plot(real(p_HD_ring_15),imag(p_HD_ring_15),'c*')
hold on
plot(real(p_HD_ring_20),imag(p_HD_ring_20),'k.')
grid on
plot(real(p_HD_ring_25),imag(p_HD_ring_25),'bx')
grid on
xlabel('Real axis')
ylabel('Imaginary axis')
legend('N=5','N=10','N=15','N=20','N=25')
ylim([-1.2 1.2])
xlim([-1 0.5])
title("Poles distribution for different number of HVs N on the ring")

%% AV Design for Lyapunov stability for N vehicles
N = 22;
%HV
alpha_1 = b_bar*k_bar;
alpha_2 = a_bar + b_bar;
alpha_3 = a_bar;

% Compute critical peak and frequency
fpeaK_PP_formula = sqrt(alpha_1*(1-2*(alpha_2/(2*sqrt(alpha_1)))^2));
dumping = alpha_2/(2*sqrt(alpha_1));
gpeak_PP_formula = 1/(2*dumping*sqrt(1-dumping^2));
[gpeak_PP,fpeak_PP] = getPeakGain(PP,1e-5);

% AV design
beta = 1;
alpha = 0.9;
delta = 23;
k_veh = abs((fpeak_PP*((gpeak_PP)^(1-N)))/(1-alpha/2))*sqrt(1/(1-(gpeak_PP)^(2-2*N))); % N>4 and first order AV design
%k_veh = 15; % Uncomment this value for N=4 or less

beta_1 = k_veh*alpha/delta;
beta_2 = k_veh*(1-alpha/2);
beta_3 = beta_2;

%AV first order
AV_N1 = beta_3;
AV_D1 = [1 beta_2];
AV1 = tf(AV_N1,AV_D1);

%AV second order (PI with saturation proposed in [16]) 
AV_N = [beta_3 beta_1];
AV_D = [1 beta_2 beta_1];
AV = tf(AV_N,AV_D);

%% Bode
figure()
bode(AV1,AV,PP)
legend('AV first order', 'PI with saturation proposed in [16]', 'HV')

%% Theorem 2 analysis
% AV first order
figure()
bode(AV1,PP)
legend('AV first order','HV')
title('Bode AV first order and HVs')

figure()
bode(PP^(N-1)*AV1)
legend('bode of N-1 HVs and AV first order')
title('AV first order action on N-1 HVs')

[p_AV1_ring,z_AV1_ring] = pzmap((PP_d)/(1-AV1*PP^(N-1)));
figure()
plot(p_AV1_ring,'d')
grid on
title('ring poles with AV first order')

%% Theorem 2 analysis
% PI with saturation proposed in [16]
%Checks AV2
figure()
bode(AV,PP)
legend('PI with saturation proposed in [16]','HV')
title('Bode PI with saturation proposed in [16] and HVs')

figure()
bode(PP^(N-1)*AV)
legend('bode of N-1 HVs and PI with saturation proposed in [16]')
title('PI with saturation proposed in [16] action on N-1 HVs')

[p_AV_ring,z_AV_ring] = pzmap((PP_d)/(1-AV*PP^(N-1)));

figure()
plot(p_AV_ring,'d')
grid on
title('ring poles with PI with saturation proposed in [16]')
