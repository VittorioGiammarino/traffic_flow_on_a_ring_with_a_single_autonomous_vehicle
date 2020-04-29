%% Initialization

length_standard = 4.5;
N = 4; % either 21 or 22 according to which experiment is considered
Ring_length = (260*N)/22;
v_max = 9.751; % From reference 3.

% From: "Stabilizing Traffic Flow via a Single Autonomous Vehicle:
% Possibilities and Limitations" 
% (x_j)''= f(h_j, (h_j)', v_j) = b_bar*(V(h_j)-v_j) + a_bar*(h_j)'
% if equilibrium (x_j)''= 0 --> f(h_eq, 0, v_eq) = 0 --> V(h_eq) = v_eq

h_eq = (Ring_length/N);

% at equilibrium the acceleration is zero and therefore v_j+1 = v_j
d_s = 6; % reference 3 depicts: max(2sec*(v_j+1 - v_j), 4m) 

V = @(h,l_v,d_s) v_max*((tanh(h-l_v-d_s) + tanh(l_v + d_s))/(1 + tanh(l_v + d_s))); % h = headway, l_v = Veh length, d_s = safe distance
dV = @(h,l_v,d_s) (v_max/(1+tanh(l_v + d_s)))*(1-(tanh(h - l_v -d_s))^2); % dV/dh derivative of V in the headway

v_eq = V(h_eq, length_standard, d_s); % from here we can assume v_eq=9.7 m/s
k_bar = dV(h_eq, length_standard, d_s); 

Human_Driver_tf = @(a_bar, b_bar, k_bar) tf([a_bar b_bar*k_bar],[1 (a_bar+b_bar) b_bar*k_bar]); % Human Driver OV-FTL transfer function
Autonomous_Veh_tf = @(k_veh, alpha) tf([0 k_veh*(1-(alpha/2))],[1 k_veh*(1-(alpha/2))]); % Autonomous Vehicle OV-FTL transfer function
%%
a=20;
b=0.5;

Autonomous_Veh_second_order = @(k_veh, alpha, beta) tf([k_veh*(1-(alpha/2)) (k_veh*alpha*beta)/23],[1 k_veh*(1-(alpha/2)) (k_veh*alpha*beta)/23]); % Autonomous Vehicle OV-FTL second order transfer function

% calculation of critical N on a string
a_bar = a/(h_eq*h_eq);
b_bar = b;
% transfer functions of Pd and P
Pd_num = [1 0];
Pd_den = [1 (a_bar+b_bar) b_bar*k_bar];
P_num = [a_bar b_bar*k_bar];
P_den = [1 (a_bar+b_bar) b_bar*k_bar];
PP = tf(P_num,P_den);
PP_d = tf(Pd_num,Pd_den);
% get the peak frequency of P and compare it with Pd at the same frequency
%%

% for i=22:23
% figure(1)
% bode(PP_d/(1-PP^i))
% hold on
% end
% 
% [p1,z] = pzmap((PP_d)/(1-PP^3));
% [p2,z2] = pzmap((PP_d)/(1-PP^10));
% [p3,z3] = pzmap((PP_d)/(1-PP^20));
% [p4,z4] = pzmap((PP_d)/(1-PP^30));
% [p5,z5] = pzmap((PP_d)/(1-PP^50));
% [p6,z6] = pzmap((PP_d)/(1-PP^100));
% 
% figure()
% plot(p1,'rx')%,p2,'bx',p3,'kx',p4,'mx')
% hold on
% plot(p2,'bx')
% hold on
% plot(p3,'kx')
% % hold on
% % plot(p4,'yx')
% % % hold on
% % % plot(p5,'gx')
% % % hold on
% % % plot(p6,'x')
% grid on

figure()
sigma((PP_d)/(1-PP^N))
hold on


figure()
sigma((PP_d)/(1-PP^N))
hold on
sigma(((PP_d)*PP)/(1-PP^N))
hold on
sigma(((PP_d)*PP^2)/(1-PP^N))
hold on
sigma(((PP_d)*PP^6)/(1-PP^N))
legend('F^{(N)}N','F^{(N)}N-1','F^{(N)}N-2','F^{(N)}N-6')
grid on


figure()
sigma((PP_d)/(1-PP^N),{10^-1,10})
hold on
sigma(((PP_d)*PP)/(1-PP^N),{10^-1,10})
hold on
sigma(((PP_d)*PP^2)/(1-PP^N),{10^-1,10})
hold on
sigma(((PP_d)*PP^6)/(1-PP^N),{10^-1,10})
legend('F^{(N)}_N','F^{(N)}_{N-1}','F^{(N)}_{N-2}','F^{(N)}_{N-6}')
grid on

%%

[peakF_omega_f, omega_f] = getPeakGain((PP_d)/(1-PP^N));
[peakGamma_omega_gamma, omega_gamma] = getPeakGain(PP);
[peakF_omega_gamma, omega_gamma_2] = getPeakGain((PP_d)/(1-PP^N),0.01,[omega_gamma-0.01, omega_gamma+0.01]);

N_bar = ceil((log(peakF_omega_f)-log(peakF_omega_gamma))/log(peakGamma_omega_gamma));



%%
% 
% [peakF_1_omega_gamma, omega_gamma_2] = getPeakGain((PP_d*PP^(N-1))/(1-PP^N),0.01,[omega_gamma-0.01, omega_gamma+0.01]);
% [peakF_1_omega_f_2, omega_f_2] = getPeakGain((PP_d*PP^(N-1))/(1-PP^N),0.01,[omega_f-0.01, omega_f+0.01]);
% 
% [peakGamma_omega_f, omega_f3] = getPeakGain(PP,0.01,[omega_f-0.01, omega_f+0.01]);


[p1,z] = pzmap((PP_d)/(1-PP^N));
[p_gamma, z_gamma] = pzmap(PP);
%[p2,z2] = pzmap((PP_d)/(1-AV*PP^N));
%[p3,z3] = pzmap((PP_d)/(1-PP^22));
%[p4,z4] = pzmap((PP_d)/(1-AV^22));
% [p5,z5] = pzmap((PP_d)/(1-AV*PP^50));
% [p6,z6] = pzmap((PP_d)/(1-AV*PP^100));

figure()
plot(p1,'kx')%,p2,'bx',p3,'kx',p4,'mx')
hold on
% plot(z,'ko')
% hold on
grid on
xlabel('Real axis')
ylabel('Imaginary axis')
xlim([-0.7 0.2])
ylim([-1.5 1.5])

figure()
plot(p_gamma, 'rx')
hold on
plot(z_gamma,0 , 'bo')
hold on
grid on

% figure()
% plot(p2,'bx')
% hold on
% plot(p3,'kx')
% hold on
% plot(p4,'go')
% hold on
% % plot(p5,'gx')
% % hold on
% % plot(p6,'x')
grid on


figure()
sigma((PP_d)/(1-PP^N))
hold on
sigma(PP)
hold on
legend('F^{(N)}N','\Gamma')
grid on
% bode(PP^2)

% bode((PP_d*PP^2)/(1-PP^N))
% hold on
% bode(PP)
% hold on
% bode(PP_d)
% hold on
% bode(1/(1-PP))
% hold on
% %legend('(P_{N})/(1-\Gamma^{3})','(P_{N}*\Gamma)/(1-\Gamma^{3})','(P_{N}*\Gamma^{2})/(1-\Gamma^{3})') %,'(P_{N}*\Gamma^{6})/(1-\Gamma^{22})')
% 
% figure()
% bode(PP)
% hold on
% bode(1-PP)
% hold on

%%
A = zeros(2*N,2*N);

for i=2:2:(2*N-2)
    A(i-1,i) = 1;
    A(i,i-1) = -b_bar*k_bar;
    A(i,i) = -a_bar - b_bar;
    A(i,i+1) = b_bar*k_bar;
    A(i,i+2) = a_bar;
end

A(2*N-1,2*N) = 1;
A(2*N,1) = b_bar*k_bar;
A(2*N,2) = a_bar;
A(2*N,2*N-1) = -b_bar*k_bar;
A(2*N,2*N) = -a_bar - b_bar;

e = eig(A);


%% AV

a=20;
b=0.5;
c = 0.5;
beta = 1;
alpha = 0.9;
delta = 23;

a_bar = a/(h_eq*h_eq);
b_bar = b;

alpha_1 = b_bar*k_bar;
alpha_2 = a_bar + b_bar;
alpha_3 = a_bar;

% transfer functions of Pd and P
Pd_num = [1 0];
Pd_den = [1 (a_bar+b_bar) b_bar*k_bar];
P_num = [a_bar b_bar*k_bar];
P_den = [1 (a_bar+b_bar) b_bar*k_bar];
PP = tf(P_num,P_den);
PP_d = tf(Pd_num,Pd_den);

fpeaK_PP_formula = sqrt(alpha_1*(1-2*(alpha_2/(2*sqrt(alpha_1)))^2));
dumping = alpha_2/(2*sqrt(alpha_1));
gpeak_PP_formula = 1/(2*dumping*sqrt(1-dumping^2));
[gpeak_PP,fpeak_PP] = getPeakGain(PP,1e-5);
%k_veh = abs((fpeak_PP*((gpeak_PP)^(1-N)))/(1-alpha/2))*sqrt(1/(1-(gpeak_PP)^(2-2*N)));
k_veh = 15;

%AV
beta_1 = k_veh*alpha/delta;
beta_2 = k_veh*(1-alpha/2);
beta_3 = beta_2;

AV_N = [beta_3 beta_1];
AV_D = [1 beta_2+c beta_1];

AV_N1 = beta_3;
AV_D1 = [1 beta_2];

AV = tf(AV_N,AV_D);
[gpeak,fpeak] = getPeakGain(AV,1e-5);

AV1 = tf(AV_N1,AV_D1);

c=0.5;
AV_N1_c = beta_3;
AV_D1_c = [1 beta_2+c];
AV1_c = tf(AV_N1_c,AV_D1_c);

%%
c=0.5;

% K_alpha = (c + sqrt(((gpeak_PP)^(2*N-2))*(c^2+fpeak_PP^2) - fpeak_PP^2))/((gpeak_PP)^(2*N-2)-1); 
% k_veh2 =K_alpha/(1-alpha/2);

K_alpha = k_veh*(1-alpha/2);
k_veh2 =K_alpha/(1-alpha/2);

beta_1 = k_veh2*alpha/delta;
beta_2 = K_alpha;
beta_3 = beta_2;

%AV First Order with c
AV_N_c = [beta_3 beta_1];
AV_D_c = [1 beta_2+c beta_1];
AV_c = tf(AV_N_c,AV_D_c);

%% Figures
% figure()
% bode(AV)
% hold on
% bode(PP)
% hold on
% figure()
% bode(PP^21*AV)

% figure()
% bode(AV1)
% hold on
% bode(PP)
% hold on
% figure()
% bode(PP^21*AV1)
 
figure()
sigma((PP_d)/(1-AV*PP^(N-1)))
hold on
sigma(PP)
grid on
legend('F_N^{AV(N)} with K_{veh}=15','\Gamma')


figure()
sigma((PP_d)/(1-AV*PP^(N-1)))
hold on
sigma(PP)
hold on
sigma((PP_d)/(1-AV_c*PP^(N-1)))
grid on
legend('F_N^{AV(N)} with K_{veh}=15','\Gamma','F_N^{AV(N)} with K_{veh}=0.8723')
% bode((PP_d*PP^1)/(1-AV*PP^21))
% hold on
% bode((PP_d*PP^2)/(1-AV*PP^21))
% hold on
% bode((PP_d*PP^6)/(1-AV*PP^21))
% grid on
% legend('(P_{N})/(1-\Gamma^{3})','(P_{N}*\Gamma)/(1-\Gamma^{3})','(P_{N}*\Gamma^{2})/(1-\Gamma^{3})')
% % 
% 
% 
% figure(2)
% pzmap((PP_d)/(1-AV*PP^22))
% hold on

figure()
sigma((PP_d)/(1-AV1*PP^(N-1)))
hold on
sigma((PP_d*PP^(1))/(1-AV1*PP^(N-1)))
hold on
sigma((PP_d*PP^(2))/(1-AV1*PP^(N-1)))
hold on
sigma((PP_d*PP^(3))/(1-AV1*PP^(N-1)))
grid on
legend('AV-veh N', 'veh N-1', 'veh N-2', 'veh N-3')

figure()
sigma((PP_d)/(1-AV1_c*PP^(N-1)))
hold on
sigma((PP_d*PP^(1))/(1-AV1_c*PP^(N-1)))
hold on
sigma((PP_d*PP^(2))/(1-AV1_c*PP^(N-1)))
hold on
sigma((PP_d*PP^(3))/(1-AV1_c*PP^(N-1)))
grid on
legend('AV-veh N', 'veh N-1', 'veh N-2', 'veh N-3')

figure()
sigma((PP_d)/(1-AV*PP^(N-1)))
hold on
sigma((PP_d*PP^(1))/(1-AV*PP^(N-1)))
hold on
sigma((PP_d*PP^(2))/(1-AV*PP^(N-1)))
hold on
sigma((PP_d*PP^(3))/(1-AV*PP^(N-1)))
grid on
legend('AV-veh N', 'veh N-1', 'veh N-2', 'veh N-3')

figure()
sigma((PP_d)/(1-AV_c*PP^(N-1)),{10^-1,10})
hold on
sigma((PP_d*PP^(1))/(1-AV_c*PP^(N-1)),{10^-1,10})
hold on
sigma((PP_d*PP^(2))/(1-AV_c*PP^(N-1)),{10^-1,10})
hold on
sigma((PP_d*PP^(3))/(1-AV_c*PP^(N-1)),{10^-1,10})
grid on
legend('F^{AV(N)}_N','F^{AV(N)}_{N-1}','F^{AV(N)}_{N-2}','F^{AV(N)}_{N-3}')


[p1,z] = pzmap((PP_d)/(1-AV*PP^(N-1)));
[p2,z2] = pzmap((PP_d)/(1-PP^N));
[p3,z3] = pzmap((PP_d)/(1-AV*PP^3));
[p4,z4] = pzmap((PP_d)/(1-AV*PP^4));
[p5,z5] = pzmap((PP_d)/(1-AV^22));
[p6,z6] = pzmap((PP_d)/(1-AV*PP^21));

 figure()
plot(p1,'rx')%,p2,'bx',p3,'kx',p4,'mx')
hold on
%plot(z,'ko')
% plot(p2,'bo')
% hold on
% plot(p3,'k*')
% hold on
% plot(p4,'+')
% hold on
% plot(z6,'.')
% hold on
% plot(p6,'d')
grid on

figure()
plot(p1,'kx')
grid on
xlabel('Real axis')
ylabel('Imm axis')
%xlim([-1 0.5])
ylim([-1.5 1.5])

% figure()
% bode((PP_d)/(1-AV*PP^21))
% hold on
% bode(AV)

%%
x=-0.2519;
y=1.031;

absolute_H = sqrt((alpha_1+alpha_3*x)^2+(alpha_3*y)^2)/...
    sqrt((alpha_1 + alpha_2*(x)+(x)^2-(y)^2)^2+(alpha_2*(y)+2*(x)*(y))^2);

absolute_AV = sqrt((beta_1+beta_3*x)^2+(beta_3*y)^2)/...
    sqrt((beta_1 + beta_2*(x)+(x)^2-(y)^2)^2+(beta_2*(y)+2*(x)*(y))^2);

absolute_H^(21/22)*absolute_AV^(1/22)

%% Diagonalization HD-platoon

En = [1 1 1; 1 exp(1j*2*pi/N) exp(1j*2*2*pi/N); 1 exp(1j*2*2*pi/N) exp(1j*4*2*pi/N)];

%%
I2 = eye(2);

F_T = kron(En',I2);
F = kron(En, I2);

A = zeros(2*N,2*N);

for i=2:2:(2*N-2)
    A(i-1,i) = 1;
    A(i,i-1) = -b_bar*k_bar;
    A(i,i) = -a_bar - b_bar;
    A(i,i+1) = b_bar*k_bar;
    A(i,i+2) = a_bar;
end

A(2*N-1,2*N) = 1;
A(2*N,1) = b_bar*k_bar;
A(2*N,2) = a_bar;
A(2*N,2*N-1) = -b_bar*k_bar;
A(2*N,2*N) = -a_bar - b_bar;

Sigma = F_T*A*F;

%% Diagonalization AV-platoon

En = [1 1; 1 exp(1j*2*pi/N)];
I2 = eye(2);

N = 3;

F_T = kron(En',I2);
F = kron(En, I2);

A_AV = zeros(2*N,2*N);

for i=2:2:(2*N-2)
    A_AV(i-1,i) = 1;
    A_AV(i,i-1) = -b_bar*k_bar;
    A_AV(i,i) = -a_bar - b_bar;
    A_AV(i,i+1) = b_bar*k_bar;
    A_AV(i,i+2) = a_bar;
end

        A_AV(2*N-1,2*N) = 1;
        A_AV(2*N,1) = beta_1;
        A_AV(2*N,2) = beta_3;
        A_AV(2*N,2*N-1) = -beta_1;
        A_AV(2*N,2*N) = -beta_2;
        
%Sigma_AV = F_T*A_AV*F;

[V1,D] = eig(A_AV)




%%
a=20;
b=0.5;
k_veh = 10;
%beta = 1;
alpha = 0.9;
delta = 23;

a_bar = a/(h_eq)^2;
b_bar = b;
v = 0:0.001:1;

alpha_1 = b_bar*k_bar;
alpha_2 = a_bar + b_bar;
alpha_3 = a_bar;

beta_1 = k_veh*alpha/delta;
beta_2 = k_veh*(1-alpha/2);
beta_3 = beta_2;

AV_design = zeros(length(v),1);
HV_design = zeros(length(v),1);

for i=1:length(v)
AV_design(i,1) = (beta_1^2 + beta_3^2*v(i))/((beta_1 - v(i)^2)^2+beta_2^2*v(i)^2);
HV_design(i,1) = (((alpha_1-v(i)^2)^2 + alpha_2^2*v(i)^2)/(alpha_1^2+alpha_3^2*v(i)^2))^(N-1);
end
figure()
plot(v,AV_design,v,HV_design)


%% Ring Stability with AV

a = 100;%0.14*h_eq*h_eq;
b = 0.5;
k_veh = 1;
beta = 1;
alpha = 0.9;

% calculation of critical N on a string
a_bar = a/(h_eq*h_eq);
b_bar = b;
% transfer functions of Pd and P
Pd_num = [1 0];
Pd_den = [1 (a_bar+b_bar) b_bar*k_bar];
P_num = [a_bar b_bar*k_bar];
P_den = [1 (a_bar+b_bar) b_bar*k_bar];
PP = tf(P_num,P_den);
PP_d = tf(Pd_num,Pd_den);

AV_N = [k_veh*(1-(alpha/2)) (k_veh*alpha*beta)/23];
AV_D = [1 k_veh*(1-(alpha/2)) (k_veh*alpha*beta)/23];

AV = tf(AV_N,AV_D);
[gpeak,fpeak] = getPeakGain(AV,1e-5)

figure()
bode(AV)

figure()
bode((PP_d*PP)/(1-AV*PP^21))
hold on
bode((PP_d*PP^2)/(1-AV*PP^21))
hold on
bode((PP_d*PP^4)/(1-AV*PP^21))
hold on
bode((PP_d*PP^21)/(1-AV*PP^21))
grid on
legend('(P_{N}*\Gamma)/(1-\Gamma_{AV}*\Gamma^{21})','(P_{N}*\Gamma^{2})/(1-\Gamma_{AV}*\Gamma^{21})','(P_{N}*\Gamma^{4})/(1-\Gamma_{AV}*\Gamma^{21})','(P_{N}*\Gamma^{21})/(1-\Gamma_{AV}*\Gamma^{21})')



