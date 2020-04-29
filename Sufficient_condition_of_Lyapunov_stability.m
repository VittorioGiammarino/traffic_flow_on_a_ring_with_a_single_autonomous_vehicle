%% Initialization
length_standard = 4.5;
N = 22; % either 21 or 22 according to which experiment is considered
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

% b_bar = 0.5;
% a_bar = 20/(h_eq)^2;

%%
a=140;
b=0.1;

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

[p_gamma, z_gamma] = pzmap(PP);
[p_HD_ring,z_HD_ring] = pzmap((PP_d)/(1-PP^N));

%HD
alpha_1 = b_bar*k_bar;
alpha_2 = a_bar + b_bar;
alpha_3 = a_bar;

figure()
bode(PP)
grid on

% figure()
% plot(p_gamma, 'rx')
% hold on
% plot(z_gamma,0 , 'bo')
% hold on
% grid on

figure()
plot(p_HD_ring,'kx')
grid on

absolute = (sqrt(alpha_1^2 + (alpha_3.*imag(p_HD_ring)).^2))./(sqrt((alpha_1-(imag(p_HD_ring)).^2).^2+(alpha_1.*imag(p_HD_ring)).^2));

figure()
plot(abs(p_HD_ring))
grid on

%%
[x,y] = meshgrid(-3 : 0.05: 3);
s = x + 1i*y;
z=abs((alpha_1 + alpha_3*s)./(s.^2 + alpha_2.*s + alpha_1));
figure;
surf(x, y, z)
hold on
surf(x,y,1)

%%

[x,y] = meshgrid(0 : 0.05: 3);
s = x + 1i*y;
z=abs((alpha_1 + alpha_3*s)./(s.^2 + alpha_2.*s + alpha_1));
figure;
surf(x, y, z)

%% tanh
x=0:0.001:6;
y=2*ones(1,length(x));
figure()
plot(x,y,'k--')
hold on
plot(x,tanh(x-2)+tanh(2))
grid on
ylim([0 2.1])
xlabel('\Deltax')
ylabel('V(\Deltax)')












