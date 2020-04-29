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

a1=100;
b1=0.1;
% calculation of critical N on a string
a1_bar = a1/(h_eq*h_eq);
b1_bar = b1;
P1_num = [a1_bar b1_bar*k_bar];
P1_den = [1 (a1_bar+b1_bar) b1_bar*k_bar];
PP1 = tf(P1_num,P1_den);

figure()
sigma(PP_d,{10^-1,10})
hold on
sigma(((PP_d)*PP),{10^-1,10})
hold on
sigma(((PP_d)*PP*PP1),{10^-1,10})
hold on
sigma(((PP_d)*PP*PP1*PP1),{10^-1,10})
legend('P_4','P_3(a=20,b=0.5)','P_2(a=100,b=0.1)','P_1(a=100,b=0.1)')
title('strong string unstable but weak string stable for i in {1,2,3}')
grid on

%% Non Linear simulation N cars with proposed AV and Disturbance on the AV 
a=20;
b=0.5;
a1=100;
b1=0.1;
v_max = 9.751;
l_v = 4.5;
d_s = 6;
N = 4; %from here I change the numbers of vehicles into the Ring
Ring_length = 260*N/22;
spacing = Ring_length/N;
V = @(x) v_max*((tanh(x-l_v-d_s) + tanh(l_v + d_s))/(1 + tanh(l_v + d_s))); 
Set_of_initial_conditions = v_eq; 
t_f=40; 

for z=1:length(Set_of_initial_conditions)
    
opts = odeset('MaxStep',1e-2);    
y0 = Initial_velocity_and_space_conditions(Set_of_initial_conditions(z),N,spacing);

[t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_String_with_dist_Review(t,y,N,a,b,a1,b1,Ring_length,V),[0 t_f],y0,opts);
[displacement_AV, velocity_AV] = Plot_Displacement_velocity(t,y_N_cars,N,Ring_length);
    
% Uncomment to plot also the positions
% figure()
% plot(t,displacement_AV(:,end),'ro',t,displacement_AV(:,1:end-1),'ko')
% xlabel('time [s]')
% ylabel('position [m]')
% xlim([100 140])
% legend('AV','HVs')
% %title({'displacement N cars with AV', ['Initial velocity condition y_0 =' num2str(Set_of_initial_conditions(z))]})

figure()
plot(t,velocity_AV(:,1:end-1),'k',t,velocity_AV(:,end),'r')%,t,velocity_AV(:,4),t,velocity_AV(:,1))
xlabel('time [s]')
ylabel('velocity [m/s]')
%legend('AV','HVs')
title('weak string stable for i in {1,2,3} with disturbance on 4-th veh')
legend('veh-1(a=100,b=0.1)','veh-2(a=100,b=0.1)','veh-3(a=20,b=0.5)','veh-4')
ylim([8.2 9.5])
xlim([0 t_f])
end

size_v = size(velocity_AV);
energy = zeros(size_v(1)-1,size_v(2));
Energy = zeros(1,N);

for veh=1:N
    for instant=1:length(energy)
        energy(instant,veh) = (velocity_AV(instant,veh)-v_eq)^2*(t(instant+1)-t(instant));
    end
    Energy(1,veh) = sum(energy(:,veh));
end

y = fliplr(Energy);

figure()
l = categorical({'veh-4','veh-3','veh-2','veh-1'});
x = reordercats(l,{'veh-4','veh-3','veh-2','veh-1'});
bar(x,y)
ylabel('Energy [m^2/s^2]')
xlabel('Vehicles')
title('Energy of vehicles wrt v-equilibrium in weak string stable platoon for i in {1,2,3}')
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

a1=20;
b1=0.5;
% calculation of critical N on a string
a1_bar = a1/(h_eq*h_eq);
b1_bar = b1;
P1_num = [a1_bar b1_bar*k_bar];
P1_den = [1 (a1_bar+b1_bar) b1_bar*k_bar];
PP1 = tf(P1_num,P1_den);

figure()
sigma(PP_d,{10^-1,10})
hold on
sigma(((PP_d)*PP),{10^-1,10})
hold on
sigma(((PP_d)*PP*PP1),{10^-1,10})
hold on
sigma(((PP_d)*PP*PP1*PP1),{10^-1,10})
legend('P_4','P_3(a=20,b=0.5)','P_2(a=20,b=0.5)','P_1(a=20,b=0.5)')
title('strong and weak string unstable')
grid on

%% Non Linear simulation N cars with proposed AV and Disturbance on the AV 
a=20;
b=0.5;
a1=20;
b1=0.5;
v_max = 9.751;
l_v = 4.5;
d_s = 6;
N = 4; %from here I change the numbers of vehicles into the Ring
Ring_length = 260*N/22;
spacing = Ring_length/N;
V = @(x) v_max*((tanh(x-l_v-d_s) + tanh(l_v + d_s))/(1 + tanh(l_v + d_s))); 
Set_of_initial_conditions = v_eq; 
t_f=40; 

for z=1:length(Set_of_initial_conditions)
    
opts = odeset('MaxStep',1e-2);    
y0 = Initial_velocity_and_space_conditions(Set_of_initial_conditions(z),N,spacing);

[t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_String_with_dist_Review(t,y,N,a,b,a1,b1,Ring_length,V),[0 t_f],y0,opts);
[displacement_AV, velocity_AV] = Plot_Displacement_velocity(t,y_N_cars,N,Ring_length);
    
% Uncomment to plot also the positions
% figure()
% plot(t,displacement_AV(:,end),'ro',t,displacement_AV(:,1:end-1),'ko')
% xlabel('time [s]')
% ylabel('position [m]')
% xlim([100 140])
% legend('AV','HVs')
% %title({'displacement N cars with AV', ['Initial velocity condition y_0 =' num2str(Set_of_initial_conditions(z))]})

figure()
plot(t,velocity_AV(:,1:end-1),'k',t,velocity_AV(:,end),'r')%,t,velocity_AV(:,4),t,velocity_AV(:,1))
xlabel('time [s]')
ylabel('velocity [m/s]')
%legend('AV','HVs')
title('strong and weak string unstable with disturbance on 4-th veh')
legend('veh-1(a=20,b=0.5)','veh-2(a=20,b=0.5)','veh-3(a=20,b=0.5)','veh-4')
ylim([8.2 9.5])
xlim([0 t_f])
end

size_v = size(velocity_AV);
energy = zeros(size_v(1)-1,size_v(2));
Energy = zeros(1,N);

for veh=1:N
    for instant=1:length(energy)
        energy(instant,veh) = (velocity_AV(instant,veh)-v_eq)^2*(t(instant+1)-t(instant));
    end
    Energy(1,veh) = sum(energy(:,veh));
end

y = fliplr(Energy);

figure()
l = categorical({'veh-4','veh-3','veh-2','veh-1'});
x = reordercats(l,{'veh-4','veh-3','veh-2','veh-1'});
bar(x,y)
ylabel('Energy [m^2/s^2]')
xlabel('Vehicles')
title('Energy of vehicles wrt v-equilibrium in strong and weak string unstable platoon')
%%
a=200;
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

a1=200;
b1=0.5;
% calculation of critical N on a string
a1_bar = a1/(h_eq*h_eq);
b1_bar = b1;
P1_num = [a1_bar b1_bar*k_bar];
P1_den = [1 (a1_bar+b1_bar) b1_bar*k_bar];
PP1 = tf(P1_num,P1_den);

figure()
sigma(PP_d,{10^-1,10})
hold on
sigma(((PP_d)*PP),{10^-1,10})
hold on
sigma(((PP_d)*PP*PP1),{10^-1,10})
hold on
sigma(((PP_d)*PP*PP1*PP1),{10^-1,10})
legend('P_4','P_3(a=200,b=0.5)','P_2(a=200,b=0.5)','P_1(a=200,b=0.5)')
title('strong string stable')
grid on

%% Non Linear simulation N cars with proposed AV and Disturbance on the AV 
a=200;
b=0.5;
a1=200;
b1=0.5;
v_max = 9.751;
l_v = 4.5;
d_s = 6;
N = 4; %from here I change the numbers of vehicles into the Ring
Ring_length = 260*N/22;
spacing = Ring_length/N;
V = @(x) v_max*((tanh(x-l_v-d_s) + tanh(l_v + d_s))/(1 + tanh(l_v + d_s))); 
Set_of_initial_conditions = v_eq; 
t_f=40; 

for z=1:length(Set_of_initial_conditions)
    
opts = odeset('MaxStep',1e-2);    
y0 = Initial_velocity_and_space_conditions(Set_of_initial_conditions(z),N,spacing);

[t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_String_with_dist_Review(t,y,N,a,b,a1,b1,Ring_length,V),[0 t_f],y0,opts);
[displacement_AV, velocity_AV] = Plot_Displacement_velocity(t,y_N_cars,N,Ring_length);
    
% Uncomment to plot also the positions
% figure()
% plot(t,displacement_AV(:,end),'ro',t,displacement_AV(:,1:end-1),'ko')
% xlabel('time [s]')
% ylabel('position [m]')
% xlim([100 140])
% legend('AV','HVs')
% %title({'displacement N cars with AV', ['Initial velocity condition y_0 =' num2str(Set_of_initial_conditions(z))]})

figure()
plot(t,velocity_AV(:,1:end-1),'k',t,velocity_AV(:,end),'r')%,t,velocity_AV(:,4),t,velocity_AV(:,1))
xlabel('time [s]')
ylabel('velocity [m/s]')
%legend('AV','HVs')
title('strong string stable with disturbance on 4-th veh')
legend('veh-1(a=200,b=0.5)','veh-2(a=200,b=0.5)','veh-3(a=200,b=0.5)','veh-4')
ylim([8.2 9.5])
xlim([0 t_f])
end

size_v = size(velocity_AV);
energy = zeros(size_v(1)-1,size_v(2));
Energy = zeros(1,N);

for veh=1:N
    for instant=1:length(energy)
        energy(instant,veh) = (velocity_AV(instant,veh)-v_eq)^2*(t(instant+1)-t(instant));
    end
    Energy(1,veh) = sum(energy(:,veh));
end

y = fliplr(Energy);

figure()
l = categorical({'veh-4','veh-3','veh-2','veh-1'});
x = reordercats(l,{'veh-4','veh-3','veh-2','veh-1'});
bar(x,y)
ylabel('Energy [m^2/s^2]')
xlabel('Vehicles')
title('Energy of vehicles wrt v-equilibrium in strong string stable platoon')

%%
load("Fig_5_stability_in_the_a-b_plane.mat")

figure()
imagesc(b_bar,a_bar*h_eq*h_eq,Stability_matrix)
hold on
plot(b_bar_bound,a_bar_bound,'w-o','MarkerSize',1.2)
hold on
plot(0.5,200,'kx','MarkerSize',20)
hold on
plot(0.1,100,'bx','MarkerSize',20)
hold on
plot(0.5,20,'rx','MarkerSize',20)
caxis([1 9])
colorbar('Ticks',[1,2,3,4,5,6,7,8,9],...
         'TickLabels',{'3-vehicles','5-vehicles','10-vehicles',...
         '20-vehicles','60-vehicles','100-vehicles','200-vehicles','500-vehicles',...
         'Stability'})
colormap('jet')
legend({'Bound for ||\Gamma||_{\infty}\leq 1'},'TextColor','w','FontSize',19)
legend('boxoff')
xlabel('b')
ylabel('a')
set(gcf, 'Position',  [100, 100, 800, 400])
