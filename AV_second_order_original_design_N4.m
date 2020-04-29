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
N = 4; % Number of vehicles in the platoon
Ring_length = (260*N)/22;
v_max = 9.751; % Max speed in m/s
h_eq = (Ring_length/N); % equilibrium headway
d_s = 6; % safe distance
V = @(h,l_v,d_s) v_max*((tanh(h-l_v-d_s) + tanh(l_v + d_s))/(1 + tanh(l_v + d_s))); %Velocity non-linear function. h = headway, l_v = Veh length, d_s = safe distance
dV = @(h,l_v,d_s) (v_max/(1+tanh(l_v + d_s)))*(1-(tanh(h - l_v -d_s))^2); % dV/dh derivative of V in the headway
v_eq = V(h_eq, length_standard, d_s); % equilibrium velocity of the ring
k_bar = dV(h_eq, length_standard, d_s); % HV param
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
PP = tf(P_num,P_den);
PP_d = tf(Pd_num,Pd_den);
% Poles of the ring setup for N HVs
[p_HD_ring,z_HD_ring] = pzmap((PP_d)/(1-PP^N));
% Bode of a HV
figure()
bode(PP)
grid on
title('Bode of a HV')
% poles of the ring for N HVs
figure()
plot(p_HD_ring,'kx')
xlabel('Real axis')
ylabel('Imaginary axis')
title('poles of the ring for N HVs')
grid on
%% This analysis shows how steering the ring to a higher velocity can improve from a ring stability pov
% Note, this part does not appear in the paper 
h_eq1=255/21;
k_bar1 = dV(h_eq1, length_standard, d_s); 
v_eq1 = V(h_eq1, length_standard, d_s);
a=20;
b=0.5;
% calculation of critical N on a string
a_bar = a/(h_eq1*h_eq1);
b_bar = b;
% transfer functions of P_N and Gamma_i for the new equilibrium
Pd_num = [1 0];
Pd_den = [1 (a_bar+b_bar) b_bar*k_bar1];
P_num = [a_bar b_bar*k_bar1];
P_den = [1 (a_bar+b_bar) b_bar*k_bar1];
PP1 = tf(P_num,P_den);
PP_d = tf(Pd_num,Pd_den);

figure()
sigma(PP1)
hold on
sigma(PP)
grid on
legend('\Gamma with v_{eq}=v_{AV*}','\Gamma with v_{eq}=v_{*}')
title('bode for HVs transfer function at different equilibrium')

%% Second order (PI with saturation proposed in [16]) 
%HV
alpha_1 = b_bar*k_bar;
alpha_2 = a_bar + b_bar;
alpha_3 = a_bar;

% Compute critical peak and frequency
fpeaK_PP_formula = sqrt(alpha_1*(1-2*(alpha_2/(2*sqrt(alpha_1)))^2));
dumping = alpha_2/(2*sqrt(alpha_1));
gpeak_PP_formula = 1/(2*dumping*sqrt(1-dumping^2));
[gpeak_PP,fpeak_PP] = getPeakGain(PP,1e-5);

%AV
alpha = 0.9;
delta = 23;
k_veh = abs((fpeak_PP*((gpeak_PP)^(1-N)))/(1-alpha/2))*sqrt(1/(1-(gpeak_PP)^(2-2*N)));

beta_1 = k_veh*alpha/delta;
beta_2 = k_veh*(1-alpha/2);
beta_3 = k_veh*(1-alpha/2);
beta_2^2 - beta_3^2 - 2*beta_1

%AV in [16]
AV_N = [beta_3 beta_1];
AV_D = [1 beta_2 beta_1];
AV = tf(AV_N,AV_D);

figure()
sigma(AV)
grid on
title('Bode of PI with saturation proposed in [16]')
%% AV Second Order with c
% c comes from (39)
c=0.5;
% We use an approximation to compute (40) since alpha/delta~0
K_alpha = (c + sqrt(((gpeak_PP)^(2*N-2))*(c^2+fpeak_PP^2) - fpeak_PP^2))/((gpeak_PP)^(2*N-2)-1); 
k_veh = K_alpha/(1-alpha/2);
K_alpha = k_veh*(1-alpha/2);
beta = k_veh*alpha;
beta_1 = beta/delta;
beta_2 = K_alpha;
beta_3 = beta_2;

%Proposed AV Second Order with c
AV_N_c = [beta_3 beta_1];
AV_D_c = [1 beta_2+c beta_1];
AV_c = tf(AV_N_c,AV_D_c);

% Bode
figure()
bode(AV,AV_c)
legend('PI with saturation proposed in [16]','Proposed AV controller design')
%% Check on the PI with saturation proposed in [16] for classical stability
%Checks AV
figure()
bode(AV,PP)
legend('PI with saturation in [16]','HV')
% Theorem 2
figure()
bode(PP^(N-1)*AV)
title('resulting bode for N-1 HVs and the PI with saturation in [16]')

[p_AV_ring,z_AV_ring] = pzmap((PP_d)/(1-AV*PP^(N-1)));
[p_HV_ring,z_HV_ring] = pzmap((PP_d)/(1-PP^(N)));

figure()
plot(p_HV_ring,'kx')
hold on
plot(p_AV_ring,'rd')
grid on
xlabel('Real axis')
ylabel('Imaginary axis')
xlim([-1 0.5])
ylim([-1.5 1.5])
legend('poles of F^{(4)}_i','poles of F^{AV(4)}_i with AV in [16]')

%% %% Check on the Proposed AV for classical stability
% Checks AV_c
figure()
bode(AV_c,PP)
legend('Proposed AV','HV')
% Theorem 2
figure()
bode(PP^(N-1)*AV_c)
title('resulting bode for N-1 HVs and the proposed AV')

[p_AV_c_ring,z_AV_c_ring] = pzmap((PP_d)/(1-AV_c*PP^(N-1)));
[p_HV_ring,z_HV_ring] = pzmap((PP_d)/(1-PP^(N)));

figure()
plot(p_HV_c_ring,'kx')
hold on
plot(p_AV_ring,'rd')
grid on
xlabel('Real axis')
ylabel('Imaginary axis')
xlim([-1 0.5])
ylim([-1.5 1.5])
legend('poles of F^{(4)}_i','poles of F^{AV(4)}_i with proposed AV')
%% Non Linear simulation N cars with proposed AV and Disturbance on the AV 
a=20;
b=0.5;
v_max = 9.751;
l_v = 4.5;
d_s = 6;
N = 4; %from here I change the numbers of vehicles into the Ring
Ring_length = 260*N/22;
spacing = Ring_length/N;

k_veh2 =15;
V = @(x) v_max*((tanh(x-l_v-d_s) + tanh(l_v + d_s))/(1 + tanh(l_v + d_s)));
        
Set_of_initial_conditions = 0; 
t_f=400; 

for z=1:length(Set_of_initial_conditions)
    
opts = odeset('MaxStep',1e-2);    
y0 = Initial_velocity_and_space_conditions(Set_of_initial_conditions(z),N,spacing);

[t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_PrePrint_white_noise(t,y,N,a,b,Ring_length,V,k_veh2,2),[0 t_f],y0,opts);
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
title({'velocity N cars with proposed AV',['Initial velocity condition y_0 =' num2str(Set_of_initial_conditions(z))]})
ylim([5 14])
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
l = categorical({'AV','veh-3','veh-2','veh-1'});
x = reordercats(l,{'AV','veh-3','veh-2','veh-1'});
bar(x,y)
ylabel('Energy [m^2/s^2]')
xlabel('Vehicles')
title('Energy of the vehicles wrt equilibrium in sim with disturbance on AV')

figure()
plot(t,velocity_AV(:,1:end-1),'k',t,velocity_AV(:,end),'r')%,t,velocity_AV(:,4),t,velocity_AV(:,1))
xlabel('time [s]')
ylabel('velocity [m/s]')
hold on
ylim([0 10])
title('simulation on the ring with the proposed AV for N=22 and disturbance on AV')


%% Uncomment this part if you want make a video of the simulation
% Advise: Reduce the simulation time or increase the number of frame per
% sec

% [x,y] = Ring_coordinates(Ring_length,N,displacement_AV);
% [a_x,a_y,b_x,b_y,c_x,c_y,d_x,d_y] = Draw_rectangle(Ring_length,N,displacement_AV);
% k=1;
% curve=animatedline;
% 
% for i=1:25:max(size(x))  
% % if mod(i,10)==0
% %    fprintf('Current time = %g/%g\n', i, max(size(x)));
% % end 
% %figure(1)
% %plot(x(i,1:N-1),y(i,1:N-1),'x',x(i,N),y(i,N),'d')
% % 
% %     subplot(2,1,1);
% %     plot(t(i),velocity(i,N),'r--o')
% %     hold on
% %     ylim([0 10])
% %     xlim([0 200])
% 
% figure(1)
% plot([a_x(i,1) b_x(i,1) d_x(i,1) c_x(i,1) a_x(i,1)],...
%     [a_y(i,1) b_y(i,1) d_y(i,1) c_y(i,1) a_y(i,1)],'b-',...
%     [a_x(i,2) b_x(i,2) d_x(i,2) c_x(i,2) a_x(i,2)],...
%     [a_y(i,2) b_y(i,2) d_y(i,2) c_y(i,2) a_y(i,2)],'b-',...
%     [a_x(i,3) b_x(i,3) d_x(i,3) c_x(i,3) a_x(i,3)],...
%     [a_y(i,3) b_y(i,3) d_y(i,3) c_y(i,3) a_y(i,3)],'b-',...
%     [a_x(i,4) b_x(i,4) d_x(i,4) c_x(i,4) a_x(i,4)],...
%     [a_y(i,4) b_y(i,4) d_y(i,4) c_y(i,4) a_y(i,4)],'b-',...
%     [a_x(i,5) b_x(i,5) d_x(i,5) c_x(i,5) a_x(i,5)],...
%     [a_y(i,5) b_y(i,5) d_y(i,5) c_y(i,5) a_y(i,5)],'b-',...
%     [a_x(i,6) b_x(i,6) d_x(i,6) c_x(i,6) a_x(i,6)],...
%     [a_y(i,6) b_y(i,6) d_y(i,6) c_y(i,6) a_y(i,6)],'b-',...
%     [a_x(i,7) b_x(i,7) d_x(i,7) c_x(i,7) a_x(i,7)],...
%     [a_y(i,7) b_y(i,7) d_y(i,7) c_y(i,7) a_y(i,7)],'b-',...
%     [a_x(i,8) b_x(i,8) d_x(i,8) c_x(i,8) a_x(i,8)],...
%     [a_y(i,8) b_y(i,8) d_y(i,8) c_y(i,8) a_y(i,8)],'b-',...
%     [a_x(i,9) b_x(i,9) d_x(i,9) c_x(i,9) a_x(i,9)],...
%     [a_y(i,9) b_y(i,9) d_y(i,9) c_y(i,9) a_y(i,9)],'b-',...
%     [a_x(i,10) b_x(i,10) d_x(i,10) c_x(i,10) a_x(i,10)],...
%     [a_y(i,10) b_y(i,10) d_y(i,10) c_y(i,10) a_y(i,10)],'b-',...
%     [a_x(i,11) b_x(i,11) d_x(i,11) c_x(i,11) a_x(i,11)],...
%     [a_y(i,11) b_y(i,11) d_y(i,11) c_y(i,11) a_y(i,11)],'b-',...
%     [a_x(i,12) b_x(i,12) d_x(i,12) c_x(i,12) a_x(i,12)],...
%     [a_y(i,12) b_y(i,12) d_y(i,12) c_y(i,12) a_y(i,12)],'b-',...
%     [a_x(i,13) b_x(i,13) d_x(i,13) c_x(i,13) a_x(i,13)],...
%     [a_y(i,13) b_y(i,13) d_y(i,13) c_y(i,13) a_y(i,13)],'b-',...
%     [a_x(i,14) b_x(i,14) d_x(i,14) c_x(i,14) a_x(i,14)],...
%     [a_y(i,14) b_y(i,14) d_y(i,14) c_y(i,14) a_y(i,14)],'b-',...
%     [a_x(i,15) b_x(i,15) d_x(i,15) c_x(i,15) a_x(i,15)],...
%     [a_y(i,15) b_y(i,15) d_y(i,15) c_y(i,15) a_y(i,15)],'b-',...
%     [a_x(i,16) b_x(i,16) d_x(i,16) c_x(i,16) a_x(i,16)],...
%     [a_y(i,16) b_y(i,16) d_y(i,16) c_y(i,16) a_y(i,16)],'b-',...
%     [a_x(i,17) b_x(i,17) d_x(i,17) c_x(i,17) a_x(i,17)],...
%     [a_y(i,17) b_y(i,17) d_y(i,17) c_y(i,17) a_y(i,17)],'b-',...
%     [a_x(i,18) b_x(i,18) d_x(i,18) c_x(i,18) a_x(i,18)],...
%     [a_y(i,18) b_y(i,18) d_y(i,18) c_y(i,18) a_y(i,18)],'b-',...
%     [a_x(i,19) b_x(i,19) d_x(i,19) c_x(i,19) a_x(i,19)],...
%     [a_y(i,19) b_y(i,19) d_y(i,19) c_y(i,19) a_y(i,19)],'b-',...
%     [a_x(i,20) b_x(i,20) d_x(i,20) c_x(i,20) a_x(i,20)],...
%     [a_y(i,20) b_y(i,20) d_y(i,20) c_y(i,20) a_y(i,20)],'b-',...
%     [a_x(i,21) b_x(i,21) d_x(i,21) c_x(i,21) a_x(i,21)],...
%     [a_y(i,21) b_y(i,21) d_y(i,21) c_y(i,21) a_y(i,21)],'b-',...
%     [a_x(i,N) b_x(i,N) d_x(i,N) c_x(i,N) a_x(i,N)],...
%     [a_y(i,N) b_y(i,N) d_y(i,N) c_y(i,N) a_y(i,N)],'r-')
% xlim([-45 45])
% ylim([-45 45])
% 
% drawnow
% F(k) = getframe(gcf);
% k=k+1;
% end
% 
% video = VideoWriter('Modified_PI_with_disturbance_on_AV.avi', 'Motion JPEG AVI');
% %video.FrameRate = 10;
% open(video)
% writeVideo(video, F);
% close(video)


