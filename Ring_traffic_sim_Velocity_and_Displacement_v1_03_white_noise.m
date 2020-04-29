%% Initialization

Lengths = [5.22; 5.15; 4.84; 4.87; 5.15; 5.15; 4.86; 4.92; 5.09; 4.86; 4.86; 5.69; 5.21; 5.15; 4.87; 5.15; 4.86; 4.87; 5.15; 5.70; 4.44; 5.15];
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

% v_eq=zeros(N,1);
% k_bar=zeros(N,1);

v_eq = V(h_eq, length_standard, d_s); % from here we can assume v_eq=9.7 m/s
k_bar = dV(h_eq, length_standard, d_s); 

%   v_eq(i)=V(h_eq, Lengths(i), d_s); % from here we can assume v_eq=9.7 m/s
%   k_bar(i)=dV(h_eq, Lengths(i), d_s);

Human_Driver_tf = @(a_bar, b_bar, k_bar) t_f([a_bar b_bar*k_bar],[1 (a_bar+b_bar) b_bar*k_bar]); % Human Driver OV-FTL transfer function
Autonomous_Veh_tf = @(k_veh, alpha) t_f([0 k_veh*(1-(alpha/2))],[1 k_veh*(1-(alpha/2))]); % Autonomous Vehicle OV-FTL transfer function

%% Manual Time selection for AV control

% change t_switch to change the moment in which the AV is activated
% change tf for the total time simulation

v_max = 9.751;
l_v = 4.5;
d_s = 6;

a = 140;
b = 0.1;

N = 22; %from here I change the numbers of vehicles into the Ring
Ring_length = 260*N/22;
spacing = Ring_length/N;

t_f = 200;
t_switch = 0;

V = @(x) v_max*((tanh(x-l_v-d_s) + tanh(l_v + d_s))/(1 + tanh(l_v + d_s)));

opts = odeset('MaxStep',1e-1); 
y0 = Initial_velocity_and_space_conditions(v_eq-1,N,spacing);
%[t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_with_disturbance(t,y,N,a,b,Ring_length,V),[0 t_f],y0,opts);
[t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_String_with_disturbance(t,y,N,a,b,Ring_length,V),[0 t_f],y0,opts);
[displacement, velocity] = Plot_Displacement_velocity(t,y_N_cars,N,Ring_length);

% figure()
% plot(t,displacement, '*')
% xlabel('time [s]')
% ylabel('displacement [m]')
% % title(['displacement N cars with manual AV activation',' t_{switch} = ', num2str(t_switch)])

%v_interest = velocity(:,18:22);

figure()
plot(t,velocity(:,:))%1),t,velocity(:,14),t,velocity(:,22))
xlabel('time [s]')
ylabel('velocity [m/s]')
ylim([3 10])
%legend('veh 1','veh 14', 'veh 22')
%title(['velocity N cars with manual AV activation', ' t_{switch} = ', num2str(t_switch)])

% figure()
% surf(t,displacement,velocity)

%% plot vehicles

[x,y] = Ring_coordinates(Ring_length,N,displacement);

for i=10000:3:max(size(x))
plot(x(i,1:22),y(i,1:22),'x')
F(i) = getframe(gcf) ;
drawnow
end
  % create the video writer with 1 fps
  writerObj = VideoWriter('myVideo.avi');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


%% Velocity and Displacement simulation N cars with second order model AV RING INSTABILITY

v_max = 9.751;
l_v = 4.5;
d_s = 6;
N = 5; %from here I change the numbers of vehicles into the Ring
Ring_length = 260*N/22;
spacing = Ring_length/N;

a = 0.14*h_eq*h_eq;
b = 0.5;
k_veh = 0.5;
%k_veh = 0.00000001;
alpha = 0.9;
delta = 23;
%c = -k_bar*h_eq +v_eq;
beta = 1;%(b*c*delta)/(alpha*7);
V = @(x) v_max*((tanh(x-l_v-d_s) + tanh(l_v + d_s))/(1 + tanh(l_v + d_s)));
        
Set_of_initial_conditions = v_eq-5; %[0 1 2 3 4 5 6 7 8 9 v_eq(1)]; 
t_f=50; % Final simulation time

for z=1:length(Set_of_initial_conditions)
    
opts = odeset('MaxStep',1e-2);    
y0 = Initial_velocity_and_space_conditions(Set_of_initial_conditions(z),N,spacing);
% [t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_second_order_AV(t,y,N,a,b,Ring_length,V,alpha,k_veh,beta,delta),[0 t_f],y0,opts);
% [t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_new_AV_second_order(t,y,N,a,b,Ring_length,V,alpha,k_veh,beta,delta),[0 t_f],y0,opts);
% [t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general(t,y,N,a,b,Ring_length,V,alpha,k_veh),[0 t_f],y0,opts);
% [t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_new_AV_first_order(t,y,N,a,b,Ring_length,V,alpha,k_veh),[0 t_f],y0,opts);
 [t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_String_with_disturbance(t,y,N,a,b,Ring_length,V,alpha,k_veh,beta,delta),[0 t_f],y0,opts);
 
[displacement_AV, velocity_AV] = Plot_Displacement_velocity(t,y_N_cars,N,Ring_length);
    
figure()
plot(t,displacement_AV,'x')
xlabel('time [s]')
ylabel('displacement [m]')
title({'displacement N cars with AV', ['Initial velocity condition y_0 =' num2str(Set_of_initial_conditions(z))]})

figure()
plot(t,velocity_AV(:,:))%,t,velocity_AV(:,4),t,velocity_AV(:,1))
xlabel('time [s]')
ylabel('velocity [m/s]')
grid on
%title({'velocity N cars with AV',['Initial velocity condition y_0 =' num2str(Set_of_initial_conditions(z))]})
ylim([0 10])
%legend('veh 1','veh 2','veh 3', 'AV')
end

size_v = size(velocity_AV);
energy = zeros(size_v(1)-1,size_v(2));
Energy = zeros(1,N);

for veh=1:N
    for instant=1:length(energy)
        energy(instant,veh) = (velocity_AV(instant,veh)-v_eq(1))^2*(t(instant+1)-t(instant));
    end
    Energy(1,veh) = sum(energy(:,veh));
end

y = fliplr(Energy);

% figure()
% c = categorical({'AV','veh-21','veh-20','veh-19','veh-18','veh-17','veh-16','veh-15','veh-14','veh-13','veh-12','veh-11','veh-10','veh-9','veh-8','veh-7','veh-6','veh-5','veh-4','veh-3','veh-2','veh-1'});
% x = reordercats(c,{'AV','veh-21','veh-20','veh-19','veh-18','veh-17','veh-16','veh-15','veh-14','veh-13','veh-12','veh-11','veh-10','veh-9','veh-8','veh-7','veh-6','veh-5','veh-4','veh-3','veh-2','veh-1'});
% bar(x,y)
% ylabel('Energy [m^2/s^2]')
% xlabel('Vehicles')
%% plot vehicles

[x,y] = Ring_coordinates(Ring_length,N,displacement_AV);

for i=1:10:max(size(x))
    
% if mod(i,10)==0
%    fprintf('Current time = %g/%g\n', i, max(size(x)));
% end 
    
figure(1)
plot(x(i,1:N-1),y(i,1:N-1),'x',x(i,N),y(i,N),'rd')
xlim([-1.2 1.2])
ylim([-1.2 1.2])

% F(i) = getframe(gcf) ;
% drawnow
end
  % create the video writer with 1 fps
%   writerObj = VideoWriter('myVideo.avi');
%   writerObj.FrameRate = 500;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

%% Velocity and Displacement simulation N cars with second order model AV RING STABILITY

v_max = 9.751;
l_v = 4.5;
d_s = 6;
N = 22; %from here I change the numbers of vehicles into the Ring
Ring_length = 260*N/22;
spacing = Ring_length/N;

a = 100;%0.14*h_eq*h_eq;
b = 0.5;
k_veh = 1;
alpha = 0.9;
delta = 23;
c = -k_bar*h_eq +v_eq;
beta = 1;%(b*c*delta)/(alpha*7);
V = @(x) v_max*((tanh(x-l_v-d_s) + tanh(l_v + d_s))/(1 + tanh(l_v + d_s)));
        
Set_of_initial_conditions = v_eq-1; %[0 1 2 3 4 5 6 7 8 9 v_eq(1)]; 
t_f=400; % Final simulation time

for z=1:length(Set_of_initial_conditions)
    
opts = odeset('MaxStep',1e-2);    
y0 = Initial_velocity_and_space_conditions(Set_of_initial_conditions(z),N,spacing);
[t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_second_order_AV(t,y,N,a,b,Ring_length,V,alpha,k_veh,beta,delta),[0 t_f],y0,opts);
[displacement_AV, velocity_AV] = Plot_Displacement_velocity(t,y_N_cars,N,Ring_length);
    
    
% figure()
% plot(t,displacement_AV,'x')
% xlabel('time [s]')
% ylabel('displacement [m]')
%title({'displacement N cars with AV', ['Initial velocity condition y_0 =' num2str(Set_of_initial_conditions(z))]})

figure()
plot(t,velocity_AV(:,15:22))
xlabel('time [s]')
ylabel('velocity [m/s]')
%title({'velocity N cars with AV',['Initial velocity condition y_0 =' num2str(Set_of_initial_conditions(z))]})
ylim([0 10])
end

size_v = size(velocity_AV);
energy = zeros(size_v(1)-1,size_v(2));
Energy = zeros(1,N);

for veh=1:N
    for instant=1:length(energy)
        energy(instant,veh) = (velocity_AV(instant,veh)-v_eq(1))^2*(t(instant+1)-t(instant));
    end
    Energy(1,veh) = sum(energy(:,veh));
end

y = fliplr(Energy);

figure()
c = categorical({'AV','veh-21','veh-20','veh-19','veh-18','veh-17','veh-16','veh-15','veh-14','veh-13','veh-12','veh-11','veh-10','veh-9','veh-8','veh-7','veh-6','veh-5','veh-4','veh-3','veh-2','veh-1'});
x = reordercats(c,{'AV','veh-21','veh-20','veh-19','veh-18','veh-17','veh-16','veh-15','veh-14','veh-13','veh-12','veh-11','veh-10','veh-9','veh-8','veh-7','veh-6','veh-5','veh-4','veh-3','veh-2','veh-1'});
bar(x,y)
ylabel('Energy [m^2/s^2]')
xlabel('Vehicles')

