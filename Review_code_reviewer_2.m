%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Source Code for Paper: Traffic flow on a ring with a single autonomous
%%%                        vehicle: an interconnected perspective
%%% Journal: Transactions on Intelligent transportation Systems, 2019
%%% Authors: Vittorio Giammarino, Simone Baldi, Paolo Frasca and Maria
%%%          Laura Delle Monache
%%%
%%% Copyright Vittorio Giammarino, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non Linear simulation N cars with proposed AV and Noise on all the HVs
v_max = 9.751;
l_v = 4.5;
d_s = 6;
N = 22; %from here I change the numbers of vehicles into the Ring
Ring_length = 260*N/22;
spacing = Ring_length/N;

K_alpha = (c + sqrt(((gpeak_PP)^(2*N-2))*(c^2+fpeak_PP^2) - fpeak_PP^2))/((gpeak_PP)^(2*N-2)-1); 
k_veh2 =K_alpha/(1-alpha/2);
beta_c = k_veh*alpha;
V = @(x) v_max*((tanh(x-l_v-d_s) + tanh(l_v + d_s))/(1 + tanh(l_v + d_s)));
        
Set_of_initial_conditions = 0; 
t_f=400; 

for z=1:length(Set_of_initial_conditions)
    
opts = odeset('MaxStep',1e-2);    
y0 = Initial_velocity_and_space_conditions(Set_of_initial_conditions(z),N,spacing);

[t,y_N_cars] = ode45(@(t,y) ODE_Non_linear_simulation_general_PrePrint_white_noise_Review(t,y,N,a,b,Ring_length,V,k_veh2,2),[0 t_f],y0,opts);
[displacement_AV, velocity_AV] = Plot_Displacement_velocity(t,y_N_cars,N,Ring_length);
    
% Uncomment to plot also the positions
figure()
plot(t,displacement_AV(:,end),'ro',t,displacement_AV(:,1:end-1),'ko')
xlabel('time [s]')
ylabel('position [m]')
xlim([100 140])
legend('AV','HVs')
%title({'displacement N cars with AV', ['Initial velocity condition y_0 =' num2str(Set_of_initial_conditions(z))]})

figure()
plot(t,velocity_AV(:,1:end-1),'k',t,velocity_AV(:,end),'r')%,t,velocity_AV(:,4),t,velocity_AV(:,1))
xlabel('time [s]')
ylabel('velocity [m/s]')
%title({'velocity N cars with AV and noise on a sigle HV',['Initial velocity condition y_0 =' num2str(Set_of_initial_conditions(z))]})
ylim([5 10])
xlim([0 400])
set(gcf,'Position',[100 100 268*5 85*5])
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
l = categorical({'AV','veh-21','veh-20','veh-19','veh-18','veh-17','veh-16','veh-15','veh-14','veh-13','veh-12','veh-11','veh-10','veh-9','veh-8','veh-7','veh-6','veh-5','veh-4','veh-3','veh-2','veh-1'});
x = reordercats(l,{'AV','veh-21','veh-20','veh-19','veh-18','veh-17','veh-16','veh-15','veh-14','veh-13','veh-12','veh-11','veh-10','veh-9','veh-8','veh-7','veh-6','veh-5','veh-4','veh-3','veh-2','veh-1'});
bar(x,y)
ylabel('Energy [m^2/s^2]')
xlabel('Vehicles')
title('Energy of the vehicles wrt equilibrium in simulation with noise')

figure()
plot(t,velocity_AV(:,1:end-1),'k',t,velocity_AV(:,end),'r')
xlabel('time [s]')
ylabel('velocity [m/s]')
title('simulation on the ring with the proposed AV for N=22 and white noise on AV')
ylim([0 10])