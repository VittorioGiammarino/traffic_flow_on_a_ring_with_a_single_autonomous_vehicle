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
[gpeak,fpeak] = getPeakGain(PP,1e-5);
[mag,phase,wout] = bode(PP_d,fpeak);
N_crit_string = log(1/mag)/log(gpeak)


N_max = 22; %max number of vehicles, j = number of vehs
j_peak = zeros(1,N_max);
N_crit_i = (N_max+10)*ones(1,N_max);
F_crit_i = zeros(1,N_max);
conv_num = 1;
conv_den = 1;
for N = 1:N_max
    
    % transfer function of 1/(1-P^N)
    conv_num = conv(conv_num,P_num);
    conv_den = conv(conv_den,P_den);
    PP_N_den = conv_den-[zeros(1,length(conv_den)-length(conv_num)) conv_num];
    PP_N_den = PP_N_den(1:end-1);  % remove the zero which gets canceled with the numerator of Pd
   
    % transfer function of Pd/(1-P^N)
    Total1_num = conv_den;
    Total1_den = conv(PP_N_den,Pd_den);
    Total1_N = tf(Total1_num,Total1_den); %(Pd_den)^j/((PP_den^j-PP_num^j)*D)
    [gpeak_tot1,fpeak_tot1] = getPeakGain(Total1_N,1e-5);
        
    Total2_num = Total1_num;
    Total2_den = Total1_den;
    for i = 1:N
    
         % transfer function of Pd*P^i/(1-P^N)
        Total2_num = conv(Total2_num,P_num);
        Total2_den = conv(Total2_den,P_den);
        Total2_N = tf(Total2_num,Total2_den); %(Pd_den)^j/((PP_den^j-PP_num^j)*D)
    
        [gpeak_tot2,fpeak_tot2] = getPeakGain(Total2_N,1e-5);
    
        if (N_crit_i(N) == (N_max+10)) && (N>1) && (gpeak_tot2 > gpeak_tot1)
            N_crit_i(N) = i;
            F_crit_i(N) = fpeak_tot2;
    %         break
        end
    end
    
%     % get the peak frequency of P and compare it with Pd/(1-P^N) at the same frequency
%     % (here we want to check uniformity of Pd*P^i/(1-P^N) with respect to i)
%     [gpeak,fpeak] = getPeakGain(PP,1e-5);
%     [gpeak_tot1,fpeak_tot1] = getPeakGain(Total1_N,1e-5);
%     [mag,phase,wout] = bode(Total1_N,fpeak);
%     N_crit_i(j) = log(gpeak_tot1/mag)/log(gpeak);

    % transfer function of Pd*P^N/(1-P^N)  (and its peak)
    Total_num = conv_num; 
    Total_den = conv(PP_N_den,Pd_den);
    Total_N = tf(Total_num,Total_den);
    [gpeak,fpeak] = getPeakGain(Total_N,1e-5);
    j_peak(N) = gpeak;
    
    % check if the peak is larger than before
    % (here we want to check uniformity of Pd*P^N/(1-P^N) with respect to N)
    if (N>1) && (j_peak(N) > j_peak(N-1))
        N_crit_ring = N
%         break
    end
    
    if j == 22
        conv_num = 1;
        conv_den = 1;
        for k=1:j
            conv_num = conv(conv_num,P_num);
            conv_den = conv(conv_den,P_den);
            
            Total_2_num = conv(Total1_num,conv_num);
            Total_2_den = conv(Total1_den,conv_den);
            Total_2_tf = tf(Total_2_num,Total_2_den);
        
            if (k==1)||(k==6)||(k==7)
            figure(2)
            bode(Total_2_tf)
            hold on
            end
    
        end
    
    end
end

figure
plot(1:N_max,j_peak)
xlabel('vehicles')
ylabel('Peak')
% This plot says that till N = 15 the peak of Pd*P^N/(1-P^N) will decrease,
% then we will see an increase at increasing N


figure
plot(1:N_max,N_crit_i,'b',1:N_max,1:N_max,'k:')
xlabel('vehicles')
ylabel('i critic')
% This plot says that till N = 6 the peak of Pd*/(1-P^N) is already greater
% than 1, so ring instability (increasing gain) appears already at Pd*P/(1-P^N)
% For N=7, I need five cars to see increasing gain, i.e. Pd*P^5/(1-P^N)
% For N>7, I need a number of cars that is larger than N itself before
% seeing increasing gain. This is the reason why for N=22 we do not see
% string instability, because we need 215 cars to see it (why is impossible)

figure
plot(1:N_max,F_crit_i)
xlabel('vehicles')
ylabel('f critic')
%%
for i=22:23
figure(1)
bode(PP_d/(1-PP^i))
hold on
end

figure(1)
pzmap((PP_d)/(1-PP^30))
hold on


figure(3)
sigma((PP_d)/(1-PP^22))
hold on
sigma((PP_d)/(1-PP^30))
hold on
sigma((PP_d)/(1-PP^40))
hold on
legend('(P_{d})/(1-\Gamma^{22})','(P_{d})/(1-\Gamma^{30})','(P_{d})/(1-\Gamma^{40})')


figure(4)
bode((PP_d*PP)/(1-PP^22))
hold on
bode((PP_d*PP^2)/(1-PP^22))
hold on
bode((PP_d*PP^4)/(1-PP^22))
hold on
bode((PP_d*PP^6)/(1-PP^22))
legend('(P_{N}*\Gamma)/(1-\Gamma^{22})','(P_{N}*\Gamma^{2})/(1-\Gamma^{22})','(P_{N}*\Gamma^{4})/(1-\Gamma^{22})','(P_{N}*\Gamma^{6})/(1-\Gamma^{22})')



N=22
figure(3)
sigma(PP_d/(1-PP^7))
hold on
sigma(PP_d/(1-PP^22))
legend('P_{d}/(1-\Gamma^7)','P_{d}/(1-\Gamma^{22})')

figure()
sigma(1/(1-PP^50))
legend('1/(1-\Gamma^N)')

figure()
sigma(PP)
legend('\Gamma')
grid on

figure()
sigma(PP_d)
legend('P_{d}')

figure()
sigma(PP_d*PP^5/(1-PP^N))
legend('F_{1}^{(7)}')


[gpeak,fpeak] = getPeakGain((PP_d*PP^7)/(1-PP^N),1e-5)
[gpeak,fpeak] = getPeakGain((PP_d*PP^4)/(1-PP^N),1e-5)
[gpeak,fpeak] = getPeakGain(PP^30,1e-5)
[gpeak,fpeak] = getPeakGain((PP_d)/(1-PP^N),1e-5)
%% Ring inStability with AV

a=20;
b=0.5;
k_veh = 0.001;
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
bode((PP_d*PP^6)/(1-AV*PP^21))
grid on
legend('(P_{N}*\Gamma)/(1-\Gamma_{AV}*\Gamma^{21})','(P_{N}*\Gamma^{2})/(1-\Gamma_{AV}*\Gamma^{21})','(P_{N}*\Gamma^{4})/(1-\Gamma_{AV}*\Gamma^{21})','(P_{N}*\Gamma^{6})/(1-\Gamma_{AV}*\Gamma^{21})')



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

















