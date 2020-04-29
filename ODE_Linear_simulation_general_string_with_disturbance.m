 function dydt = ODE_Linear_simulation_general_string_with_disturbance(t,y,N,a_bar,b_bar,Ring_length,V,k_bar,alpha,k_veh,t_a,t_b,t_c)
%{
INPUTS:

- t: time for ODE,
- y: vector(2*N,1). It represents the dynamics of the N vehicles. For each
vehicle i, first is added the variation of displacement dynamics then the
variation of the velocity.
- N: Number of vehicles in the ring.
- a: [m^2/s] Variable of the human driver.
- b: [1/s^(-1)] Variable of the human driver.
- Ring_lenght: Length of the path.
- V: Function for the velocity saturation. See Optimal Velocity model.

- alpha: variable of the AV.
- k_veh: variable of the AV.

- t_a: time switch to activate AV.

- t_b: time switch to deactivate AV.
- t_c: time switch to reactivate AV.
%}

%ye = Equilibrium_conditions(v_eq,N,h_eq);

dydt=zeros(2*N,1);

h_eq = Ring_length/N;
v_eq = V(h_eq);
c = -(k_bar*h_eq)+v_eq;

switch nargin 
    case 8
        for i=1:2:(2*(N-1)-1)
          dydt(i,1)= y(i+1);
          dydt(i+1,1)= a_bar*(y(i+3) - y(i+1)) + b_bar*(k_bar*(y(i+2)-y(i))+c-y(i+1));
        end
% %           if (t <= 1.05) && (t>=0.95)
%         if (t == 0)
% %           dydt(2*N-1,1)= y(2*N) + 0.1;
% %           dydt(2*N,1) = a_bar*(y(2) - y(2*N)) + b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N));
%           dydt(2*N-1,1)= y(2*N);
%           dydt(2*N,1) = 0 + 100;
%         else
% %           dydt(2*N-1,1)= y(2*N);
% %           dydt(2*N,1) = a_bar*(y(2) - y(2*N)) + b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N));
%           dydt(2*N-1,1)= y(2*N);
%           dydt(2*N,1) = 0;
%         end
%           if (t <= 1.05) && (t>=0.95)
        if (t >= 0) && (t < 10)
%           dydt(2*N-1,1)= y(2*N) + 0.1;
          dydt(2*N,1) = a_bar*(y(2) - y(2*N)) + b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N)) + 100;
          dydt(2*N-1,1)= y(2*N);
%           dydt(2*N,1) = 0 + 100;
%         elseif (t >= 10) && (t < 20)
% %           dydt(2*N-1,1)= y(2*N) + 0.1;
%           dydt(2*N,1) = a_bar*(y(2) - y(2*N)) + b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N)) + 200;
%           dydt(2*N-1,1)= y(2*N);
% %           dydt(2*N,1) = 0 + 200;
%         elseif (t >= 20) && (t < 30)
% %           dydt(2*N-1,1)= y(2*N) + 0.1;
%           dydt(2*N,1) = a_bar*(y(2) - y(2*N)) + b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N)) + 100;
%           dydt(2*N-1,1)= y(2*N);
% %           dydt(2*N,1) = 0 + 100;
        else
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = 0; %a_bar*(y(2) - y(2*N)) + b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N));
%          dydt(2*N-1,1)= y(2*N);
%          dydt(2*N,1) = 0;
        end
%           dydt(2*N-1,1)= y(2*N);
%           dydt(2*N,1) = a_bar*(y(2) - y(2*N)) + b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N));
    case 10
        for i=1:2:(2*(N-1)-1)
          dydt(i,1)= y(i+1);
          dydt(i+1,1)= a_bar*(y(i+3) - y(i+1)) + b_bar*(k_bar*(y(i+2)-y(i))+c-y(i+1));
        end
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 - y(2*N));
    case 11
        if t <= t_a
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a_bar*(y(i+3) - y(i+1)) + b_bar*(k_bar*(y(i+2)-y(i))+c-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a_bar*(y(2) - y(2*N)) + b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N));
        else
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a_bar*(y(i+3) - y(i+1)) + b_bar*(k_bar*(y(i+2)-y(i))+c-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 - y(2*N));
        end
    case 12
        if t_b<=t_a
            error('Final time greater than initial time');
        end
         if t <= t_a
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a_bar*(y(i+3) - y(i+1)) + b_bar*(k_bar*(y(i+2)-y(i))+c-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a_bar*(y(2) - y(2*N)) + b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N));
         end
         if (t >= t_a)&&(t <= t_b)
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a_bar*(y(i+3) - y(i+1)) + b_bar*(k_bar*(y(i+2)-y(i))+c-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 - y(2*N));
         end
         if t>=t_b
             for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a_bar*(y(i+3) - y(i+1)) + b_bar*(k_bar*(y(i+2)-y(i))+c-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a_bar*(y(2) - y(2*N))+ b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N));
         end
    case 13
                 if t <= t_a
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a_bar*(y(i+3) - y(i+1)) + b_bar*(k_bar*(y(i+2)-y(i))+c-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a_bar*(y(2) - y(2*N)) + b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N));
         end
         if (t >= t_a)&&(t <= t_b)
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a_bar*(y(i+3) - y(i+1))+ b_bar*(k_bar*(y(i+2)-y(i))+c-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 - y(2*N));
         end
         if (t >= t_b)&&(t <= t_c)
             for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a_bar*(y(i+3) - y(i+1)) + b_bar*(k_bar*(y(i+2)-y(i))+c-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a_bar*(y(2) - y(2*N)) + b_bar*(k_bar*(Ring_length+y(1)-y(2*N-1))+c-y(2*N));
         end
         if (t >= t_c)
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a_bar*(y(i+3) - y(i+1)) + b_bar*(k_bar*(y(i+2)-y(i))+c-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 - y(2*N));
         end
    otherwise
        error('Wrong number of input arguments');
end