function dydt = ODE_Non_linear_simulation_general_with_disturbance_SIM(t,y,N,a,b,Ring_length,V,alpha,k_veh,t_a,t_b,t_c)
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


dydt=zeros(2*N,1);

switch nargin 
    case 7 %only human vehicles
        for i=1:2:(2*(N-1)-1)
          dydt(i,1)= y(i+1);
          dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
        end
%         if (t <= 1.05) && (t>=0.95)
%         if (t <= 1e-8) && (t>=0)
        if (t == 0)
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N)) + 100;
        else
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
        end
    case 9 %with AV
        for i=1:2:(2*(N-1)-1)
          dydt(i,1)= y(i+1);
          dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
        end
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 - y(2*N));
    case 10
        if t <= t_a
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
        else
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 - y(2*N));
        end
    case 11
        if t_b<=t_a
            error('Final time greater than initial time');
        end
         if t <= t_a
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
         end
         if (t >= t_a)&&(t <= t_b)
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 - y(2*N));
         end
         if t>=t_b
             for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
         end
    case 12
                 if t <= t_a
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
         end
         if (t >= t_a)&&(t <= t_b)
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 - y(2*N));
         end
         if (t >= t_b)&&(t <= t_c)
             for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
         end
         if (t >= t_c)
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 - y(2*N));
         end
    otherwise
        error('Wrong number of input arguments');
end