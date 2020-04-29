function dydt = ODE_Non_linear_simulation_general_more_communication_v2(t,y,N,a,b,Ring_length,V,alpha,k_veh,N_c,t_a,t_b,t_c)
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

- N_c: is the number of vehicles communicating with the AV

- t_a: time switch to activate AV.

- t_b: time switch to deactivate AV.
- t_c: time switch to reactivate AV.
%}


dydt=zeros(2*N,1);

switch nargin 
    case 7
        for i=1:2:(2*(N-1)-1)
          dydt(i,1)= y(i+1);
          dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
        end
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
    case 10
        for i=1:2:(2*(N-1)-1)
          dydt(i,1)= y(i+1);
          dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
        end
          dydt(2*N-1,1)= y(2*N);
          if N_c > (N-1)
              error('Number of communication required is higher than the total number of vehicles');
          end
          v = 0;
          for i=2:2:(2*N_c)
              v = v + y(i);
          end
          dydt(2*N,1) = k_veh*((1-alpha)*y(2*N_c)+alpha*(y(2*N)+v)/(N_c+1) - y(2*N));
    case 11
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
                if N_c > (N-1)
                    error('Number of communication required is higher than the total number of vehicles');
                end
                v = 0;
                for i=2:2:(2*N_c)
                v = v + y(i);
                end
                dydt(2*N,1) = k_veh*((1-alpha)*y(2*N_c)+alpha*(y(2*N)+v)/(N_c+1) - y(2*N));
        end
    case 12
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
                if N_c > (N-1)
                    error('Number of communication required is higher than the total number of vehicles');
                end
                v = 0;
                for i=2:2:(2*N_c)
                v = v + y(i);
                end
                dydt(2*N,1) = k_veh*((1-alpha)*y(2*N_c)+alpha*(y(2*N)+v)/(N_c+1) - y(2*N));
         end
         if t>=t_b
             for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
         end
    case 13
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
                if N_c > (N-1)
                    error('Number of communication required is higher than the total number of vehicles');
                end
                v = 0;
                for i=2:2:(2*N_c)
                v = v + y(i);
                end
                dydt(2*N,1) = k_veh*((1-alpha)*y(2*N_c)+alpha*(y(2*N)+v)/(N_c+1) - y(2*N));
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
                if N_c > (N-1)
                    error('Number of communication required is higher than the total number of vehicles');
                end
                v = 0;
                for i=2:2:(2*N_c)
                v = v + y(i);
                end
                dydt(2*N,1) = k_veh*((1-alpha)*y(2*N_c)+alpha*(y(2*N)+v)/(N_c+1) - y(2*N));
         end
    otherwise
        error('Wrong number of input arguments');
end