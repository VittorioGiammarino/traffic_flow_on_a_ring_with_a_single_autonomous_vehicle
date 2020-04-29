function dydt = ODE_Non_linear_simulation_general_new_AV_second_order(t,y,N,a,b,Ring_length,V,K_alpha,c,beta,delta,t_a,t_b,t_c)
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
- beta: variable of the AV.

- t_a: time switch to activate AV.

- t_b: time switch to deactivate AV.
- t_c: time switch to reactivate AV.
%}


dydt=zeros(2*N,1);

h_eq = Ring_length/N;
v_eq = V(h_eq);

switch nargin 
    case 7
        for i=1:2:(2*(N-1)-1)
          dydt(i,1)= y(i+1);
          dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
        end
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
    case 11
        for i=1:2:(2*(N-1)-1)
           if (t>=1620) && (t<=1621)
          dydt(i,1)= y(i+1);
          dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1)) - 2;
          else
          dydt(i,1)= y(i+1);
          dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
          end
        end
        if (t>=40) && (t<=41)
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = K_alpha*(y(2)-y(2*N))+(beta/delta)*(Ring_length+y(1)-y(2*N-1)-7)/delta + c*(v_eq - y(2*N)) - 0.05;% - b*y(2*N);
        else
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = K_alpha*(y(2)-y(2*N))+(beta/delta)*(Ring_length+y(1)-y(2*N-1)-7)/delta + c*(v_eq - y(2*N)); %k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 +(alpha*beta)*(Ring_length+y(1)-y(2*N-1)-7)/delta-y(2*N)) + c*(v_eq - y(2*N));% - b*y(2*N);
        end
    case 12
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
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 +(alpha*beta)*(y(1)-y(2*N-1)-7)/23-y(2*N));
        end
    case 13
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
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 +(alpha*beta)*(y(1)-y(2*N-1)-7)/23-y(2*N));
         end
         if t>=t_b
             for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
         end
    case 14
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
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 +(alpha*beta)*(y(1)-y(2*N-1)-7)/23-y(2*N));
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
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 +(alpha*beta)*(y(1)-y(2*N-1)-7)/23-y(2*N));
         end
    otherwise
        error('Wrong number of input arguments');
end