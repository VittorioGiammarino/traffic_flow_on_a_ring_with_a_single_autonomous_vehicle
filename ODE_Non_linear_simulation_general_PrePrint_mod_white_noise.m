function dydt = ODE_Non_linear_simulation_general_PrePrint_mod_white_noise(t,y,N,a,b,Ring_length,V,beta,gamma,t_a,t_b,t_c)
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

- beta: variable of the AV.
- gamma: variable of the AV.

- t_a: time switch to activate AV.

- t_b: time switch to deactivate AV.
- t_c: time switch to reactivate AV.
%}

lambda=0.16;
dydt=zeros(2*N+1,1);

v_eq = V(Ring_length/N);

switch nargin 
    case 7 %Only human cars
        for i=1:2:(2*(N-1)-1)
          dydt(i,1)= y(i+1);
          dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
        end
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
    case 9
            for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                Vd=y(2*N+1);
                v_target = Vd + 1*min(max((Ring_length+y(1)-y(2*N-1)-7)/23,0),1);
                Delta_x_s = max(2*(y(2) - y(2*N)),4);
                alpha = min(max(((Ring_length+y(1)-y(2*N-1)-Delta_x_s)/gamma),0),1);
                dydt(2*N,1) = beta*(alpha*v_target+(1-alpha)*y(2)-y(2*N));
                dydt(2*N+1,1) = -lambda*y(2*N+1) + lambda*y(2*N); 
    case 10 %N-1 human cars, 1 AV, t_a is the the activation time for the AV
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
                Vd=v_eq;
                v_target = Vd + 1*min(max((Ring_length+y(1)-y(2*N-1)-7)/23,0),1);
                Delta_x_s = max(2*(y(2) - y(2*N)),4);
                alpha = min(max(((Ring_length+y(1)-y(2*N-1)-Delta_x_s)/gamma),0),1);
                dydt(2*N,1) = beta*(alpha*v_target+(1-alpha)*y(2)-y(2*N));
        end
    case 11   % as above with t_a activation and t_b deactivation for AV   
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
                Vd=v_eq;
                v_target = Vd + 1*min(max((Ring_length+y(1)-y(2*N-1)-7)/23,0),1);
                Delta_x_s = max(2*(y(2) - y(2*N)),4);
                alpha = min(max(((Ring_length+y(1)-y(2*N-1)-Delta_x_s)/gamma),0),1);
                dydt(2*N,1) = beta*(alpha*v_target+(1-alpha)*y(2)-y(2*N));
         end
         if t>=t_b
             for i=1:2:(2*(N-1)-1)
                dydt(i,1)= y(i+1);
                dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
            end
                dydt(2*N-1,1)= y(2*N);
                dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
         end
    case 12       % as above with t_a activation and t_b deactivation and t_c reactivation 
        if (t_b<=t_a) || (t_c<=t_a) || (t_c<=t_b)
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
                Vd=v_eq;
                v_target = Vd + 1*min(max((Ring_length+y(1)-y(2*N-1)-7)/23,0),1);
                Delta_x_s = max(2*(y(2) - y(2*N)),4);
                alpha = min(max(((Ring_length+y(1)-y(2*N-1)-Delta_x_s)/gamma),0),1);
                dydt(2*N,1) = beta*(alpha*v_target+(1-alpha)*y(2)-y(2*N));
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
                Vd=v_eq;
                v_target = Vd + 1*min(max((Ring_length+y(1)-y(2*N-1)-7)/23,0),1);
                Delta_x_s = max(2*(y(2) - y(2*N)),4);
                alpha = min(max(((Ring_length+y(1)-y(2*N-1)-Delta_x_s)/gamma),0),1);
                dydt(2*N,1) = beta*(alpha*v_target+(1-alpha)*y(2)-y(2*N));
         end
    otherwise
        error('Wrong number of input arguments');
end