function dydt = ODE_Non_linear_simulation_general_String_with_dist_Review(t,y,N,a,b,a1,b1,Ring_length,V,alpha,k_veh,beta,delta,t_a,t_b,t_c)
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

v_eq = V(Ring_length/N);
% a2=200;
% b2=0.1;

dydt=zeros(2*N,1);

switch nargin 
    case 9 %only human vehicles
        for i=1:2:(2*(N-1)-1)
            dydt(i,1)= y(i+1);
            dydt(i+1,1)= a1*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b1*(V(y(i+2)-y(i))-y(i+1));
          if i == (2*(N-1)-1)
            dydt(i,1)= y(i+1);
            dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
          end
        end
        if (t <= 3) && (t>=1)
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = 0.6*(v_eq-y(2*N)) - 0.5; %a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
        else
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) =0.6*(v_eq-y(2*N)); %a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
        end
        
        
    otherwise
        error('Wrong number of input arguments');
end