function dydt = ODE_Non_linear_simulation_general_H_inf_controller_with_noise(t,y,N,a,b,k_bar,Ring_length,V,K,t_a,t_b,t_c)
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
b_bar =b;
a_bar =a/(h_eq*h_eq);

switch nargin 
    case 7
        for i=1:2:(2*(N-1)-1)
          dydt(i,1)= y(i+1);
          dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
        end
          dydt(2*N-1,1)= y(2*N);
          dydt(2*N,1) = a*((y(2) - y(2*N))/(Ring_length+y(1) - y(2*N-1))^2) + b*(V(Ring_length+y(1)-y(2*N-1))-y(2*N));
    case 9
        for i=1:2:(2*(N-1)-1)
           dydt(i,1)= y(i+1);
           dydt(i+1,1)= awgn(a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1)),25);
          if i == (2*(N-1)-1) %- (N-2))
              if (t>10) && (t<11)
                    dydt(i,1)= y(i+1);
                    dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1)) - 0.5;
              else
                    dydt(i,1)= y(i+1);
                    dydt(i+1,1)= a*((y(i+3) - y(i+1))/(y(i+2) - y(i))^2) + b*(V(y(i+2)-y(i))-y(i+1));
              end
          end
        end
        if (t>410) && (t<411) %(t>10) && (t<11)
        position = 0;
        velocity = 0;
        for i=3:(2*N-2)
            if mod(i,2)==1
                position = position + K(i-2)*(y(i)-y(i-2)-h_eq);
            end
            if mod(i,2)==0
                velocity = velocity + K(i-2)*(y(i)-y(i-2));
            end
        end
        AV = (K(2*N-3)+b_bar*k_bar)*(y(2*N-1)-y(2*N-3)-h_eq)+(K(2*N-2)+a_bar)*(y(2*N)-y(2*N-2)) -b_bar*(y(2*N-2)-9.3);
        dydt(2*N-1,1)= y(2*N);
        dydt(2*N,1) = position + velocity + AV-100;%1000;
        else
        position = 0;
        velocity = 0;
        for i=3:(2*N-2)
            if mod(i,2)==1
                position = position + K(i-2)*(y(i)-y(i-2)-h_eq);
            end
            if mod(i,2)==0
                velocity = velocity + K(i-2)*(y(i)-y(i-2));
            end
        end
        AV = (K(2*N-3)+b_bar*k_bar)*(y(2*N-1)-y(2*N-3)-h_eq)+(K(2*N-2)+a_bar)*(y(2*N)-y(2*N-2)) -b_bar*(y(2*N-2)-v_eq);
        dydt(2*N-1,1)= y(2*N);
        dydt(2*N,1) = position + velocity + AV;
        end
  
%         if (t>100) && (t<103) %(t>10) && (t<13) 
%           dydt(2*N-1,1)= y(2*N);
%           dydt(2*N,1) = ((K(1))*(y(3)-y(1)-h_eq)+(K(2))*(y(4)-y(2))+(K(3))*(y(5)-y(3)-h_eq)+...
%               (K(4))*(y(6)-y(4)) + (K(5))*(y(7)-y(5)+6.5) + (K(6))*(y(8)-9.3-y(6))) - 5;% - b*y(2*N);
%         else
%           dydt(2*N-1,1)= y(2*N);
%           dydt(2*N,1) = ((K(1))*(y(3)-y(1)-h_eq)+(K(2))*(y(4)-y(2))+(K(3))*(y(5)-y(3)-h_eq)+...
%               (K(4))*(y(6)-y(4)) + (K(5))*(y(7)-y(5)+6.5) + (K(6))*(y(8)-9.3-y(6)));% - b*y(2*N);
%         end
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
                dydt(2*N,1) = k_veh*((1-alpha)*y(2)+alpha*(y(2*N)+y(2))/2 +(alpha*beta)*(y(1)-y(2*N-1)-7)/23-y(2*N));
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