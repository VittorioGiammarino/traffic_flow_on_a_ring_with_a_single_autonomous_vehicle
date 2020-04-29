T = [-1 0 1 0 0 0 0 0;...
    0 -1 0 1 0 0 0 0;...
    0 0 -1 0 1 0 0 0;...
    0 0 0 -1 0 1 0 0;...
    0 0 0 0 -1 0 1 0;...
    0 0 0 0 0 -1 0 1];

T_1 = pinv(T);
%% Initialization
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

v_eq = V(h_eq, length_standard, d_s); % from here we can assume v_eq=9.7 m/s
k_bar = dV(h_eq, length_standard, d_s); 

b_bar = 0.5;
a_bar = 20/(h_eq)^2;
%% Linearization 1

%a_bar = 0.14;
b_bar = 0.5;

A = zeros(2*N,2*N);

% \chi=[x, v]^T 
% for i=2:2:(2*N-2)
%     A(i-1,i) = 1;
%     A(i,i-1) = -b_bar*k_bar;
%     A(i,i) = -a_bar - b_bar;
%     A(i,i+1) = b_bar*k_bar;
%     A(i,i+2) = a_bar;
% end
% 
% A(2*N-1,2*N) = 1;
% A(2*N,1) = b_bar*k_bar;
% A(2*N,2) = a_bar;
% A(2*N,2*N-1) = -b_bar*k_bar;
% A(2*N,2*N) = -a_bar - b_bar;

% \chi=[h, v]^T  modified
for i=2:2:(2*N-2)
    A(i-1,i) = -1;
    A(i-1,i+2) = 1;
    A(i,i-1) = +b_bar*k_bar;
    A(i,i) = -a_bar - b_bar;
    A(i,i+2) = a_bar;
end

A(2*N-1,2*N) = -1;
A(2*N-1,2) = 1;
A(2*N,2) = a_bar;
A(2*N,2*N-1) = 0;
A(2*N,2*N) = -a_bar - b_bar;

for i =1:(2*N-2)
    
    if mod(i,2)==1
        A(end,i)=-b_bar*k_bar;
    end
    
end

[V_l,D_l,W_l] = eig(A);

rank(A)

%%
% \chi=[h, v]^T not modified

A = zeros(2*N,2*N);

for i=2:2:(2*N-2)
    A(i-1,i) = -1;
    A(i-1,i+2) = 1;
    A(i,i-1) = +b_bar*k_bar;
    A(i,i) = -a_bar - b_bar;
    A(i,i+2) = a_bar;
end

A(2*N-1,2*N) = -1;
A(2*N-1,2) = 1;
A(2*N,2) = a_bar;
A(2*N,2*N-1) = +b_bar*k_bar;
A(2*N,2*N) = -a_bar - b_bar;

[V,D,W] = eig(A);

rank(A)

%% Linearization 1

A_N = zeros(2*N-2,2*N-2);

for i=2:2:(2*N-4)
    A_N(i-1,i) = 1;
    A_N(i,i-1) = -b_bar*k_bar;
    A_N(i,i) = -a_bar - b_bar;
    A_N(i,i+1) = b_bar*k_bar;
    A_N(i,i+2) = a_bar;
end

for i =1:(2*N-2)
    
    if mod(i,2)==1
        A_N(end,i)=-b_bar*k_bar;
    end
    if mod(i,2)==0
        A_N(end,i)=-a_bar;
    end
    
end

A_N(2*N-3,2*N-2) = 1;
A_N(2*N-2,2*N-3) = -2*b_bar*k_bar;
A_N(2*N-2,2*N-2) = -2*a_bar - b_bar;


[V_N,D_N,W_N] = eig(A_N);
rank(A_N)



%% Projection for \chi=[x, v]^T 

T = zeros(2*N-2,2*N);
for i=1:2*N-2
    
T(i,i)=-1;
T(i,i+2)=1;

end

T_1=pinv(T);

A_new = T*A*T_1;

%% Projection for \chi=[h, v]^T 

T = zeros(2*N-2,2*N);
for i=1:2:2*N-2
T(i,i)=1;
T(i+1,i+1)=-1;
T(i+1,i+3)=1;
end


T_1=pinv(T);

A_new = T*A*T_1;

[V_new,D_new,W_new]=eig(A_new);


