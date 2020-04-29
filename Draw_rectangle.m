function [a_x,a_y,b_x,b_y,c_x,c_y,d_x,d_y] = Draw_rectangle(Ring_length,N,position)

time = max(size(position)); 
y = zeros(time,N);
x = zeros(time,N);
P1_x = zeros(time,N);
P1_y = zeros(time,N);
P2_x = zeros(time,N);
P2_y = zeros(time,N);
a_x = zeros(time,N);
a_y = zeros(time,N);
b_x = zeros(time,N);
b_y = zeros(time,N);
c_x = zeros(time,N);
c_y = zeros(time,N);
d_x = zeros(time,N);
d_y = zeros(time,N);
radius = Ring_length/(2*pi);


for i=1:N
    for j=1:max(size(position)) 
        x(j,i) = radius*cos(2*pi*position(j,i)/Ring_length);
        y(j,i) = radius*sin(2*pi*position(j,i)/Ring_length);
        P1_x(j,i) = 0*(cos(pi/2 + 2*pi*position(j,i)/Ring_length)) + x(j,i);
        P1_y(j,i) = 0*(sin(pi/2 + 2*pi*position(j,i)/Ring_length)) + y(j,i);
        P2_x(j,i) = 2.5*(-cos(pi/2 + 2*pi*position(j,i)/Ring_length)) + x(j,i);
        P2_y(j,i) = 2.5*(-sin(pi/2 + 2*pi*position(j,i)/Ring_length)) + y(j,i);
        a_x(j,i) = cos(2*pi*position(j,i)/Ring_length) + P1_x(j,i);
        a_y(j,i) = sin(2*pi*position(j,i)/Ring_length) + P1_y(j,i);
        b_x(j,i) = -cos(2*pi*position(j,i)/Ring_length) + P1_x(j,i);
        b_y(j,i) = -sin(2*pi*position(j,i)/Ring_length) + P1_y(j,i);
        c_x(j,i) = cos(2*pi*position(j,i)/Ring_length) + P2_x(j,i);
        c_y(j,i) = sin(2*pi*position(j,i)/Ring_length) + P2_y(j,i);
        d_x(j,i) = -cos(2*pi*position(j,i)/Ring_length) + P2_x(j,i);
        d_y(j,i) = -sin(2*pi*position(j,i)/Ring_length) + P2_y(j,i);
    end    
end
end