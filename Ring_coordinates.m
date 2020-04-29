function [x,y] = Ring_coordinates(Ring_length,N,position)

time = max(size(position)); 
y = zeros(time,N);
x = zeros(time,N);
radius = Ring_length/(2*pi);

for i=1:N
    for j=1:max(size(position)) 
        x(j,i) = radius*cos(2*pi*position(j,i)/Ring_length);
        y(j,i) = radius*sin(2*pi*position(j,i)/Ring_length);
    end    
end
end