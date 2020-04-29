function [displacement, velocity] = Plot_Displacement_velocity(t,y,N,Ring_length)

velocity = zeros(length(t),N);
displacement = zeros(length(t),N);
j=1;
k=1;

for i=1:2*N
    if mod(i,2)==1
    displacement(:,j)= y(:,i);
    j=j+1;
    else
        if mod(i,2)==0
            velocity(:,k)=y(:,i);
            k=k+1;
        end
    end
end

size_d=size(displacement);

for i=1:size_d(2)
    for j=1:size_d(1)
        if displacement(j,i)>Ring_length
            displacement(j,i)=displacement(j,i)-(Ring_length*floor(displacement(j,i)/Ring_length));
        end
    end
end