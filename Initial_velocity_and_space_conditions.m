function y0 = Initial_velocity_and_space_conditions(v0,N,spacing)

y0 = zeros(2*N,1);
for i=2:2*N
    if mod(i,2)==1
    y0(i,1)= y0(i-2,1) + spacing;
    else
        if mod(i,2)==0
            y0(i,1)=v0;
        end
    end
end

