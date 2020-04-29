function y0 = Initial_velocity_and_Random_space_conditions(v0,N,spacing,d_s)

y0 = zeros(2*N,1);
for i=2:2*N
    if mod(i,2)==1
    y0(i,1)= y0(i-2,1) + d_s + rand*spacing;
    else
        if mod(i,2)==0
            y0(i,1)=v0;
        end
    end
end