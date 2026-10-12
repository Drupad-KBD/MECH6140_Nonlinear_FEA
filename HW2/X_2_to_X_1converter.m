function [converted_X_cap] = X_2_to_X_1converter(x_cap)
% function to convert x,y coordinates from x,y format to x;y format

[n,m] = size(x_cap); %getting size to run for loop
converted_X_cap = zeros(n*m, 1);% dummy array 
a = 1;

for i = 1 : n
    for j = 1 : m
        converted_X_cap(a) = x_cap(i,j);
        a = a + 1;
    end
end

end