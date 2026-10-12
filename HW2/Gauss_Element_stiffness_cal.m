function [K] = Gauss_Element_stiffness_cal(ShapeFun_derivative,x_cap, D, n, thickness)

[Guss_location,Weight] = gauss_location_weight(n);

converted_X_cap = X_2_to_X_1converter(x_cap);
K = 0;
 for i = 1:n      %Looping for number of gauss points 2 in our case                                                                                                                 (2 for 2D 1 in case of 1D)?
     for j = 1:n
           J = Jacobian_cal(ShapeFun_derivative, [Guss_location(i), Guss_location(j)] , converted_X_cap);
           B = B_Matrix(ShapeFun_derivative, [Guss_location(i), Guss_location(j)], converted_X_cap);
          
% performing the summation of K matrix
% Element_Stiffness = K*thickness; % multiplying with thickness of element to get element stiffness
            K = K + B'*D'*B*det(J)*Weight(i)*Weight(j)*thickness;
      end
 end
 
end
