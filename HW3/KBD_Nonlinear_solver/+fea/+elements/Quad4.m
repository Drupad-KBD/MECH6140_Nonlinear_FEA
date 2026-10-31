classdef Quad4 < fea.elements.ElementBase
    methods
        function obj=Quad4(X,material,thickness)
            obj@fea.elements.ElementBase(X,material,thickness,4,2);
        end
        function [N,dN]=ShapeFunctions(~,xi)
            r=xi(1); s=xi(2);
            N=0.25*[(1-r)*(1-s),(1+r)*(1-s),(1+r)*(1+s),(1-r)*(1+s)];
            dN=0.25*[-(1-s),(1-s),(1+s),-(1+s);-(1-r),-(1+r),(1+r),(1-r)];
        end
        function [N,dN]=ShapeFun_4Noded(obj,xi)
            [N,dN]=obj.ShapeFunctions(xi);
        end
        function dN=Derivative_Of_ShapeFun_4Noded(obj,xi)
            [~,dN]=obj.ShapeFunctions(xi);
        end
    end
end
