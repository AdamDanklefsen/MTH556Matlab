classdef BC_rect
    
    % a*u + b*du/dn = L
    properties
        Dirichlet = nan;        % True for Dirichlet, false for Robin/Neumann
        a = nan;                
        L = nan;
        b = nan;
    end
    
    methods
        function obj = BC_rect(a,L,b)
            obj.Dirichlet = true;
            obj.a = a;
            obj.L = L;
            if nargin==3
                obj.Dirichlet = false;
                obj.b = b;
            end
        end
    end
end

