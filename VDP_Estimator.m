classdef VDP_Estimator < CDKF
    
    %add an additional parameter to the object
    properties(Nontunable = true)
        mu = 1; %scaling parameter in VDP equations
    end
    
    methods (Access = protected)
        
        %modify the default validation function to check our new property:
        function validatePropertiesImpl(obj)
            coreValidate(obj) %perform validation of core CDKF properties
            %(ensure consistent dimensions etc.)           
            
            if ~isnumeric(obj.mu) || numel(obj.mu) ~= 1
                error('Mu must be a scalar value')
            end
        end
        
        %overwrite the stateFcn and outputFcn for this specific CDKF. Note
        %that is must be vectorised to handle the sigma points (each row is
        %a set of sigma points for that state)
        function x_new = stateFcn(obj,x,u,w)
            %solve continuous Van der Pol equation, then approximate in
            %discrete time.
            dx = x(2,:);
            dx(2,:) = -x(1,:) + obj.mu.*(1 - x(1,:).^2).*x(2,:);
            
            x_new = x + dx*obj.dt + w;
        end
        
        function y = outputFcn(obj,x,u,v)
            y = x(1,:) .*(1+v);
        end
    end
end