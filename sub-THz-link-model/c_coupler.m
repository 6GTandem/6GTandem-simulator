% The Mode parameter can be:
% 
classdef c_coupler < c_component

    properties (Access = public)
        Damping = 0; % in dB
    end
    methods
        function out=c_coupler(SimPars)
            out = out@c_component;
        end
        
        function xout=run(obj,x)
            xout=x*db2mag(obj.Damping);
        end
    end
end