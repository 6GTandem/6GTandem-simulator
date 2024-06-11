% The Mode parameter can be:
% 
classdef c_fiber < c_component

    properties (Access = public)
        Length = 5; % m
        DampingPerMeter = 5; % dB/m
        FS = 15e9;
%        Filter = [1;0.1;-0.01*1i];
        Filter = 1;
    end
    properties (Dependent)
        Delay    % delay in number of samples
        Damping  % in dB
    end
    methods
        function out=c_fiber(SimPars)
            out = out@c_component;
            if nargin==0, SimPars.PA=struct; end
            if isfield(SimPars,'Fiber')
                if isfield(SimPars.Fiber,'Mode'), out.Mode=SimPars.Fiber.Mode; end
                if isfield(SimPars.Fiber,'Length'), out.Length=SimPars.Fiber.Length; end
                if isfield(SimPars.Fiber,'DampingPerMeter'), out.DampingPerMeter=SimPars.Fiber.DampingPerMeter; end
                if isfield(SimPars.Fiber,'Filter'), out.Filter=SimPars.Fiber.Filter; end
            end
        end
        
        function xout=run(obj,x)
            import Usefulfunctions.*
            xout=delay(filter(obj.Filter,1,x),obj.Delay);
            xout=setdbm(xout,getdbm(x)-obj.Damping);
        end

        function D=get.Delay(obj)
            D=obj.Length*1.5/3e8*obj.FS;
        end    
        function Damp=get.Damping(obj)
            Damp=obj.Length*obj.DampingPerMeter;
        end    

    end
end