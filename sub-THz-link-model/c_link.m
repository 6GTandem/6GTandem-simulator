% c_link: A link in 6GTANDEM.
% A link is defined as a length of plastic fiber (PMF), followed by a
% repeater, which consists of an input coupler, an amplifier, and an output coupler.
% Usage:> L = c_link;
%       > y = L.run(x);
classdef c_link < c_component
    properties (Access = public)
        % components
        Amp=c_amplifier;
        Fiber=c_fiber;
        COUP_IN=c_coupler;
        COUP_OUT=c_coupler;
    end

    methods
        function out=c_link
            out = out@c_component;

            out.Amp=c_amplifier;
            out.Fiber=c_fiber;
            out.COUP_IN=c_coupler;
            out.COUP_OUT=c_coupler;
        end

        function yout=run(obj,y)
            z=obj.Fiber.run(y);
            x1=obj.COUP_IN.run(z);
            x=obj.Amp.run(x1);
            yout=obj.COUP_OUT.run(x);
        end
    end
end
