% Definition of a transmitter
% Components: c_oscillator, c_iqmod, c_pa, c_dac, c_antenna
% Usage:> TX = c_transmitter;
%       > y = TX.run(x);
classdef c_transmitter < c_component
    properties (Access = public)
        OSC = c_oscillator
        IQMOD = c_iqmodem
        PA = c_amplifier
        DAC = c_dac
        Delay=0;
    end
    
    methods
        function out=c_transmitter(Pars)
            if nargin==0, Pars.TX=struct; end
            out = out@c_component;
            out.DAC=c_dac(Pars.TX);
            out.OSC=c_oscillator(Pars.TX);
            out.IQMOD=c_iqmodem(Pars.TX);
            out.PA=c_amplifier(Pars.TX);
            if isfield(Pars.TX,'Delay'), out.Delay=Pars.TX.Delay; end
        end
        
        function xout=run(obj,x,phasor)
            if nargin<3, phasor=obj.OSC.run(length(x)); end
            x1 = obj.DAC.run(x);
            x2 = obj.IQMOD.run(x1, phasor);
            xout = obj.PA.run(x2);
            if obj.Delay~=0
                xout=Usefulfunctions.delay(xout,obj.Delay);
            end
        end
    end
end