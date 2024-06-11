classdef c_dac < c_component
    properties (Access = public)
        Trunclevel=inf      % limits output to +/-Trunclevel
        Nobits=inf          % 
    end
    methods
        function out=c_dac(SimPars)
            out = out@c_component;
            if nargin==0, SimPars.DAC=struct; end
            if isfield(SimPars,'DAC')
                if isfield(SimPars.DAC,'TruncLevel'), out.Trunclevel=SimPars.DAC.TruncLevel; end
                if isfield(SimPars.DAC,'NoBits'), out.Nobits=SimPars.DAC.NoBits; end
            end
        end
        
        function yout=run(obj,yin)
            if obj.Trunclevel>1e98
                yout=yin;
                return;
            end
            if obj.Nobits>20
                yout = Usefulfunctions.limiter(yin,obj.Trunclevel);
            else
                Step=2*obj.TruncLevel/2^obj.Nobits;
                yout = Usefulfunctions.limiter(round(yin/Step)*Step,obj.Trunclevel);
            end
        end
        
    end
end