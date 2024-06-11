% Mode can be 'ideal', 'filter', 'static'

classdef c_iqmodem < c_component
    properties (Access = public)
        IQICoef = 0
        IQIFilter = 1
        IQIDelayImbalance = 0
        DCOffset = 0
    end
    methods
        function out=c_iqmodem(SimPars)
            out = out@c_component;
            if nargin==0, SimPars.IQModem=struct; end
            if isfield(SimPars,'IQModem')
                if isfield(SimPars.IQModem,'Mode'), out.Mode=SimPars.IQModem.Mode; end
                if isfield(SimPars.IQModem,'IQICoef'), out.IQICoef=SimPars.IQModem.IQICoef; end
                if isfield(SimPars.IQModem,'IQIDelayImbalance'), out.IQIDelayImbalance=SimPars.IQModem.IQIDelayImbalance; end
                if isfield(SimPars.IQModem,'IQIFilter'), out.IQIFilter=SimPars.IQModem.IQIFilter; end
                if isfield(SimPars.IQModem,'DCOffset'), out.DCOffset=SimPars.IQModem.DCOffset; end
            end
        end

        function yout=run(obj,yin,phasor)
            import Usefulfunctions.delay
            assert(length(yin)==length(phasor));
            switch obj.Mode
                case 'ideal'
                    yout=yin;
                case 'filter'
                    Del=size((obj.IQIFilter-1)/2);
                    XC=filter(obj.IQIFilter,1,[conj(yin);zeros(Del,1)]);
                    XC=XC(Del+1:end);
                    yout = yin+XC+obj.DCOffset;
                case 'static'
                    yout=delay(real(yin),-obj.IQIDelayImbalance/2)+1i*delay(imag(yin),obj.IQIDelayImbalance/2);
                    yout = yout+obj.IQICoef*conj(yout)+obj.DCOffset;
                otherwise
                    error('Wrong parameter Mode of c_iqmodem class');
            end
            yout=yout.*phasor;
        end

    end
end