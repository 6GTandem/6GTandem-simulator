% c_amplifier Amplifier class.
% The c_amplifier implements a nonlinear function with low-signal gain
% 1, max output amplitude 1 (13 dBm).
% Y = G*f(X+W)
% where X is the input signal, W is noise. f() is a nonlinear function with
% max amplitude 1 V. G is the output gain.
% The gain, nonlinearity function and noise power level can be controlled
% Example usage:
% pa=c_pa;
% pa.Mode='tanh';
% pa.Gain=3.2;
% y=pa.run(x);
%
% The Mode parameter can be:
% 'ideal','linear','atan', 'tanh', 'poly3','poly3_pm',
% 'poly5','limiter','softlimiter','weblab','6gtandem'
% Desription: 
classdef c_amplifier < c_component
    properties (Access = public)
        Gain=1; % The low-signal gain of the amplifier
        MaxOutputAmplitude=1;
        NoiseVar = 0; % The variance of AWGN noise [V^2] per channel
        Smoothness = 1; % Used in Mode 'softlimiter'
    end
    methods
        function out=c_amplifier(SimPars)
            out = out@c_component;
            if nargin==0, SimPars.PA=struct; end
            if isfield(SimPars,'PA')
                if isfield(SimPars.PA,'Mode'), out.Mode=SimPars.PA.Mode; disp(out.Mode);end
                if isfield(SimPars.PA,'Smoothness'), out.Smoothness=SimPars.PA.Smoothness; end
                if isfield(SimPars.PA,'Gain'), out.Gain=SimPars.PA.Gain; end
            end
        end
        
        function xout=run(obj,x)
            x=obj.Gain/obj.MaxOutputAmplitude*(x + sqrt(obj.NoiseVar) * (randn(size(x))+1i*randn(size(x))));
            switch obj.Mode
                case {'ideal','linear'}
                    xout=x;
                case 'atan'
                    %alpha=0.6340;% This factor gives 1dB compression at x=1
                    alpha=2/pi;% This factor gives Amax=1;
                    xout = 1/alpha*atan(alpha*abs(x)).*exp(1i*angle(x));
                case 'tanh'
                    %alpha=0.6125;% This factor gives 1dB compression at x=1
                    alpha=1;% This factor gives Amax=1;
                    xout = 1/alpha*tanh(alpha*abs(x)).*exp(1i*angle(x));
                case 'poly3'
                    %alpha=1-10^(-0.05*1);% This factor gives 1dB compression at x=1
                    alpha=4/27; % This factor gives Amax=1;
                    xout = x.*(1-alpha*abs(x).^2);
                    peakx=sqrt(1/3/alpha); % dont allow for negative gain
                    ind=(abs(x)>peakx);
                    xout(ind)=peakx*(1-alpha*abs(peakx).^2)*exp(1i*angle(xout(ind)));
                case 'poly3_pm'
                    %alpha=(1-10^(-0.05*1))*exp(1i*0.2);% This factor gives 1dB compression at x=1
                    alpha=4/27.*exp(1i*0.2);% This factor Amax=1
                    xout = x.*(1-alpha*abs(x).^2);
                    peakx=sqrt(1/3/abs(alpha)); % dont allow for negative gain
                    ind=(abs(x)>peakx);
                    xout(ind)=peakx*(1-alpha*abs(peakx).^2)*exp(1i*angle(xout(ind)));
                case 'poly5'
                    alpha=10/9-sqrt(20/9*10^-0.05-80/81); %This factor gives 1dB compression at x=1
                    beta=9*alpha^2/20; % gives minimum derivative=0, at x=sqrt(2/3/alpha)
                    xout = obj.MaxOutputAmplitude*x.*(1-alpha*abs(x).^2+beta*abs(x).^4);
                case 'limiter'
                    xout=x;
                    ind=abs(xout)>1;  % this gives saturation at 1 V
                    xout(ind)=xout(ind)./abs(xout(ind));
                case 'softlimiter'
                    xout=Usefulfunctions.softlimiter(x,obj.Smoothness);
                case '6gtandem'
                    xout=Usefulfunctions.softlimiter(x,obj.Smoothness);
            end
            xout=obj.MaxOutputAmplitude*xout;
        end

        function setnoisevar(obj,T,B,NFDB)
            k=1.3806E-23;
            NoiseDensity=k*T*db2pow(NFDB);
            NoisePower=NoiseDensity*B;
            % P = U^2/R
            VoltagePower=NoisePower*50;
            NoiseVariancePerChannel=VoltagePower/2;
            obj.NoiseVar=NoiseVariancePerChannel;
        end

        function setmaxpower(obj,maxPdbm)
            obj.MaxOutputAmplitude=sqrt(50*db2pow(maxPdbm)*1e-3);
        end

        function setaveragepower(obj,Pin,Pout)   % will give the desired output power if no nonlinearity, some discrepancy with NL
            obj.Gain=db2mag(Pout-Pin);
        end
    end
end

