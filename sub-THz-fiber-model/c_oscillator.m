% c_oscillator class.
% example usage:
% X=c_oscillator;
% pn=X.run_phase(1000); % generates 1000 phase noise samples
% modes of pn generation 'ideal', 'cfo', 'model', 'spectrum'
classdef c_oscillator < c_component
    properties (Access = public)
        % common parameters for all modes
        FS=1.5e9        % sampling frequency [Hz]
        CFO=0           % carrier frequency offset (only receiver) [Hz]
        CFO_std=0       % carrier frequency random offset (only receiver) [Hz]
        
        % specific param for mode='model'
        L100_db=-200    % pn level at 100 kHz [dB]
        Linf_db=-300    % white pn level [dB]
        F3db=0          % cutoff frequency [Hz]
        
        
        % specific param for mode='spectrum'
        Freq=[]        % Definition of spectrum response [Hz]
        Spec=[]        % Definition of spectrum response [dB]
    end
    properties (Dependent)
        L0_db
        a
        pnvar
        pnmean
        whitepnvar
    end
    properties (Access = public)
        % 
        currentphase
        LastPhasor=0
    end
    methods
        % constructor
        function out=c_oscillator(SimPars)
            out = out@c_component;
            out.Mode='ideal';
            if nargin==0, SimPars.PN=struct; end
            if isfield(SimPars,'FS'), out.FS=SimPars.FS; end
            if isfield(SimPars,'PN')
                if isfield(SimPars.PN,'Mode'), out.Mode=SimPars.PN.Mode; end
                if isfield(SimPars.PN,'CFO'), out.CFO=SimPars.PN.CFO; end
                if isfield(SimPars.PN,'CFO_STD'), out.CFO_std=SimPars.PN.CFO_STD; end
                if isfield(SimPars.PN,'SpectrumFreqs'), out.Freq=SimPars.PN.SpectrumFreqs; end
                if isfield(SimPars.PN,'SpectrumVals'), out.Spec=SimPars.PN.SpectrumVals; end
                if isfield(SimPars.PN,'L100'), out.L100_db=SimPars.PN.L100; end
                if isfield(SimPars.PN,'F3DB'), out.F3db=SimPars.PN.F3DB; end
                if isfield(SimPars.PN,'LINF'), out.Linf_db=SimPars.PN.LINF; end
            end
           
            out.currentphase=rand*2*pi;
        end
        
        % functions to access dependent data
        function l0=get.L0_db(obj)
            l0=pow2db(1e10*(db2pow(obj.L100_db)-db2pow(obj.Linf_db))/obj.F3db^2);
        end
        function a=get.a(obj)
            a=exp(-2*pi*obj.F3db/obj.FS);
        end    
        function pnvar=get.pnvar(obj)
            a2=obj.a^2;
            if 1-a2<1e-10
                pnvar=4*pi^2*1e10*db2pow(obj.L100_db)/obj.FS;
            else
                pnvar=(1-a2)*pi*1e10*db2pow(obj.L100_db)/obj.F3db;
            end
        end
        function pnmean=get.pnmean(obj)
            pnmean=2*pi*(obj.CFO+obj.CFO_std*randn)/obj.FS; % extra variance to the CFO offset
        end
        function whitepnvar=get.whitepnvar(obj)
            whitepnvar=db2pow(obj.Linf_db)*obj.FS;
        end
        
        % functions to run the oscillator
        function fi=run_phase(obj,nosamples)
            switch obj.Mode
                case 'ideal'
                    fi=zeros(nosamples,1);
                case 'cfo'
                    fi=obj.run_phase_cfo(nosamples);
                case 'model'
                    fi=obj.run_phase_model(nosamples);
                case 'spectrum'
                    fi=obj.run_phase_spectrum(nosamples);
                case 'iid'
                    fi=obj.run_phase_iid(nosamples);
                otherwise
                    error('');
            end
        end        
        function X=run(obj,nosamples)
            X=exp(1i*obj.run_phase(nosamples));
            if strcmp(obj.Mode,'model') % add AWGN
                X=X+sqrt(obj.whitepnvar)*(randn(size(X))+1i*randn(size(X)) );
            end
            obj.LastPhasor=X;
        end
        
        % help functions
        function v=variance_phase_spectrum(obj)
            N=2^23;
            F=[eps,obj.Freq,obj.FS]';
            S=[obj.Spec(1),obj.Spec,obj.Spec(end)]';
            Flin=linspace(0,obj.FS/2,N+1)';Flin=Flin(2:end);
            X=db2pow(interp1(log(F),S,log(Flin))/2);
            X=X*sqrt(Flin(2)-Flin(1))*N*2;
            v=sum(abs(X(2:end)).^2)/2/length(X)^2;
        end
    end
    methods (Access=private) % private, called from other functions
        function fi=run_phase_spectrum(obj,nosamples)
            N=nosamples/2;
            F=([eps;obj.Freq';obj.FS]);
            S=[obj.Spec(1);obj.Spec';obj.Spec(end)];
            Flin=linspace(0,obj.FS/2,N+1)';Flin=Flin(2:end);
            X=db2pow(interp1(log(F),S,log(Flin))/2);
            X=X*sqrt(Flin(2)-Flin(1))*N*2;
            X=X*sqrt(0.5).*(Usefulfunctions.randn_c(size(X)));
            %                    semilogx(Flin,pow2db(abs2(X)));
            fi=real(ifft([0;X(2:end);0;flipud(conj(X(2:end)))]));
            fi=fi+obj.run_phase_cfo(size(fi));
            obj.currentphase=fi(end);
        end
        function fi=run_phase_model(obj,nosamples)
            u=sqrt(obj.pnvar)*randn(nosamples,1)+obj.pnmean;
            fi=filter(1,[1 -obj.a],[obj.currentphase;u]);
            fi=fi(2:end);
            obj.currentphase=fi(end);
        end
        function fi=run_phase_cfo(obj,nosamples)
            fi=cumsum(obj.pnmean*ones(nosamples,1))+obj.currentphase;
            obj.currentphase=fi(end);
        end
    end
end