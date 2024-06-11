% c_radiostripe
% This class implements a radiostripe in the 6GTANDEM project.
% The stripe consists of a transmitter component, and a set of links.
% Usage:> RS = c_radiostripe;
%       > Y = RS.run(x); % Y is a matrix 
classdef c_radiostripe < c_component

    properties (Access = public)
        TX=c_transmitter;    % A transmitter;
        Link=c_link;         % a vector of links
        BW = 5e9;
        OS = 5;
    end

    methods
        function out=c_radiostripe(nolinks)
            if nargin<1, nolinks=3; end
            out = out@c_component;

            MaxPower=10; % dBm
            AveragePower=5; % dBm
            out.TX=c_transmitter;
            out.TX.PA.Mode='6gtandem';
            out.TX.PA.setmaxpower(MaxPower);
            out.TX.PA.setaveragepower(15,AveragePower);

            out.Link(nolinks)=c_link;

            for i=1:nolinks
                out.Link(i).Amp.Mode='6gtandem';
                out.Link(i).Amp.setmaxpower(MaxPower);
                out.Link(i).Amp.setaveragepower(AveragePower-out.Link(i).Fiber.Damping-out.Link(i).COUP_IN.Damping-out.Link(i).COUP_OUT.Damping,AveragePower);
                out.Link(i).Amp.setnoisevar(300,out.BW*out.OS,10);
            end
        end

        function y=run(obj,x)
            y=zeros(size(x,1),1+length(obj.Link));
            y(:,1)=obj.TX.run(x);
            for i=1:length(obj.Link)
                y(:,i+1)=obj.Link(i).run(y(:,i));
            end
        end

        function calibrate(obj,x, DesiredAmplifierDBM)
        % function calibrate
        % This functions sets the small-signal gain of the link amplifiers
        % to give an approximately constant power (DesiredAmplifierDBM) at the output of each
        % link. The calibration is valid for a given input signal, and
        % recalibration must be performed if the signal statistics changes.
            import Usefulfunctions.*

            for j=1:3
                z=obj.TX.run(x);
                scale=db2mag(DesiredAmplifierDBM-Usefulfunctions.getdbm(z));
                obj.TX.PA.Gain=obj.TX.PA.Gain*scale;
            end
            z=obj.TX.run(x);
            for i=1:length(obj.Link)
                for j=1:3
                    z2=obj.Link(i).run(z);
                    scale=db2mag(DesiredAmplifierDBM-Usefulfunctions.getdbm(z2));
                    obj.Link(i).Amp.Gain=obj.Link(i).Amp.Gain*scale;
                end
                z=obj.Link(i).run(z);
            end
        end


    end
end
