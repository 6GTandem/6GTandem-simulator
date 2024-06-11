classdef (Abstract) c_component < handle
    properties (Access = public)
%        Mode {mustBeTextScalar} ='ideal';     
        Mode  ='ideal';     
% other things: FS, Delay, drawing functions, 
    end
    methods
        function obj=c_component(mode, varargin)
            if nargin==1, obj.Mode=mode; end
        end
    end
    methods (Abstract) % these methods must be redefined in derived classes.
        y=run(obj,x);
    end
end
