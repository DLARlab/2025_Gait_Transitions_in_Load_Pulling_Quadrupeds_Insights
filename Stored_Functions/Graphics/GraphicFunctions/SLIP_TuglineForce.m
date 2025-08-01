classdef SLIP_TuglineForce < OutputCLASS
    properties
        fig;                   % Figure handle
        ax;                    % Axis handle
        LineExp;               % Line handle for experimental force (mean)
        RangeExp;              % Fill handle for experimental force (range)
        LineSim;               % Line handle for simulation force
        F_exp;                 % Experimental tugline force (array or struct)
        F_sim;                 % Simulation tugline force (N x 2 array)
        Number_of_Strides;     % Number of strides (integer number)
    end

    methods
        function obj = SLIP_TuglineForce(varargin)
            % Input format:
            % SLIP_TuglineForce(F_exp, F_sim, PlotPosition, PlotTarget)
            % SLIP_TuglineForce(F_sim, PlotPosition, PlotTarget)

            obj.slowDown = 1;
            obj.rate = 0.05;

            if nargin == 5
                obj.F_sim = varargin{1};
                obj.F_exp = varargin{2};
                obj.Number_of_Strides = varargin{3};
                PlotPosition = varargin{4};
                PlotTarget = varargin{5};
            elseif nargin == 4
                obj.F_sim = varargin{1};
                obj.F_exp = [];
                obj.Number_of_Strides = varargin{2};
                PlotPosition = varargin{3};
                PlotTarget = varargin{4};
            else
                error('SLIP_TuglineForce expects 4 or 5 inputs.');
            end

            if isa(PlotTarget, 'matlab.graphics.axis.Axes')
                obj.ax = PlotTarget;
                obj.fig = ancestor(obj.ax, 'figure');
            elseif isa(PlotTarget, 'matlab.ui.Figure') || isa(PlotTarget, 'matlab.graphics.figure.Figure')
                obj.fig = PlotTarget;
                figure(obj.fig); clf(obj.fig);
                obj.ax = axes('Parent', obj.fig);
            else
                error('PlotTarget must be an axis or figure handle.');
            end

            set(obj.ax, 'Position', PlotPosition, 'Box', 'on');

            obj = obj.InitializePlots();
        end

        function obj = InitializePlots(obj)

            cla(obj.ax); hold(obj.ax, 'on'); 
            obj.LineExp = [];  obj.RangeExp = []; obj.LineSim = [];
        
            % Time vector based on simulation length
            if isprop(obj, 'Number_of_Strides') && ~isempty(obj.Number_of_Strides)
                N = obj.Number_of_Strides;
            else
                N = 1;  % Default to 1 stride if not specified
            end

            T = linspace(0, N, length(obj.F_sim));
            F_sim_val = obj.F_sim;
        
            % Plot simulation tugline force
            obj.LineSim = plot(obj.ax, T, F_sim_val, '-', 'LineWidth', 2, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Sim');
            xlim(obj.ax, [0 T(end)]);
        
            % Plot experimental force if provided
            if ~isempty(obj.F_exp)
                if isstruct(obj.F_exp)
                    fns = fieldnames(obj.F_exp);
                    LF_mean = []; LF_std = [];
                    for i = 1:numel(fns)
                        fn = fns{i}; vec = obj.F_exp.(fn);
                        if numel(vec) ~= numel(T)
                            tNorm = linspace(0, T(end), numel(vec));
                            vec = interp1(tNorm, vec, T, 'linear', 'extrap');
                        end
                        lname = lower(fn);
                        if contains(lname, 'mean')
                            LF_mean = vec;
                        elseif contains(lname, 'std') || contains(lname, 'dev')
                            LF_std = vec;
                        end
                    end
                    
                    if ~isempty(LF_mean)
                        obj.LineExp = plot(obj.ax, T, LF_mean, ':', 'LineWidth', 2.5, 'Color', [0.6350 0.0780 0.1840], 'DisplayName','Exp Mean');
                    end

         
                    if ~isempty(LF_mean) && ~isempty(LF_std)
                        T = T(:);
                        LF_mean = LF_mean(:);
                        LF_std = LF_std(:);
                        obj.RangeExp = fill(obj.ax, [T; flipud(T)], [LF_mean+LF_std; flipud(LF_mean-LF_std)], [0.6350 0.0780 0.1840], 'FaceAlpha',0.3, 'EdgeColor','none', 'DisplayName','Exp \pm1\sigma');
                    end

                    legend(obj.ax, 'Sim','Exp Mean', 'Exp Range', 'Location', 'best');

                else
                    F_exp_val = DataLengthResampler(obj.F_exp, T);
                    obj.LineExp = plot(obj.ax, T, F_exp_val, ':','LineWidth',2.5, 'DisplayName','Exp');
                    legend(obj.ax, 'Sim','Exp','Location', 'best');
                end
            end
        
            xlabel(obj.ax, 'Stride Time  $[\sqrt{l_0/g}]$', 'Interpreter', 'LaTeX');
            ylabel(obj.ax, 'Leash Force  $[mg]$', 'Interpreter', 'LaTeX');
            title(obj.ax, 'Loading Force Along the Rope');
        
            strideTicks = 0:0.25:N;
            xticks(obj.ax, strideTicks);
            xticklabels(obj.ax, compose('%.2f', strideTicks));
        end

        function obj = update(obj, F_exp_new, F_sim_new)
            obj.F_exp = F_exp_new;
            obj.F_sim = F_sim_new;
            obj = obj.InitializePlots();
        end
    end
end
