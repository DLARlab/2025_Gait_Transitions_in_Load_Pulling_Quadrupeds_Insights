classdef SLIP_FootfallSequence_Quad < OutputCLASS
    properties
        fig;           % Figure handle
        ax;            % Axis handle
        Patches;       % Patch handles for stance phase (exp + sim)
        ffseq_exp;     % Experimental footfall sequence (struct or N x 8 array)
        ffseq_sim;     % Simulation footfall sequence (N x 8 array)
    end

    methods
        function obj = SLIP_FootfallSequence_Quad(varargin)
            % Usage:
            % SLIP_FootfallSequence_Quad(FFSeq_Sim, PlotPosition, PlotTarget)
            % SLIP_FootfallSequence_Quad(FFSeq_Exp, FFSeq_Sim, PlotPosition, PlotTarget)

            narginchk(3,4);

            if nargin == 3
                FFSeq_Sim = varargin{1};
                FFSeq_Exp = [];
                PlotPosition = varargin{2};
                PlotTarget = varargin{3};
            else
                FFSeq_Sim = varargin{1};
                FFSeq_Exp = varargin{2};
                PlotPosition = varargin{3};
                PlotTarget = varargin{4};
            end

            % Initialize OutputCLASS properties
            obj.slowDown = 1;
            obj.rate = 0.05;

            % Handle axes or figure input
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

            set(obj.ax, 'Position', PlotPosition, 'Box', 'off');

            obj.ffseq_exp = FFSeq_Exp;
            obj.ffseq_sim = FFSeq_Sim;

            obj = obj.InitializePlots();
        end

        function obj = InitializePlots(obj)
            cla(obj.ax); hold(obj.ax, 'on'); box(obj.ax, 'on');

            % Limb-specific colors: [RH RF LF LH]
            colors = [
                8 118 179;
                239 229 72;
                214  99  8;
                8 163 119
            ] / 255;

            hw = 0.25;
            strideTicks = [];
            patchIdx = 0;
            
            % Plot experimental (offset up)
            if ~isempty(obj.ffseq_exp)
                if isstruct(obj.ffseq_exp)
                    fns = fieldnames(obj.ffseq_exp);
                    seq_mean = obj.ffseq_exp.(fns{contains(lower(fns),'exp')});
                    if any(contains(lower(fns),'std') | contains(lower(fns),'dev'))
                        seq_std = obj.ffseq_exp.(fns{contains(lower(fns),'std') | contains(lower(fns),'dev')});
                    else
                        seq_min = obj.ffseq_exp.(fns{contains(lower(fns),'min')});
                        seq_max = obj.ffseq_exp.(fns{contains(lower(fns),'max')});
                        seq_std = (seq_max - seq_min) / 2;
                    end

                    N = size(seq_mean,1);
                    strideTicks = 0:0.25:N;
                    centers_exp = (4:-1:1) + 0.15;

                    for i = 1:N
                        for limb = 1:4
                            m0 = seq_mean(i,2*limb-1); m1 = seq_mean(i,2*limb);
                            s0 = seq_std(i,2*limb-1); s1 = seq_std(i,2*limb);
                            yc = centers_exp(limb);
                            rectX = [m0 m1 m1 m0];
                            rectY = yc + hw/2 * [-1 -1 1 1];

                            patchIdx = patchIdx + 1;
                            obj.Patches(patchIdx) = fill(obj.ax, rectX, rectY, 'w', 'EdgeColor', 'k', 'LineWidth', 1.5);

                            plot(obj.ax, [m0-s0, m0+s0], [yc, yc], 'k-', 'LineWidth', 1.5);
                            plot(obj.ax, [m0-s0, m0-s0], [yc-hw/8, yc+hw/8], 'k-', 'LineWidth', 1);
                            plot(obj.ax, [m0+s0, m0+s0], [yc-hw/8, yc+hw/8], 'k-', 'LineWidth', 1);

                            plot(obj.ax, [m1-s1, m1+s1], [yc, yc], 'k-', 'LineWidth', 1.5);
                            plot(obj.ax, [m1-s1, m1-s1], [yc-hw/8, yc+hw/8], 'k-', 'LineWidth', 1);
                            plot(obj.ax, [m1+s1, m1+s1], [yc-hw/8, yc+hw/8], 'k-', 'LineWidth', 1);
                        end
                    end
                else
                    N = size(obj.ffseq_exp, 1);
                    strideTicks = 0:0.25:N;
                    centers_exp = (4:-1:1) + 0.15;

                    for i = 1:N
                        for limb = 1:4
                            m0 = obj.ffseq_exp(i, 2*limb-1);
                            m1 = obj.ffseq_exp(i, 2*limb);
                            yc = centers_exp(limb);
                            rectX = [m0 m1 m1 m0];
                            rectY = yc + hw/2 * [-1 -1 1 1];

                            patchIdx = patchIdx + 1;
                            obj.Patches(patchIdx) = fill(obj.ax, rectX, rectY, 'w', 'EdgeColor', 'k', 'LineWidth', 1.5);
                        end
                    end
                end
            end

            % Plot simulation (offset down)
            if ~isempty(obj.ffseq_sim)
                N = size(obj.ffseq_sim, 1);
                if isempty(strideTicks)
                    strideTicks = 0:0.25:N;
                end
                centers_sim = (4:-1:1) - 0.15;

                for i = 1:N
                    for limb = 1:4
                        m0 = obj.ffseq_sim(i, 2*limb-1);
                        m1 = obj.ffseq_sim(i, 2*limb);
                        yc = centers_sim(limb);
                        rectX = [m0 m1 m1 m0];
                        rectY = yc + hw/2 * [-1 -1 1 1];

                        patchIdx = patchIdx + 1;
                        obj.Patches(patchIdx) = fill(obj.ax, rectX, rectY, colors(limb,:), 'EdgeColor', 'none');
                    end
                end
            end

            xlim(obj.ax, [0 N]);
            ylim(obj.ax, [0.5 4.5]);
            xticks(obj.ax, strideTicks);
            yticks(obj.ax, 1:4);  
            yticklabels(obj.ax, {'RH','RF','LF','LH'});  % bottom-to-top order


            xlabel(obj.ax, 'Stride Time  [%]');
            title(obj.ax, 'Footfall Sequence');
            pbaspect(obj.ax, [2*N 1 1]);


            % ----- Create Legend Handles -----
            legendHandles = [];
            legendLabels  = [];
            
            % Simulation limbs (RH, RF, LF, LH)
            ytl = yticklabels(obj.ax);  % Get tick labels once
            for limb = 1:4
                h = patch('Parent', obj.ax, 'XData', nan, 'YData', nan, ...
                            'FaceColor', colors(limb,:), 'EdgeColor', 'none');
                legendHandles(end+1) = h; %#ok<AGROW>
                legendLabels{end+1} = ['Sim - ', ytl{limb}]; %#ok<AGROW>
            end

            legendHandles = flip(legendHandles);
            legendLabels = flip(legendLabels);
            
            % Add experimental legend entry if applicable
            if ~isempty(obj.ffseq_exp)
                hExp = patch('Parent', obj.ax, 'XData', nan, 'YData', nan, ...
                               'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5);
                legendHandles = [hExp, legendHandles];
                legendLabels  = ['Exp', legendLabels];
            end
            

            if ~isempty(obj.ax) && isvalid(obj.ax)
                legend(obj.ax, legendHandles, legendLabels, 'Location', 'northeastoutside', 'Box', 'off');
            else
                warning('Axes handle is invalid at legend creation.');
            end

            % Display legend
            if N > 1
            else
                legend(obj.ax, legendHandles, legendLabels, 'Location', 'best');
            end
    
        end

        function obj = update(obj, ffseq_exp_, ffseq_sim_)
            obj.ffseq_exp = ffseq_exp_;
            obj.ffseq_sim = ffseq_sim_;
            obj = obj.InitializePlots();
        end
    end
end
