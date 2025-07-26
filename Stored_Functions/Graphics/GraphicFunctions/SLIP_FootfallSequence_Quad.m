classdef SLIP_FootfallSequence_Quad < OutputCLASS
    properties
        fig;       % Figure handle
        ax;        % Axis handle
        Patches;   % Patch handles for stance phase
    end

    methods
        function obj = SLIP_FootfallSequence_Quad(FFSeq, PlotPosition, PlotTarget)
            % FFSeq: N x 8 matrix [RH_start RH_end RF_start RF_end LF_start LF_end LH_start LH_end]

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

            set(obj.ax, 'Position', PlotPosition, 'Box', 'Off')

            obj = obj.InitializePlots(FFSeq);
        end

        function obj = InitializePlots(obj, FFSeq)
            cla(obj.ax); hold(obj.ax, 'on'); box(obj.ax, 'on');

            % Limb-specific colors: [RH RF LF LH]
            colors = [
                8 118 179;
                239 229 72;
                214  99  8;
                8 163 119
            ] / 255;

            N = size(FFSeq,1);  % number of strides
            strideTicks = 0:0.25:N;
            centers = 1:4;
            hw = 0.35;

            for i = 1:N
                for limb = 1:4
                    m0 = FFSeq(i, 2*limb-1);  % start
                    m1 = FFSeq(i, 2*limb);    % end
                    yc = centers(limb);
                    rectX = [m0 m1 m1 m0];
                    rectY = yc + hw/2 * [-1 -1 1 1];

                    obj.Patches(i, limb) = fill(obj.ax, rectX, rectY, colors(limb,:), ...
                        'EdgeColor', 'none');
                end
            end

            xlim(obj.ax, [0 N]);
            ylim(obj.ax, [0.5 4.5]);
            xticks(obj.ax, strideTicks);
            yticks(obj.ax, 1:4);
            yticklabels(obj.ax, {'RH','RF','LF','LH'});
            xlabel(obj.ax, 'Stride Cycle [%]');
            title(obj.ax, 'Footfall Sequence');
            pbaspect(obj.ax, [2*N 1 1]);
        end

        function obj = update(obj, FFSeq_)
            [N, ~] = size(FFSeq_);

            for i = 1:N
                for limb = 1:4
                    m0 = FFSeq_(i, 2*limb-1);
                    m1 = FFSeq_(i, 2*limb);
                    yc = limb;
                    hw = 0.35;
                    rectX = [m0 m1 m1 m0];
                    rectY = yc + hw/2 * [-1 -1 1 1];
                    set(obj.Patches(i, limb), 'XData', rectX, 'YData', rectY);
                end
            end
        end
    end
end
