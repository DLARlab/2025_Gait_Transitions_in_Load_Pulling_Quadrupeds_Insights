classdef SLIP_LegTrajectories_Quad < OutputCLASS
    % Plots quadruped leg angular velocities with limb-specific colors
    % Lines order: [LH, LF, RF, RH]
    
    properties
        fig;       % Figure handle (classic or UIFigure)
        ax;        % Axis handle
        Lines;     % Line handles: [LH, LF, RF, RH]
    end

    methods
        function obj = SLIP_LegTrajectories_Quad(T, Y, Number_of_Strides, PlotPosition, PlotTarget)
            obj.slowDown = 1;
            obj.rate = 0.05;

            % Handle single PlotTarget (either axis or figure)
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
            pbaspect(obj.ax, [2*Number_of_Strides 1 1]);

            T = linspace(0, Number_of_Strides, length(T));

            obj = obj.InitializePlots(T, Y);

        end

        function obj = InitializePlots(obj, T, Y)
            FontSize = 12;

            % ------------------------------------------------------------
            % Leg velocity indices (column indices in Y):
            % Required order for plotting lines: [LH, LF, RF, RH]
            % User-specified indices: 8 10 14 12
            % ------------------------------------------------------------
            leg_vel_indices = [8, 10, 14, 12]; % [LH, LF, RF, RH]
            leg_vel_data = Y(:, leg_vel_indices);

            % ------------------------------------------------------------
            % Limb-specific colors provided by user:
            % [RH RF LF LH]
            % We need to map them to line order [LH, LF, RF, RH]
            % so the color rows applied are [LH<-4, LF<-3, RF<-2, RH<-1]
            % ------------------------------------------------------------
            colors_RH_RF_LF_LH = [
                8   118 179;  % RH
                239 229 72;   % RF
                214 99  8;    % LF
                8   163 119   % LH
            ] / 255;

            % Map to [LH, LF, RF, RH]
            color_order = [4, 3, 2, 1];
            colors = colors_RH_RF_LF_LH(color_order, :);

            cla(obj.ax); hold(obj.ax, 'on');
            obj.Lines = gobjects(1,4);
            for i = 1:4
                obj.Lines(i) = plot(obj.ax, T, leg_vel_data(:, i), '-', ...
                    'LineWidth', 1.5, 'Color', colors(i, :));
            end

            xlabel(obj.ax, 'Stride Time [%]', 'FontSize', FontSize);
            ylabel(obj.ax, 'Angular Velocity $[\sqrt{g/l_0}]$', 'Interpreter', 'LaTeX', 'FontSize', FontSize);
%             legend(obj.ax, { ...
%                 '$\dot{\alpha}_{LH}$', ...
%                 '$\dot{\alpha}_{LF}$', ...
%                 '$\dot{\alpha}_{RF}$', ...
%                 '$\dot{\alpha}_{RH}$' }, ...
%                 'Interpreter', 'LaTeX', 'Location', 'best');
            title(obj.ax, 'Leg Angular Velocities', 'FontSize', FontSize);
            xlim(obj.ax, [0 T(end)]);

            y_min = min(leg_vel_data, [], 'all');
            y_max = max(leg_vel_data, [], 'all');
            if ~isfinite(y_min) || ~isfinite(y_max) || (y_min == y_max)
                y_min = -5; y_max = 5; % sensible fallback
            end
            padding = 0.05 * (y_max - y_min);
            ylim(obj.ax, [min(y_min - padding, -5), max(y_max + padding, 5)]);
        end

        function obj = update(obj, T_, Y_)
            % Update leg velocity lines using the same indices/order
            leg_vel_indices = [8, 10, 14, 12]; % [LH, LF, RF, RH]
            for i = 1:4
                obj.Lines(i).XData = T_;
                obj.Lines(i).YData = Y_(:, leg_vel_indices(i));
            end
        end
    end
end
