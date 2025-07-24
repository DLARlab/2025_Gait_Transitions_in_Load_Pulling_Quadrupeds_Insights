classdef SLIP_LegTrajectories_Quad < OutputCLASS
    properties
        fig;       % Figure handle (classic or UIFigure)
        ax;        % Axis handle
        Lines;     % Line handles: [LH, LF, RF, RH]
    end

    methods
        function obj = SLIP_LegTrajectories_Quad(T, Y, PlotTarget)
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

            obj = obj.InitializePlots(T, Y);
        end

        function obj = InitializePlots(obj, T, Y)
            FontSize = 12;

            % Leg velocity indices: [LH, LF, RF, RH] => [8, 10, 14, 12]
            leg_vel_indices = [8, 10, 14, 12];
            leg_vel_data = Y(:, leg_vel_indices);

            cla(obj.ax); hold(obj.ax, 'on'); grid(obj.ax, 'on');
            obj.Lines = plot(obj.ax, T, leg_vel_data, '-');
            xlabel(obj.ax, 'Stride Time $[\sqrt{l_0/g}]$', 'Interpreter', 'LaTeX', 'FontSize', FontSize);
            ylabel(obj.ax, 'Angular Velocity [rad/s]', 'Interpreter', 'LaTeX', 'FontSize', FontSize);
            legend(obj.ax, {'$\dot{\alpha}_{LH}$', '$\dot{\alpha}_{LF}$', '$\dot{\alpha}_{RF}$', '$\dot{\alpha}_{RH}$'}, 'Interpreter', 'LaTeX');
            title(obj.ax, 'Leg Angular Velocities', 'FontSize', FontSize);
            xlim(obj.ax, [0 T(end)]);

            y_min = min(leg_vel_data, [], 'all');
            y_max = max(leg_vel_data, [], 'all');
            padding = 0.05 * (y_max - y_min);
            ylim(obj.ax, [y_min - padding, y_max + padding]);
        end

        function obj = update(obj, T_, Y_)
            % Update leg velocity lines
            leg_vel_indices = [8, 10, 14, 12];
            for i = 1:4
                obj.Lines(i).XData = T_;
                obj.Lines(i).YData = Y_(:, leg_vel_indices(i));
            end
        end
    end
end
