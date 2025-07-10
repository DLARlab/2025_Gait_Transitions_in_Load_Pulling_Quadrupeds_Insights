classdef SLIP_PeriodicOrbit_Quad < OutputCLASS 
    properties 
        fig; % The output window
        axes;
        % Patch objects used in the graphical representation
        Orbit;
        Poincare_Section;
        Current_Position;
        Text;
    end
  
    methods
        % Constructor:
        function obj = SLIP_PeriodicOrbit_Quad(Y,plotPositions,FigOrAx,gaitcolor)
            obj.slowDown = 1;      % Run this in real time.
            obj.rate     = 0.05;   % with 25 fps
            
            if isa(FigOrAx, 'matlab.ui.Figure')
                obj.fig = FigOrAx;     clf(obj.fig);
                % Set window properties
                set(obj.fig, 'Name','Periodic Orbit');  % Window title
                set(obj.fig, 'Color','w');          % Background color
                set(obj.fig, 'Renderer','OpenGL');
                set(obj.fig, 'position', plotPositions);
                
                obj.axes = axes(obj.fig);  hold on;

            elseif isa(FigOrAx, 'matlab.graphics.axis.Axes') || isa(FigOrAx, 'matlab.ui.control.UIAxes')
                obj.axes = FigOrAx;
                obj.fig = ancestor(obj.axes, 'figure');
                set(obj.axes, 'Position', plotPositions, 'Box', 'Off')
            else
                error('Third input must be either a figure handle or axes handle.');
            end

            obj = obj.InitializePlots(Y, gaitcolor);
        end

        function obj = InitializePlots(obj, Y, gaitcolor)
            cla(obj.axes, 'reset');
            hold(obj.axes, 'on');

            obj.axes.XLim = [min(Y(:,2))*(1-sign(min(Y(:,2)))*0.1)   max(Y(:,2))*(1+sign(max(Y(:,2)))*0.1)];
            obj.axes.YLim = [min(Y(:,4))*(1-sign(min(Y(:,4)))*0.1)   max(Y(:,4))*(1+sign(max(Y(:,4)))*0.1)];
            obj.axes.ZLim = [min(Y(:,6))-0.02  max(Y(:,6))+0.02];

            xlabel(obj.axes, '$\dot{q}_x  [\sqrt{gl_0}]$', 'Interpreter', 'LaTex', 'FontSize', 15);
            ylabel(obj.axes, '$\dot{q}_z  [\sqrt{gl_0}]$', 'Interpreter', 'LaTex', 'FontSize', 15);
            zlabel(obj.axes, '$\dot{q}_{pitch}  [rad/s]$', 'Interpreter', 'LaTex', 'FontSize', 15);

            obj.axes.Title.String = 'Periodic Orbit of The Solution';

            obj.Orbit = plot3(Y(:,2),Y(:,4),Y(:,6),'LineWidth',2,'Color',[0 0 0],'Parent',obj.axes);

            obj.Poincare_Section = fill3(...
                [min(Y(:,2))*(1-sign(min(Y(:,2)))*0.1), max(Y(:,2))*(1+sign(max(Y(:,2)))*0.1), ...
                 max(Y(:,2))*(1+sign(max(Y(:,2)))*0.1), min(Y(:,2))*(1-sign(min(Y(:,2)))*0.1)], ...
                [0 0 0 0], ...
                [min(Y(:,6))-0.02, min(Y(:,6))-0.02, max(Y(:,6))+0.02, max(Y(:,6))+0.02], ...
                [0.7 0.7 0.7], 'Parent', obj.axes, 'FaceAlpha', 0.3);

            textPS = text(max(Y(:,2))*1.06, 0, (max(Y(:,6)) + min(Y(:,6)))*0.1,...
                           {'Poincare Section','$\dot{q}_z = 0$'}, 'Interpreter','LaTex', ...
                           'HorizontalAlignment','center', 'Parent',obj.axes);

            textPO = text(min(Y(:,2))*0.98, min(Y(:,4)), (max(Y(:,6))+min(Y(:,6)))/2, ...
                          'Periodic Orbit', 'HorizontalAlignment','center', 'Parent',obj.axes);

            obj.Text = struct('textPS',textPS,'textPO',textPO);

            obj.Current_Position = scatter3(Y(1,2),Y(1,4),Y(1,6),120,'filled','Parent',obj.axes,...
                                            'MarkerEdgeColor',gaitcolor,'MarkerFaceColor',gaitcolor);

            obj.axes.View = [45 45];
        end

        function obj = update(obj,y)
            if isa(obj.fig, 'matlab.ui.Figure')
                figure(obj.fig);  % only valid for classic figures
            end
            % update the scatter point regardless
            set(obj.Current_Position,'xData',y(2),'yData',y(4),'zData',y(6));
        end
    end
end
