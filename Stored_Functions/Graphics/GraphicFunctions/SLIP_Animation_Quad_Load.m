% *************************************************************************
% classdef SLIP_Model_Graphics(p) < OutputCLASS
%
% Two dimensional graphics of a SLIP model.
%
% The graphics object must be initialized with the vector of system
% parameters p.
%
%
% Properties: - NONE
% Methods:    - NONE
%
%
% Created by C. David Remy on 07/10/2011
% MATLAB 2010a - Windows - 64 bit
%
% Documentation:
%  'A MATLAB Framework For Gait Creation', 2011, C. David Remy (1), Keith
%  Buffinton (2), and Roland Siegwart (1),  International Conference on
%  Intelligent Robots and Systems, September 25-30, San Francisco, USA 
%
% (1) Autonomous Systems Lab, Institute of Robotics and Intelligent Systems, 
%     Swiss Federal Institute of Technology (ETHZ) 
%     Tannenstr. 3 / CLA-E-32.1
%     8092 Zurich, Switzerland  
%     cremy@ethz.ch; rsiegwart@ethz.ch
%
% (2) Department of Mechanical Engineering, 
%     Bucknell University
%     701 Moore Avenue
%     Lewisburg, PA-17837, USA
%     buffintk@bucknell.edu
%
%   See also OUTPUTCLASS.
%
classdef SLIP_Animation_Quad_Load< OutputCLASS 
    % Private attributes:
    properties
        States;
        fig;
        axes;
        Body;  
        Leg_BL;
        Leg_FL;
        Leg_BR;
        Leg_FR;
        COM;
        RopeLoad;
        Ground;
        Title;
        options
    end
    % Public methods:
    methods
        % Constructor:
        function obj = SLIP_Animation_Quad_Load(Y, P,plotPositions,FigOrAx,options)
            obj.slowDown = 1;      % Run this in real time.
            obj.rate     = 0.05;   % with 25 fps
            obj.options  = options;
            
            obj.States = Y;

            if isa(FigOrAx, 'matlab.ui.Figure')

                obj.fig = FigOrAx; clf(obj.fig);
                set(obj.fig, 'Name','SLIP model');
                set(obj.fig, 'Color','w');
                set(obj.fig, 'Renderer','OpenGL');
                set(obj.fig, 'Position', plotPositions);

                obj.axes = axes(obj.fig);
                outerpos = obj.axes.OuterPosition;  
                ti = obj.axes.TightInset; 
                left = outerpos(1) + ti(1);
                bottom = outerpos(2) + ti(2);
                ax_width = outerpos(3) - ti(1) - ti(3);
                ax_height = outerpos(4) - ti(2) - ti(4);
                set(obj.axes, 'Position',[left bottom ax_width ax_height],'Box','Off')
                hold(obj.axes, 'on');

            elseif isa(FigOrAx, 'matlab.graphics.axis.Axes') || isa(FigOrAx, 'matlab.ui.control.UIAxes')

                obj.axes = FigOrAx;
                obj.fig = ancestor(obj.axes, 'figure');
                set(obj.axes, 'Position', plotPositions, 'Box', 'Off')
                hold(obj.axes, 'on');

            else
                error('Third input must be either a figure handle or axes handle.');
            end


            obj = obj.InitializePlots(P);
        end

        function obj = InitializePlots(obj, P)
            
            % Clear and reset current axes
            cla(obj.axes, 'reset');
            hold(obj.axes, 'on');
            obj.axes.XAxis.Visible = 'Off';
            obj.axes.YAxis.Visible = 'Off';
            
            Y = obj.States;
            x = Y(1,1); y = Y(1,3); phi_body = Y(1,5);

            lb = P(15); l_leg_rest = P(13);

            obj.Leg_BL = DrawLegs([x-lb;y], l_leg_rest, l_leg_rest, 0.3, obj.axes);
            obj.Leg_FL = DrawLegs([x+(1-lb);y], l_leg_rest, l_leg_rest, -0.1, obj.axes); 
            obj.Ground = DrawGround(obj.axes); 
            obj.Body  = DrawBody([x,y,phi_body],lb, obj.axes);

            if ~(lb == 0.5)
                radius = 0.15;
                obj.COM = DrawCOM([x,y,phi_body],radius, obj.axes);
            end

            obj.Leg_BR = DrawLegs([x-lb;y], l_leg_rest, l_leg_rest, 0.3, obj.axes);
            obj.Leg_FR = DrawLegs([x+(1-lb);y], l_leg_rest, l_leg_rest, -0.1, obj.axes);

            obj.RopeLoad  = DrawRopeLoad( Y(1,:), obj.axes);


            axis(obj.axes, [-3 + x , 1.5+x, -0.1, 2]);
        end

        % Updated function. Is called by the integrator:
        function obj = update(obj, t, y, P)
            
            if isa(obj.fig, 'matlab.ui.Figure')
                figure(obj.fig);  % only valid for classic figures
            end

            [LegLength, LegAngle, BodyJPos, BackJPos, FrontJPos] = ComputeJoint_LegLA(t, y, P);

            % Left Legs
            SetLegs(BackJPos,  LegLength.BL, P(13), LegAngle.BL,obj.Leg_BL);
            SetLegs(FrontJPos, LegLength.FL, P(13), LegAngle.FL,obj.Leg_FL);
            % Main
            SetBody(BodyJPos, obj.Body,P(15));
            % Center of Mass
            if ~(P(15) == 0.5)
                radius = 0.15;
                SetCOM(BodyJPos,obj.COM,radius);
            end         
            % Right Legs
            SetLegs(BackJPos,  LegLength.BR, P(13), LegAngle.BR,obj.Leg_BR);
            SetLegs(FrontJPos, LegLength.FR, P(13), LegAngle.FR,obj.Leg_FR);        
            
            SetRopeLoad(y, obj.RopeLoad);
            
            axis(obj.axes, [-3 + BodyJPos(1),1.5 + BodyJPos(1),-0.1,2]);
            drawnow();
        end
            
            
 
    end
end


%% Functions that draw the parts of the body

% Create patches of main body in the Constructor. Save handles for update function.
function  BodyH = DrawBody(BodyJPos,lb,ax)

    [ x1, y1, f,v ]  = ComputeBodyGraphics(BodyJPos,lb);

    % White background (face color white, no line)
    b1  = patch('XData',x1,'YData',y1,'LineStyle','none','FaceColor',[1 1 1], 'Parent', ax); 
    % Shade lines (Grey lines)
    b2  = patch('faces', f, 'vertices', v,...
           'linewidth',3,'FaceColor',[1 1 1],'EdgeColor',0.8*[1 1 1], 'Parent', ax); 
    % Black Outline (No face color white, black outline)
    b3  = patch('XData',x1,'YData',y1,'linewidth',4,'FaceColor','none', 'Parent', ax);

    % Save all the handles for update function
    BodyH = struct('B_bg',b1,'B_sha',b2,'B_out',b3);

end


% Create patches of legs in the Constructor. Save handles for update function.
function  LegParts=DrawLegs(vecS, l_leg, l_leg_rest, gamma_leg, ax)

    [LegVertices, LegFaces] = ComputeLegGraphics(vecS, l_leg, l_leg_rest, gamma_leg);
    % Spring Part 1************************************************************
    L1   = patch('faces', LegFaces.L_Sp1, 'vertices', LegVertices.L_Sp1,...
           'linewidth',5,'EdgeColor',[245 131 58]/256, 'Parent', ax); % color blue

    % Upper Leg****************************************************************
    % 1. outline of upper leg
    L2_1  = patch('faces',LegFaces.L_UpBO, 'vertices', LegVertices.L_UpBO,...
          'LineStyle','none','FaceColor',[1 1 1], 'Parent', ax);
    % 2. Draw shaded region in the upper leg
    L3   = patch('faces', LegFaces.L_Ups,  'vertices', LegVertices.L_Ups,...
          'linewidth',3,'FaceColor','none','EdgeColor',0.8*[1 1 1], 'Parent', ax); % grey
    % 3. Outline of upper leg  
    L2_2  = patch('faces',LegFaces.L_UpBO, 'vertices', LegVertices.L_UpBO,...
          'linewidth',3,'FaceColor','none', 'Parent', ax);

    % Lower Leg ***************************************************************
    L4  = patch('faces', LegFaces.L_low, 'vertices', LegVertices.L_low,...
          'linewidth',3,'FaceColor',[1 1 1],'EdgeColor',[0 0 0], 'Parent', ax);

    % Spring Part 2 ***********************************************************
    L5   = patch('faces',LegFaces.L_Sp2, 'vertices', LegVertices.L_Sp1,...
          'linewidth',5,'EdgeColor',[245 131 58]/256, 'Parent', ax);

    % Save all the handles for update function
    LegParts = struct('L_Sp1',L1,'L_UpB',L2_1,'L_Upo',L2_2,'L_Ups',L3,'L_low',L4,'L_Sp2',L5);

end

% Create patches of Centor of mass in the Constructor. Save handles for update function.
function COMH = DrawCOM(BodyJPos,radius,ax)

    RotM = [ cos(BodyJPos(3)), -sin(BodyJPos(3));
             sin(BodyJPos(3)),  cos(BodyJPos(3))];
    % Draw Center of Mass
    alpha = linspace(0, pi*2, 40);
    % vert_x_out = sin(alpha)*0.2;
    % vert_y_out = cos(alpha)*0.2;
    vert_out = [sin(alpha)*radius*0.5
                cos(alpha)*radius*0.5];
    vert_out = RotM*vert_out + BodyJPos(1:2)'*ones(1,size(vert_out,2));


    b1 = patch(vert_out(1,:), vert_out(2,:),'white','linewidth',5, 'Parent', ax); 

    alpha = linspace(0, pi/2, 10);
    vert = [0,sin(alpha)*radius*0.75,0
                 0,cos(alpha)*radius*0.75,0];
    vert_x = [vert(1,:);vert(1,:);-vert(1,:);-vert(1,:)];
    vert_y = [vert(2,:);-vert(2,:);-vert(2,:);vert(2,:)];
    Vert = zeros(8,size(vert,2));
    for i = 1:4
        Vert(2*i-1:2*i,:) = RotM*[vert_x(i,:);vert_y(i,:)];
    end

    vert_x = [Vert(1,:);Vert(3,:);Vert(5,:);Vert(7,:)]' + BodyJPos(1);
    vert_y = [Vert(2,:);Vert(4,:);Vert(6,:);Vert(8,:)]' + BodyJPos(2);
    b2 = patch(vert_x, vert_y, cat(3,[1 0 1 0], [1 0 1 0],[1 0 1 0]),'LineWidth',3, 'Parent', ax);


    COMH = struct('B_COMOuter',b1,'B_COMInner',b2);

end

% Create patches of Ground in the Constructor. Save handles for update function.
function GroundH = DrawGround(ax)

        % Draw the ground. It reaches from -2.5 to +6.5.
        h   = 0.01; % Height of the bar at the top
        n   = 5000;  % Number of diagonal stripes in the shaded area
        s   = 0.05; % Spacing of the stripes
        w   = 0.01; % Width of the stripes
        ext = 0.1;  % Length of the stripes

        % Create vertices by shifting a predefined pattern 'n' times to the right:
        v = [     -50,0;
            repmat([     0,    -h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
            repmat([  -ext,-ext-h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
            repmat([-ext+w,-ext-h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
            repmat([     w,    -h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
            -50+s*n+w,0];
        % Connect to faces:
        f = [1,2,4*n+1,4*n+2;
             repmat([0,n,2*n,3*n],n,1) + repmat((1:n)',1,4)+1];
         
        vert_x_out = [-15 100 100 -15];
        vert_y_out = [0 0 -20 -20];
        ground = patch(vert_x_out, vert_y_out,'white', 'Parent', ax);   
        lines  = patch('faces', f, 'vertices', v, 'Parent', ax);
        
        GroundH = struct('ground',ground,'lines',lines);
end

function RopeLoad = DrawRopeLoad(y,ax)

    [Load_x,Load_y] = ComputeLoadGraphics(y);

    Load  = patch('XData',Load_x,'YData',Load_y,'FaceColor',[0 0 0],'FaceAlpha',0.3,...
                  'EdgeColor',[0 0 0],'LineWidth',2, 'Parent', ax);

    xdata_Load   = [y(1);y(1);y(15);y(15)];
    ydata_Load   = [y(3);y(3);y(17);y(17)];
    Rope  = patch('XData',xdata_Load,'YData',ydata_Load,'FaceColor',[0 0 0],'FaceAlpha',0.3,...
                  'EdgeColor',[0 0 0],'LineWidth',2, 'Parent', ax);

    % Load  = patch('XData',Load_x,'YData',Load_y,'LineStyle','none','FaceColor',[0 0 0]);
    
    RopeLoad = struct('Rope',Rope,'Load',Load);
end


%% Set Functions in updating the animation
% Set the patches of main body when figure is updated.
function SetBody(BodyJPos, Body,lb)

    [ x1, y1, f,v ] = ComputeBodyGraphics(BodyJPos,lb);

    set(Body.B_bg , 'xData', x1, 'yData',    y1);                   
    set(Body.B_sha, 'faces',  f, 'vertices',  v);
    set(Body.B_out, 'xData', x1, 'yData',    y1);   

end

% Set the patches of legs when figure is updated.
function  SetLegs( vecS, l_leg, l_leg_rest, gamma_leg, Leghandle)

    [LegVertices, LegFaces] = ComputeLegGraphics( vecS, l_leg, l_leg_rest, gamma_leg);
    % Spring Part 1 (Zigzag line)**********************************************
    set(Leghandle.L_Sp1,'faces', LegFaces.L_Sp1,  'vertices', LegVertices.L_Sp1);
    % Upper Leg****************************************************************
    % 1. Background of upper leg
    set(Leghandle.L_UpB,'faces', LegFaces.L_UpBO, 'vertices', LegVertices.L_UpBO );
    % 2. Shade the upper leg
    set(Leghandle.L_Ups,'faces', LegFaces.L_Ups,  'vertices', LegVertices.L_Ups);
    % 3. Outline the upper leg
    set(Leghandle.L_Upo,'faces', LegFaces.L_UpBO, 'vertices', LegVertices.L_UpBO );
    % Lower Leg ***************************************************************
    set(Leghandle.L_low,'faces', LegFaces.L_low,  'vertices', LegVertices.L_low );
    % Spring Part 2 (parallel)*************************************************
    set(Leghandle.L_Sp2,'faces', LegFaces.L_Sp2,  'vertices', LegVertices.L_Sp1);

end

% Set the patches of center of mass when figure is updated.
function SetCOM(BodyJPos,COM,radius)
    RotM = [ cos(BodyJPos(3)), -sin(BodyJPos(3));
             sin(BodyJPos(3)),  cos(BodyJPos(3))];
    % Draw Center of Mass
    alpha = linspace(0, pi*2, 40);
    % vert_x_out = sin(alpha)*0.2;
    % vert_y_out = cos(alpha)*0.2;
    vert_out = [sin(alpha)*radius*0.5
                cos(alpha)*radius*0.5];
    vert_out = RotM*vert_out + BodyJPos(1:2)'*ones(1,size(vert_out,2));

    alpha = linspace(0, pi/2, 10);
    vert = [0,sin(alpha)*radius*0.75,0
                 0,cos(alpha)*radius*0.75,0];
    vert_x = [vert(1,:);vert(1,:);-vert(1,:);-vert(1,:)];
    vert_y = [vert(2,:);-vert(2,:);-vert(2,:);vert(2,:)];
    Vert = zeros(8,size(vert,2));
    for i = 1:4
        Vert(2*i-1:2*i,:) = RotM*[vert_x(i,:);vert_y(i,:)];
    end

    vert_x = [Vert(1,:);Vert(3,:);Vert(5,:);Vert(7,:)]' + BodyJPos(1);
    vert_y = [Vert(2,:);Vert(4,:);Vert(6,:);Vert(8,:)]' + BodyJPos(2);

    set(COM.B_COMOuter, 'xData', vert_out(1,:), 'yData', vert_out(2,:))
    set(COM.B_COMInner, 'xData', vert_x, 'yData', vert_y)
end




function SetRopeLoad(y,RopeLoad)

    [load_x,load_y] = ComputeLoadGraphics(y);
    set(RopeLoad.Load, 'xData', load_x,'yData', load_y)


    xdata_load   = [y(1);y(1);y(15);y(15)];
    ydata_load   = [y(3);y(3);y(17);y(17)];
    set(RopeLoad.Rope, 'xData', xdata_load,'yData', ydata_load)
    
end
