function fig_handle = PlotFootfallSequence_Hildebrand_v2(varargin)
% PlotFootfallSequence_Hildebrand_v3  Plot footfall mean and ±1σ errorbars for structured inputs

    narginchk(1,3);
    obj1 = varargin{1};
    [~, seq1_mean, ~, seq1_std, N] = parseSequence(obj1, 1);
    input1IsStruct = isstruct(obj1);

    seq2_mean = [];
    seq2_std = [];
    input2IsStruct = false;
    fig_handle = [];

    if nargin >= 2
        arg2 = varargin{2};
        if isgraphics(arg2,'figure')
            fig_handle = arg2;
        else
            [~, seq2_mean, ~, seq2_std, N2] = parseSequence(arg2, 2);
            input2IsStruct = isstruct(arg2);
            if N2 ~= N
                error('Second sequence must have same stride count.');
            end
        end
    end

    if nargin == 3 && isempty(fig_handle)
        arg3 = varargin{3};
        if isgraphics(arg3,'figure')
            fig_handle = arg3;
        else
            error('Third input must be a figure handle.');
        end
    end

    % Prepare figure
    if isempty(fig_handle)
        fig_handle = figure('Name','Footfall Sequence','NumberTitle','off');
    else
        figure(fig_handle); clf;
    end
    ax = gca; hold(ax,'on');

    % Horizontal offsets
    if isempty(seq2_mean)
        hw = 0.35; centers1 = 1:4;
    else
        hw = 0.25; offset = 0.6*hw;
        centers1 = (1:4) + offset;
        centers2 = (1:4) - offset;
    end

    % Nested plot function
    function plotSeq(meanSeq, stdSeq, centers, isStruct)
        for i = 1:N
            for limb = 1:4
                yc = centers(limb);
                m0 = meanSeq(i,2*limb-1);
                m1 = meanSeq(i,2*limb);
                s0 = stdSeq(i,2*limb-1);
                s1 = stdSeq(i,2*limb);

                if isStruct
                    % Mean block (white with black edge)
                    rectX = [m0 m1 m1 m0];
                    rectY = yc + hw/2 * [-1 -1 1 1];
                    fill(rectX, rectY, 'w', 'EdgeColor', 'k', 'LineWidth', 1.5, 'Parent', ax);
                    % ±1σ errorbars in black
                    % start event
                    plot([m0-s0, m0+s0], [yc, yc], 'k-', 'LineWidth', 1.5, 'Parent', ax);
                    plot([m0-s0, m0-s0], [yc-hw/8, yc+hw/8], 'k-', 'LineWidth', 1, 'Parent', ax);
                    plot([m0+s0, m0+s0], [yc-hw/8, yc+hw/8], 'k-', 'LineWidth', 1, 'Parent', ax);
                    % end event
                    plot([m1-s1, m1+s1], [yc, yc], 'k-', 'LineWidth', 1.5, 'Parent', ax);
                    plot([m1-s1, m1-s1], [yc-hw/8, yc+hw/8], 'k-', 'LineWidth', 1, 'Parent', ax);
                    plot([m1+s1, m1+s1], [yc-hw/8, yc+hw/8], 'k-', 'LineWidth', 1, 'Parent', ax);
                else
                    % Simulation mean block in black
                    rectX = [m0 m1 m1 m0];
                    rectY = yc + hw/2 * [-1 -1 1 1];
                    colors = [
                                 8 163 119;   % LH
                               214  99   8;   % LF
                                 8 118 179;   % RH
                               239 229  72    % RF
                              ] / 255;
                    fill(rectX, rectY, colors(limb,:), 'EdgeColor', 'none', 'Parent', ax);

                end
            end
        end
    end

    % Plot sequences
    plotSeq(seq1_mean, seq1_std, centers1, input1IsStruct);
    if ~isempty(seq2_mean)
        plotSeq(seq2_mean, seq2_std, centers2, input2IsStruct);
    end

    % Axes formatting
    xlim([0 N]); xticks(0:0.25:N); xlabel('Stride Cycle [%]');
    ylim([0.5 4.5]); yticks(1:4); yticklabels({'RH','RF','LF','LH'});
    pbaspect([2*N 1 1]); box(ax,'on'); hold(ax,'off');

    % Legend
    legendEntries = {};
    if input1IsStruct
        legendEntries = {'Mean Value of Experimental Footfall Sequence', 'Deviation Value of Experimental Footfall Sequence'};
    end
    if ~isempty(seq2_mean)
        legendEntries{end+1} = 'Simulation Experimental Footfall';
    end
    if ~isempty(legendEntries)
        legend(ax, legendEntries, 'Location', 'best');
    else
        legend(ax, 'off');
    end
end

function [mi, me, ma, st, N] = parseSequence(obj, idx)
    if isnumeric(obj)
        validateattributes(obj, {'numeric'}, {'2d','ncols',8,'nonnegative','finite'});
        mi = obj; me = obj; ma = obj;
        st = zeros(size(obj));
        N = size(obj,1);
    elseif isstruct(obj)
        fns = fieldnames(obj);
        mi = obj.(fns{contains(lower(fns),'min')});
        me = obj.(fns{contains(lower(fns),'mean')});
        ma = obj.(fns{contains(lower(fns),'max')});
        if any(contains(lower(fns),'std')|contains(lower(fns),'dev'))
            st = obj.(fns{contains(lower(fns),'std')|contains(lower(fns),'dev')});
        else
            st = (ma - mi)/2; % fallback
        end
        N = size(me,1);
    else
        error('Input %d must be numeric N×8 or struct with min/mean/max (and optional std) fields.', idx);
    end
end