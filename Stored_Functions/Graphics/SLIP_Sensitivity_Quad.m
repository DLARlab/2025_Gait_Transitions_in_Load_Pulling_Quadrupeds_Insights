classdef SLIP_Sensitivity_Quad < OutputCLASS
    properties
        figs          % top-level GUI ancestor (if any)
        axes          % [2x1] axes handles: (1) line plot, (2) bar plot
        Lines         % struct of plotting handles
        Colors        % custom colors (4x3)
        Data          % cached input data (struct)
    end

    methods
        function obj = SLIP_Sensitivity_Quad(SensitivityStudyData, PlotTargets)
            % Rate settings (not used for timing here, but consistent with your OOD)
            obj.slowDown = 1;
            obj.rate     = 0.05;

            % ===== Validate / capture axes targets (exactly 2) =====
            if ~(isvector(PlotTargets) && numel(PlotTargets) == 2 && ...
                 all(arrayfun(@(ax) isa(ax,'matlab.graphics.axis.Axes'), PlotTargets)))
                error('PlotTargets must be a 2-element vector of axes handles.');
            end
            obj.axes = PlotTargets(:);

            % Ensure both axes share same top-level ancestor (useful for App Designer)
            top_ancestors = arrayfun(@(ax) ancestor(ax,'matlab.ui.Figure'), obj.axes, 'UniformOutput', false);
            if ~isequal(top_ancestors{:})
                error('Both axes must share the same top-level GUI ancestor.');
            end
            obj.figs = top_ancestors{1};

            % ===== Colors =====
            obj.Colors = [...
                  8 163 119;   % LH
                214  99   8;   % LF
                  8 118 179;   % RH
                239 229  72] / 255; % RF

            % ===== Initialize Plots =====
            obj = obj.InitializePlots(SensitivityStudyData);
        end

        function obj = InitializePlots(obj, SensitivityStudyData)
            % -------- Load required fields from structure --------
            % Expecting: C, deltas, names, percs, solution_cells, sortedDeltas, sortedNames
            % (We recompute deltas/sorted* for robustness and consistency with the plotting.)
            S = SensitivityStudyData;
            C = S.C;
            names = S.names;
            percs = S.percs;

            if isrow(percs); percs = percs(:)'; end
            nP = size(C,1);

            % Compute as in the provided script
            C0 = C(:, percs == 0);
            if isempty(C0)
                % Fallback: center column (if 0% column not stored)
                C0 = C(:, floor(size(C,2)/2)+1);
            end
            deltas = 100 * (max(C,[],2) - min(C,[],2)) ./ max(C0, eps);
            [sortedDeltas, idxSort] = sort(deltas, 'descend');
            sortedNames = names(idxSort);

            % Cache data
            obj.Data = struct('C',C, 'percs',percs, 'names',{names}, ...
                              'C0',C0, 'deltas',deltas, ...
                              'sortedDeltas',sortedDeltas, 'sortedNames',{sortedNames}, ...
                              'idxSort',idxSort);

            % ----------- Axes(1): Objective vs %-Perturbation -----------
            ax1 = obj.axes(1);
            cla(ax1); hold(ax1,'on'); grid(ax1,'on');
            FontSize = 12;

            % Plot each parameter curve with markers
            obj.Lines.curves = gobjects(nP,1);
            for i = 1:nP
                ci = obj.Colors(mod(i-1, size(obj.Colors,1))+1, :);
                y = 100*(C(i,:) - C0(i)) / max(C0(i), eps);
                obj.Lines.curves(i) = plot(ax1, percs, y, '-o', ...
                    'Color', ci, ...
                    'DisplayName', names{i}, ...
                    'LineWidth', 1.5);
            end
            xlabel(ax1, '$\mathrm{Perturbation}\ [\%]$', 'Interpreter','latex');
            ylabel(ax1, 'Relative Objective Change [\%]',       'Interpreter','latex');
            title(ax1,  'Objective vs %-Perturbation', 'FontWeight','bold','FontName','Arial');
            legend(ax1, 'Location','best','Interpreter','latex');
            set(ax1,'FontSize',12);
            axis(ax1,'tight');
            xlim(ax1, [min(percs) max(percs)]);

            % ----------- Axes(2): Sorted Percentage Objective Variation -----------
            ax2 = obj.axes(2);
            cla(ax2); grid(ax2,'on'); hold(ax2,'on');

            b = bar(ax2, sortedDeltas, 'FaceColor','flat','EdgeColor','none','BarWidth',0.6);
            b.CData = obj.Colors(idxSort, :); % color bars according to sorted order

            ax2.XTick = 1:nP;
            ax2.XTickLabel = sortedNames;
            ax2.XTickLabelRotation = 0;
            ax2.TickLabelInterpreter = 'latex';
            set(ax2,'FontSize',12);

            xlabel(ax2, '$\mathrm{Paramters}$', 'Interpreter','latex');
            ylabel(ax2, '$\Delta\,\mathrm{Objective}\ [\%]$', 'Interpreter','latex');
            title(ax2,  'Sorted Percentage Objective Variation', 'FontWeight','bold','FontName','Arial');

            ax2.GridLineStyle = '--';
            ax2.GridAlpha = 0.25;
            ylim(ax2, [0, max(sortedDeltas)*1.15]);

            % Text labels above bars
            obj.Lines.bar = b;
            obj.Lines.barTexts = gobjects(nP,1);
            for i = 1:nP
                obj.Lines.barTexts(i) = text(ax2, i, ...
                    sortedDeltas(i) + max(sortedDeltas)*0.02, ...
                    sprintf('%.2f\\%%', sortedDeltas(i)), ...
                    'HorizontalAlignment','center', ...
                    'Interpreter','latex', ...
                    'FontSize',10);
            end
        end

        function obj = update(obj, SensitivityStudyData)
            % Recompute from new data and update visuals in-place
            S = SensitivityStudyData;
            C = S.C;
            names = S.names;
            percs = S.percs;
            if isrow(percs); percs = percs(:)'; end
            nP = size(C,1);

            % Robust color cycling if nP changes
            % (Re-init plots if nP or line count changed)
            needReinit = (~isfield(obj.Lines,'curves')) || numel(obj.Lines.curves) ~= nP;

            % Compute deltas/sorted as before
            C0 = C(:, percs == 0);
            if isempty(C0)
                C0 = C(:, floor(size(C,2)/2)+1);
            end
            deltas = 100 * (max(C,[],2) - min(C,[],2)) ./ max(C0, eps);
            [sortedDeltas, idxSort] = sort(deltas, 'descend');
            sortedNames = names(idxSort);

            % Cache
            obj.Data = struct('C',C, 'percs',percs, 'names',{names}, ...
                              'C0',C0, 'deltas',deltas, ...
                              'sortedDeltas',sortedDeltas, 'sortedNames',{sortedNames}, ...
                              'idxSort',idxSort);

            if needReinit
                % If the dimensionality changed, rebuild the plots
                obj = obj.InitializePlots(SensitivityStudyData);
                return;
            end

            % -------- Update Axes(1) lines --------
            ax1 = obj.axes(1);
            for i = 1:nP
                ci = obj.Colors(mod(i-1, size(obj.Colors,1))+1, :);
                set(obj.Lines.curves(i), 'XData', percs, ...
                                         'YData', 100*(C(i,:) - C0(i))/max(C0(i), eps), ...
                                         'Color', ci, ...
                                         'DisplayName', names{i});
            end
            legend(ax1, 'Location','best','Interpreter','latex','FontSize',11);
            axis(ax1,'tight'); xlim(ax1,[min(percs) max(percs)]);

            % -------- Update Axes(2) bar & labels --------
            ax2 = obj.axes(2);
            b   = obj.Lines.bar;

            % Update bar heights and colors
            b.YData = sortedDeltas(:)';     % ensure row
            b.CData = obj.Colors(idxSort, :);

            % Update ticks and labels
            ax2.XTick = 1:nP;
            ax2.XTickLabel = sortedNames;
            ylim(ax2, [0, max(sortedDeltas)*1.15]);

            % Update text annotations (recreate if count changed)
            if numel(obj.Lines.barTexts) ~= nP || ~all(isgraphics(obj.Lines.barTexts))
                delete(obj.Lines.barTexts(isgraphics(obj.Lines.barTexts)));
                obj.Lines.barTexts = gobjects(nP,1);
                for i = 1:nP
                    obj.Lines.barTexts(i) = text(ax2, i, ...
                        sortedDeltas(i) + max(sortedDeltas)*0.02, ...
                        sprintf('%.2f\\%%', sortedDeltas(i)), ...
                        'HorizontalAlignment','center', ...
                        'Interpreter','latex', ...
                        'FontSize',10);
                end
            else
                for i = 1:nP
                    obj.Lines.barTexts(i).Position = [i, sortedDeltas(i) + max(sortedDeltas)*0.02, 0];
                    obj.Lines.barTexts(i).String   = sprintf('%.2f\\%%', sortedDeltas(i));
                end
            end
        end
    end
end
