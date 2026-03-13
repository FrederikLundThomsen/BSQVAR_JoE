% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

function plt_fmt = PlotFmt(fig,varargin)

    if isempty(varargin)||isempty(varargin{1})
        same_legend = 0;
    else
        same_legend = varargin{1};
    end

    if isempty(varargin)||length(varargin) < 2||isempty(varargin{2})
        gry_scl = 0;
    else
        gry_scl = varargin{2};
    end
    
    %Default color scheme
    dflt_col = {
        [0,50,153]/255,...
        [255,180,0]/255,...
        [255,75,0]/255,...
        [101,184,0]/255,...
        [0,177,234]/255,...
        [0,120,22]/255,...
        [129,57,198]/255,...
        [92,92,92]/255,...
        [152,161,208]/255,...
        [253,221,167]/255,...
        [246,177,131]/255,...
        [206,225,175]/255,...
        [215,238,248]/255,...
        [141,184,141]/255,...
        [174,151,199]/255,...
        [169,169,169]/255,...
        [217,217,217]/255
        };

    dflt_gry = cellfun(@rgb2gray,dflt_col,'UniformOutput',false);

    if gry_scl == 1
        clr_schm = dflt_gry;
    else
        clr_schm = dflt_col;
    end
    
    if (exist('fig') == 0)||(isempty(fig) == 1)
        fig = gcf();
    end
    
    %Reset units of plot to centimeters
    set(fig,'Units','centimeters')

    font_sz = 9;
        
    fchld = fig.Children;
    
    for c = 1:length(fchld)
        if strcmp(fchld(c).Type,'legend')

            lgnd = fchld(c);

            %Set legend properties
            if isempty(lgnd) == 0
                %Adjust legend
                lgnd.Location = 'layout';
                lgnd.Orientation = 'Horizontal';
                lgnd.Box = 'off';
                lgnd.FontName = 'Arial';
                lgnd.FontSize = font_sz;
                lgnd.TextColor = [83,83,83]/255;
                lgnd.ItemTokenSize = lgnd.ItemTokenSize/2;
                while (lgnd.Position(1)<0)||(lgnd.Position(2)<0)||(lgnd.Position(3) > 1)||(lgnd.Position(4) > 1)
                    lgnd.NumColumns = lgnd.NumColumns - 1;
                end        
                lgnd.Units = 'pixels';
            end
        end
        if strcmp(fchld(c).Type,'tiledlayout')
            fchldt = fig.Children(c).Children;
            for ct = 1:length(fchldt)
                if strcmp(fchldt(ct).Type,'axes')
                  
                    %Get axes properties
                    ax = fchldt(ct);

                    %Set general properties
                    ax.Color = 'none';
                    ax.XGrid = 'on';
                    ax.YGrid = 'on';    
                    ax.GridColor = [217,217,217]/255;
                    ax.GridAlpha = 1;
                    ax.LineWidth = 0.3;    
                    ax.Box = 'off';
                    ax.FontName = 'Arial';
                    ax.FontSize = font_sz;
                    if strcmp(ax.XLimMode,'auto')
                        ax.XLimitMethod = 'tight';
                    end
                    if strcmp(ax.YLimMode,'auto')
                        ax.YLimitMethod = 'tight';
                    end
                    ax.Title.FontSize = font_sz;
                    ax.Title.FontName = 'Arial';
                    ax.Title.Color = [83,83,83]/255;
                
                    %Number of default axes
                    yno = length(ax.YAxis);

                    %Loop through all data series and identify 3d plots
                    is2d = 1;
                    for si = 1:length(ax.Children)
                        s = length(ax.Children) - (si-1);
                        if ~(strcmp(ax.Children(s).Type,'constantline')||strcmp(ax.Children(s).Type,'bar'))
                            is2d = min(is2d,isempty(ax.Children(s).ZData));
                        end
                    end

                    %Loop through all series colours
                    plt_col = cell(0);
                    for aks = 1:yno
                        if is2d
                            if aks == 1
                                yyaxis(ax,'left');
                            else
                                yyaxis(ax,'right');
                            end
                        end
                        for si = 1:length(ax.Children)
                            s = length(ax.Children) - (si-1);
                            %To avoid changing colors for y- and x-lines
                            if ~strcmp(ax.Children(s).Type,'constantline')
                                if strcmp(ax.Children(s).Type,'line')
                                    plt_col{end+1} = ax.Children(s).Color;
                                elseif strcmp(ax.Children(s).Type,'contour')
                                    plt_col{end+1} = ax.Children(s).LineColor;
                                elseif strcmp(ax.Children(s).Type,'stair')
                                    plt_col{end+1} = ax.Children(s).Color;
                                elseif strcmp(ax.Children(s).Type,'patch')
                                    plt_col{end+1} = ax.Children(s).FaceColor;
                                    if ~strcmp(ax.Children(s).EdgeColor,'none')
                                        plt_col{end+1} = ax.Children(s).EdgeColor;
                                    end
                                elseif strcmp(ax.Children(s).Type,'bar')
                                    plt_col{end+1} = ax.Children(s).FaceColor;
                                elseif strcmp(ax.Children(s).Type,'scatter')
                                    plt_col{end+1} = ax.Children(s).CData;
                                elseif strcmp(ax.Children(s).Type,'surface')
                                    plt_col{end+1} = ax.Children(s).FaceColor;
                                    if ~strcmp(ax.Children(s).EdgeColor,'none')
                                        plt_col{end+1} = ax.Children(s).EdgeColor;
                                    end
                                end 
                            end
                        end
                    end

                    %Identify unique colours
                    u_plt_col = unique(cell2mat(plt_col(:)),'rows','stable');

                    %Map colours to defaults
                    for aks = 1:yno
                        if is2d
                            if aks == 1
                                yyaxis(ax,'left');
                            else
                                yyaxis(ax,'right');
                            end
                        end
                        for s = 1:length(ax.Children)
                            %To avoid changing colors for y- and x-lines
                            if ~strcmp(ax.Children(s).Type,'constantline')
                                if strcmp(ax.Children(s).Type,'line')
                                    ax.Children(s).Color = clr_schm{find(sum(u_plt_col == ax.Children(s).Color,2) == 3)};
                                    if ~strcmp(ax.Children(s).Marker,'none')
                                        ax.Children(s).MarkerFaceColor = ax.Children(s).Color;
                                    end
                                elseif strcmp(ax.Children(s).Type,'contour')
                                    ax.Children(s).LineColor = clr_schm{find(sum(u_plt_col == ax.Children(s).LineColor,2) == 3)};
                                elseif strcmp(ax.Children(s).Type,'stair')
                                    ax.Children(s).Color = clr_schm{find(sum(u_plt_col == ax.Children(s).Color,2) == 3)};
                                    if ~strcmp(ax.Children(s).Marker,'none')
                                        ax.Children(s).MarkerFaceColor = ax.Children(s).Color;
                                    end
                                elseif strcmp(ax.Children(s).Type,'patch')
                                    ax.Children(s).FaceColor = clr_schm{find(sum(u_plt_col == ax.Children(s).FaceColor,2) == 3)};
                                    if ~strcmp(ax.Children(s).EdgeColor,'none')
                                        ax.Children(s).EdgeColor = clr_schm{find(sum(u_plt_col == ax.Children(s).EdgeColor,2) == 3)};
                                    end
                                elseif strcmp(ax.Children(s).Type,'bar')
                                    ax.Children(s).EdgeColor = clr_schm{find(sum(u_plt_col == ax.Children(s).FaceColor,2) == 3)};
                                    ax.Children(s).FaceColor = clr_schm{find(sum(u_plt_col == ax.Children(s).FaceColor,2) == 3)};
                                elseif strcmp(ax.Children(s).Type,'scatter')
                                    ax.Children(s).CData = clr_schm{find(sum(u_plt_col == ax.Children(s).CData,2) == 3)};
                                elseif strcmp(ax.Children(s).Type,'surface')
                                    ax.Children(s).FaceColor = clr_schm{find(sum(u_plt_col == ax.Children(s).FaceColor,2) == 3)};
                                    if ~strcmp(ax.Children(s).EdgeColor,'none')
                                        ax.Children(s).EdgeColor = clr_schm{find(sum(u_plt_col == ax.Children(s).EdgeColor,2) == 3)};
                                    end
                                end
                            end
                        end
                    end

                    %Get legend
                    lgnd = fchldt(ct).Legend;

                    %Set legend properties
                    if isempty(lgnd) == 0
                        if same_legend == 1
                            %Adjust legend
                            lgnd.Location = 'layout';
                            lgnd.Orientation = 'Horizontal';
                            lgnd.Box = 'off';
                            lgnd.FontName = 'Arial';
                            lgnd.FontSize = font_sz;
                            lgnd.TextColor = [83,83,83]/255;
                            lgnd.ItemTokenSize = lgnd.ItemTokenSize/2;
                            while (lgnd.Position(1)<0)||(lgnd.Position(2)<0)||(lgnd.Position(3) > 1)||(lgnd.Position(4) > 1)
                                lgnd.NumColumns = lgnd.NumColumns - 1;
                            end        
                        else
                            %Adjust legend
                            lgnd.Location = 'northoutside';
                            lgnd.Orientation = 'Horizontal';
                            lgnd.Box = 'off';
                            lgnd.FontName = 'Arial';
                            lgnd.FontSize = font_sz;
                            lgnd.TextColor = [83,83,83]/255;
                            lgnd.ItemTokenSize = lgnd.ItemTokenSize/2;
                            while (lgnd.Position(1)-ax.OuterPosition(1)<0)||(lgnd.Position(2)-ax.OuterPosition(2)<0)||(lgnd.Position(3) > ax.OuterPosition(3))||(lgnd.Position(4) > ax.OuterPosition(4))
                                lgnd.NumColumns = lgnd.NumColumns - 1;
                            end        
                        end
                        lgnd.Units = 'pixels';
                    end
                    
                    %Set X axis properties
                    xax = ax.XAxis;
                    xax.FontName = 'Arial';
                    xax.FontSize = font_sz;
                    xax.TickDirection = 'out';
                    xax.Color = [83,83,83]/255;
                    if ~iscategorical(xax.TickValues)
                        if ~isempty(xax.Exponent)
                            xax.Exponent = 0;
                        end
                    end
                    xax.TickLabelRotation = 0;

                    %Set Z axis properties
                    zax = ax.ZAxis;
                    zax.FontName = 'Arial';
                    zax.FontSize = font_sz;
                    zax.TickDirection = 'out';
                    zax.Color = [83,83,83]/255;
                    if ~iscategorical(zax.TickValues)
                        if ~isempty(zax.Exponent)
                            zax.Exponent = 0;
                        end
                    end
                    zax.TickLabelRotation = 0;
                    
                    %Get Y axis properties
                    ylmt = ax.YLim;
                    ytck = ax.YTick;
                    
                    if is2d
                        %Call secondary Y axis
                        yyaxis(ax,'right');
                        
                        %If a seperate secondary Y axis is not specified, set it
                        %equal to the first Y axis
                        if yno == 1
                            ax.YAxis(2).Limits = ylmt;
                            ax.YAxis(2).TickValues = ytck;
                        end
                    end

                    %Set primary Y axis properties
                    ylax = ax.YAxis(1);
                    ylax.FontName = 'Arial';
                    ylax.FontSize = font_sz;
                    ylax.TickDirection = 'in';
                    ylax.Color = [83,83,83]/255;
                    if ~iscategorical(ylax.TickValues)
                        if ~isempty(ylax.Exponent)
                            ylax.Exponent = 0;
                        end
                    end

                    %Find the appropriate number of decimals
                    if ~iscategorical(ylax.TickValues)
                        if ~isempty(ylax.Exponent)
                            ydec = cellfun(@(c)strsplit(c,'.'),ylax.TickLabels,'UniformOutput',false);
                            ydeca = 0;
                            for d = 1:length(ydec)
                                if length(ydec{d}) == 2
                                    ydeca = max(ydeca,strlength(ydec{d}(2)));
                                end
                            end
                            ylax.TickLabelFormat = strcat('%,.',num2str(ydeca),'f');
                        end
                    end
                    ytcklbl = ylax.TickLabels;

                    if is2d
                        %Set secondary Y axis properties
                        yrax = ax.YAxis(2);
                        yrax.FontName = 'Arial';
                        yrax.FontSize = font_sz;
                        yrax.TickDirection = 'in';
                        yrax.Color = [83,83,83]/255;
                        if ~iscategorical(yrax.TickValues)
                            if ~isempty(yrax.Exponent)
                                yrax.Exponent = 0;
                                yrax.TickLabelFormat = strcat('%,.',num2str(ydeca),'f');
                            end
                        end
                        
                        if yno == 2
                            %Find the appropriate number of decimals
                            if ~iscategorical(yrax.TickValues)
                                if ~isempty(yrax.Exponent)
                                    ydec = cellfun(@(c)strsplit(c,'.'),yrax.TickLabels,'UniformOutput',false);
                                    ydeca = 0;
                                    for d = 1:length(ydec)
                                        if length(ydec{d}) == 2
                                            ydeca = max(ydeca,strlength(ydec{d}(2)));
                                        end
                                    end
                                    yrax.TickLabelFormat = strcat('%,.',num2str(ydeca),'f');
                                end
                            end
                        end
                        
                        if yno == 1
                            yrax.TickLabels = ytcklbl;
                        end
                    end
                    
                    yrax.TickLabelsMode = 'manual';
                    yrax.TickValuesMode = 'manual';

                    ylax.TickLabelsMode = 'manual';
                    ylax.TickValuesMode = 'manual';

                end
            end
        elseif strcmp(fchld(c).Type,'axes')
            
            %Get axes properties
            ax = fchld(c);

            %Set general properties
            ax.Color = 'none';
            ax.XGrid = 'on';
            ax.YGrid = 'on';    
            ax.GridColor = [217,217,217]/255;
            ax.GridAlpha = 1;
            ax.LineWidth = 0.3;    
            ax.Box = 'off';
            ax.FontName = 'Arial';
            ax.FontSize = font_sz;
            if strcmp(ax.XLimMode,'auto')
                ax.XLimitMethod = 'tight';
            end
            if strcmp(ax.YLimMode,'auto')
                ax.YLimitMethod = 'tight';
            end
            ax.Title.FontSize = font_sz;
            ax.Title.FontName = 'Arial';
            ax.Title.Color = [83,83,83]/255;

            %Number of default axes
            yno = length(ax.YAxis);

            %Loop through all data series and identify 3d plots
            is2d = 1;
            for si = 1:length(ax.Children)
                s = length(ax.Children) - (si-1);
                if ~(strcmp(ax.Children(s).Type,'constantline')||strcmp(ax.Children(s).Type,'bar'))
                    is2d = min(is2d,isempty(ax.Children(s).ZData));
                end
            end
            
            %Loop through all series colours
            plt_col = cell(0);
            for aks = 1:yno
                if is2d
                    if aks == 1
                        yyaxis(ax,'left');
                    else
                        yyaxis(ax,'right');
                    end
                end
                for si = 1:length(ax.Children)
                    s = length(ax.Children) - (si-1);
                    %To avoid changing colors for y- and x-lines
                    if ~strcmp(ax.Children(s).Type,'constantline')
                        if strcmp(ax.Children(s).Type,'line')
                            plt_col{end+1} = ax.Children(s).Color;
                        elseif strcmp(ax.Children(s).Type,'contour')
                            plt_col{end+1} = ax.Children(s).LineColor;
                        elseif strcmp(ax.Children(s).Type,'stair')
                            plt_col{end+1} = ax.Children(s).Color;
                        elseif strcmp(ax.Children(s).Type,'patch')
                            plt_col{end+1} = ax.Children(s).FaceColor;
                            if ~strcmp(ax.Children(s).EdgeColor,'none')
                                plt_col{end+1} = ax.Children(s).EdgeColor;
                            end
                        elseif strcmp(ax.Children(s).Type,'bar')
                            plt_col{end+1} = ax.Children(s).FaceColor;
                        elseif strcmp(ax.Children(s).Type,'scatter')
                            plt_col{end+1} = ax.Children(s).CData;
                        elseif strcmp(ax.Children(s).Type,'surface')
                            plt_col{end+1} = ax.Children(s).FaceColor;
                            if ~strcmp(ax.Children(s).EdgeColor,'none')
                                plt_col{end+1} = ax.Children(s).EdgeColor;
                            end
                        end 
                    end
                end
            end

            %Identify unique colours
            u_plt_col = unique(cell2mat(plt_col(:)),'rows','stable');

            %Map colours to defaults
            for aks = 1:yno
                if is2d
                    if aks == 1
                        yyaxis(ax,'left');
                    else
                        yyaxis(ax,'right');
                    end
                end
                for s = 1:length(ax.Children)
                    %To avoid changing colors for y- and x-lines
                    if ~strcmp(ax.Children(s).Type,'constantline')
                        if strcmp(ax.Children(s).Type,'line')
                            ax.Children(s).Color = clr_schm{find(sum(u_plt_col == ax.Children(s).Color,2) == 3)};
                            if ~strcmp(ax.Children(s).Marker,'none')
                                ax.Children(s).MarkerFaceColor = ax.Children(s).Color;
                            end
                        elseif strcmp(ax.Children(s).Type,'contour')
                            ax.Children(s).LineColor = clr_schm{find(sum(u_plt_col == ax.Children(s).LineColor,2) == 3)};
                        elseif strcmp(ax.Children(s).Type,'stair')
                            ax.Children(s).Color = clr_schm{find(sum(u_plt_col == ax.Children(s).Color,2) == 3)};
                            if ~strcmp(ax.Children(s).Marker,'none')
                                ax.Children(s).MarkerFaceColor = ax.Children(s).Color;
                            end
                        elseif strcmp(ax.Children(s).Type,'patch')
                            ax.Children(s).FaceColor = clr_schm{find(sum(u_plt_col == ax.Children(s).FaceColor,2) == 3)};
                            if ~strcmp(ax.Children(s).EdgeColor,'none')
                                ax.Children(s).EdgeColor = clr_schm{find(sum(u_plt_col == ax.Children(s).EdgeColor,2) == 3)};
                            end
                        elseif strcmp(ax.Children(s).Type,'bar')
                            ax.Children(s).EdgeColor = clr_schm{find(sum(u_plt_col == ax.Children(s).FaceColor,2) == 3)};
                            ax.Children(s).FaceColor = clr_schm{find(sum(u_plt_col == ax.Children(s).FaceColor,2) == 3)};
                        elseif strcmp(ax.Children(s).Type,'scatter')
                            ax.Children(s).CData = clr_schm{find(sum(u_plt_col == ax.Children(s).CData,2) == 3)};
                        elseif strcmp(ax.Children(s).Type,'surface')
                            ax.Children(s).FaceColor = clr_schm{find(sum(u_plt_col == ax.Children(s).FaceColor,2) == 3)};
                            if ~strcmp(ax.Children(s).EdgeColor,'none')
                                ax.Children(s).EdgeColor = clr_schm{find(sum(u_plt_col == ax.Children(s).EdgeColor,2) == 3)};
                            end
                        end
                    end
                end
            end

            %Set X axis properties
            xax = ax.XAxis;
            xax.FontName = 'Arial';
            xax.FontSize = font_sz;
            xax.TickDirection = 'out';
            xax.Color = [83,83,83]/255;
            if ~iscategorical(xax.TickValues)
                if ~isempty(xax.Exponent)
                    xax.Exponent = 0;
                end
            end
            xax.TickLabelRotation = 0;

            %Set Z axis properties
            zax = ax.ZAxis;
            zax.FontName = 'Arial';
            zax.FontSize = font_sz;
            zax.TickDirection = 'out';
            zax.Color = [83,83,83]/255;
            if ~iscategorical(zax.TickValues)
                if ~isempty(zax.Exponent)
                    zax.Exponent = 0;
                end
            end
            zax.TickLabelRotation = 0;

            %Get Y axis properties
            ylmt = ax.YLim;
            ytck = ax.YTick;
            
            if is2d
                %Call secondary Y axis
                yyaxis(ax,'right');
                
                %If a seperate secondary Y axis is not specified, set it
                %equal to the first Y axis
                if yno == 1
                    ax.YAxis(2).Limits = ylmt;
                    ax.YAxis(2).TickValues = ytck;
                end
            end

            %Set primary Y axis properties
            ylax = ax.YAxis(1);
            ylax.FontName = 'Arial';
            ylax.FontSize = font_sz;
            ylax.TickDirection = 'in';
            ylax.Color = [83,83,83]/255;
            if ~iscategorical(ylax.TickValues)
                if isempty(ylax.Exponent) == 0
                    ylax.Exponent = 0;
                end
            end

            %Find the appropriate number of decimals
            if is2d
                if ~isempty(ylax.Exponent)
                    ydec = cellfun(@(c)strsplit(c,'.'),ylax.TickLabels,'UniformOutput',false);
                    ydeca = 0;
                    for d = 1:length(ydec)
                        if length(ydec{d}) == 2
                            ydeca = max(ydeca,strlength(ydec{d}(2)));
                        end
                    end
                    ylax.TickLabelFormat = strcat('%,.',num2str(ydeca),'f');
                end
            end
            ytcklbl = ylax.TickLabels;


            if is2d
                %Set secondary Y axis properties
                yrax = ax.YAxis(2);
                yrax.FontName = 'Arial';
                yrax.FontSize = font_sz;
                yrax.TickDirection = 'in';
                yrax.Color = [83,83,83]/255;
                if ~iscategorical(yrax.TickValues)
                    if ~isempty(yrax.Exponent)
                        yrax.Exponent = 0;
                        yrax.TickLabelFormat = strcat('%,.',num2str(ydeca),'f');
                    end
                end
                
                if yno == 2
                    if ~iscategorical(yrax.TickValues)
                        if ~isempty(yrax.Exponent)
                            ydec = cellfun(@(c)strsplit(c,'.'),yrax.TickLabels,'UniformOutput',false);
                            ydeca = 0;
                            for d = 1:length(ydec)
                                if length(ydec{d}) == 2
                                    ydeca = max(ydeca,strlength(ydec{d}(2)));
                                end
                            end
                            yrax.TickLabelFormat = strcat('%,.',num2str(ydeca),'f');
                        end
                    end
                end
                
                if yno == 1
                    yrax.TickLabels = ytcklbl;
                end
            end
        end
    end

        
end
