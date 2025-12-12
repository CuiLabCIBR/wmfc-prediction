function plot_block_matrix(Net_real, p_mat, row_names, col_names, ...
    triMode, titleStr, saveName, cmap)

Net_plot = Net_real;
[nr, nc] = size(Net_real);

mask_valid = true(nr, nc);

switch lower(triMode)
    case 'upper'
        
        lower_mask = tril(true(nr, nc), -1);
        Net_plot(lower_mask)   = NaN;
        mask_valid(lower_mask) = false;

    case 'lower'
        
        upper_mask = triu(true(nr, nc), 1);
        Net_plot(upper_mask)   = NaN;
        mask_valid(upper_mask) = false;

    case 'full'
        
    otherwise
        error('triMode  ''full'', ''upper'' or ''lower''');
end

% ---------- figure ----------
figure('Position', [100 100 550 500]);
ax = axes;

% 
set(ax, 'Color', [1 1 1]);

% NaN
hImg = imagesc(Net_plot, 'Parent', ax);
colormap(ax, cmap);
set(hImg, 'AlphaData', ~isnan(Net_plot));

vals = Net_real(~isnan(Net_real));
if isempty(vals)
    maxAbs = 1;
else
    maxAbs = max(abs(vals(:)));
    if maxAbs == 0
        maxAbs = 1;
    end
end
caxis(ax, [-maxAbs, maxAbs]);


cb = colorbar('peer', ax);
cb.FontSize = 12;
cb.FontName = 'Arial';


set(ax, 'XTick', 1:nc, 'XTickLabel', col_names, ...
    'YTick', 1:nr, 'YTickLabel', row_names, ...
    'TickLength', [0 0], ...
    'FontName', 'Arial', 'FontSize', 12);
ax.XTickLabelRotation = 45;


set(ax, 'YDir', 'reverse');
axis(ax, 'square');
box(ax, 'off');   

hold(ax, 'on');


baseColor = 'k';
baseLW    = 0.5;


for i = 1:nr
    for j = 1:nc
        if ~mask_valid(i,j)
            continue;  
        end

        x_left   = j - 0.5;
        x_right  = j + 0.5;
        y_top    = i - 0.5;
        y_bottom = i + 0.5;

        switch lower(triMode)

            case 'full'
               
                line([x_left x_right x_right x_left x_left], ...
                    [y_top  y_top   y_bottom y_bottom y_top], ...
                    'Color', baseColor, 'LineWidth', baseLW, 'Parent', ax);

              
                p_val = p_mat(i,j);
                if ~isnan(p_val) && p_val < 0.05
                    text(j, i + 0.1, '*', ...
                        'Color', 'w', ...
                        'FontSize', 28, ...
                        'FontWeight', 'bold', ...
                        'FontName', 'Arial', ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'Parent', ax);
                end
            case 'lower'
             

               
                line([x_left x_left], [y_top y_bottom], ...
                    'Color', baseColor, 'LineWidth', baseLW, 'Parent', ax);
                line([x_left x_right], [y_bottom y_bottom], ...
                    'Color', baseColor, 'LineWidth', baseLW, 'Parent', ax); 

               
                if i > 1
                    line([x_left x_right], [y_top y_top], ...
                        'Color', baseColor, 'LineWidth', baseLW, 'Parent', ax);
                end

               
                if j < nc
                    line([x_right x_right], [y_top y_bottom], ...
                        'Color', baseColor, 'LineWidth', baseLW, 'Parent', ax);
                end

               
                p_val = p_mat(i,j);
                if ~isnan(p_val) && p_val < 0.05
                    text(j, i + 0.05, '*', ...
                        'Color', 'w', ...
                        'FontSize', 28, ...
                        'FontWeight', 'bold', ...
                        'FontName', 'Arial', ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'Parent', ax);
                end

            case 'upper'
                
        end
    end
end

hold(ax, 'off');


title(ax, titleStr, ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');


if nargin >= 7 && ~isempty(saveName)
   
    set(gcf, 'Renderer', 'painters');
    set(gcf, 'PaperPositionMode', 'auto');

  
    saveas(gcf, [saveName '.fig']);

    % 300 dpi PNG（位图）
    print(gcf, [saveName '.png'], '-dpng',  '-r300');

    % 300 dpi TIF（位图）
    print(gcf, [saveName '.tif'], '-dtiff', '-r300');

    % PDF（矢量，适合 Adobe）
    print(gcf, [saveName '.pdf'], '-dpdf', '-bestfit');

    % EPS（矢量，适合 Adobe Illustrator）
    print(gcf, [saveName '.eps'], '-depsc', '-painters');
end
end