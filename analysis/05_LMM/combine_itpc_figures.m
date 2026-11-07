function combine_itpc_figures(sourceFigures,labels,outputStem)
n = numel(sourceFigures);
assert(n>0 && numel(labels)==n);
set(sourceFigures(1),'Units','centimeters');
sz = get(sourceFigures(1),'Position');
w = sz(3)*n; h = sz(4);
combined = figure('Visible','off','Color','w','Renderer','painters', ...
    'Units','centimeters','Position',[1 1 w h]);
for j=1:n
    axs = findall(sourceFigures(j),'Type','axes');
    legends = findall(sourceFigures(j),'Type','legend');
    originals = [axs(:);legends(:)];
    for k=1:numel(originals); set(originals(k),'Units','normalized'); end
    copied = copyobj(originals,combined);
    for k=1:numel(copied)
        pos = get(originals(k),'Position');
        set(copied(k),'Units','normalized','Position',[(pos(1)+j-1)/n,pos(2),pos(3)/n,pos(4)]);
    end
    annotation(combined,'textbox',[(j-1+.03)/n,1-.42/h,.94/n,.30/h], ...
        'String',labels(j),'EdgeColor','none','Margin',0,'HorizontalAlignment','center', ...
        'FontName','Arial','FontSize',7,'FontWeight','bold');
end
set(combined,'PaperUnits','centimeters','PaperSize',[w h], ...
    'PaperPosition',[0 0 w h],'PaperPositionMode','manual');
drawnow;
print(combined,char(outputStem+".pdf"),'-dpdf','-painters');
close(combined); close(sourceFigures);
end
