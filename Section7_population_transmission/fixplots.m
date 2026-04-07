function  fixplots(size)

figHandles = findobj('Type', 'figure');

numfigs=numel(figHandles);

for jjjj=1:numfigs
    figure(figHandles(jjjj).Number);
    set(findall(gcf,'-property','FontSize'),'FontSize',size,'FontName','Times New Roman')
end

end

