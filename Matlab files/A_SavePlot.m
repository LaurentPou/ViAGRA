function A_SavePlot(bool_save,h,plot_title)
% Increase the size of the windows (full height, 75% width) and save it as
% pdf and png files

if bool_save
    set(h,'Units','Normalized','OuterPosition',[0 0 0.75 1]); % pos x (px) pos y (px) width x (%) width y (%)
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
    file_name = strrep(plot_title,sprintf('\n'),''); % Delete the '\n'
    file_name = strrep(file_name,'\',''); % Delete the '\'
    print(h,file_name,'-dpdf');

    print(h,file_name,'-dpng');
    
    saveas(h,file_name,'fig');
end

end