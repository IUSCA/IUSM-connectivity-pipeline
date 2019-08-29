function save_figure(gcf,path2fig,figName,dpi)

%print(gcf,fullfile(path2fig,'figure2.png'),'-dpng',['-r',num2str(600)],'-opengl') %save file 
saveas(gcf,fullfile(path2fig,sprintf('%s.fig',figName)),'fig');
print(gcf,fullfile(path2fig,sprintf('%s.eps',figName)),'-depsc2') %save file
print(gcf,fullfile(path2fig,sprintf('%s.png',figName)),'-dpng') %save file

%pause(2);
pixels = round(250*dpi/25.4);
%system(sprintf('convert -density %d -resize %dx -colorspace RGB %s %s',dpi,pixels,fullfile(path2fig,sprintf('%s.eps',figName)),fullfile(path2fig,sprintf('%s.png',figName))));
