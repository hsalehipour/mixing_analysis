function saveplot(filename)


fontname = 'Times';
fontsize = 18;
h = gcf;
ha = gca;

set(h,'PaperPositionMode','auto')
set(ha, 'FontName',fontname, 'FontSize',fontsize-2);
% grid on;

saveas(h,['../results/',filename], 'fig')
saveas(h,['../results/',filename,'.eps'], 'psc2')
saveas(h,['../results/',filename], 'png')

%print('-depsc','-painter','-r75',['results/',filename])

end