letterTriTrans = transpose(allMatrices.letterTri);
heatmap(letterTriTrans, 'Colormap', parula)
set(gca,'XData',wordlist, 'YData', wordlist)
set(gca, 'FontSize', 8)

