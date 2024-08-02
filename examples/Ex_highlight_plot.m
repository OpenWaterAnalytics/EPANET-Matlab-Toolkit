start_toolkit;

d = epanet('Net1.inp');

node_index = [5, 6];
node_index_other = [2, 4];
node_index_color = 'r';
node_index_other_color = 'k';

link_index = [2, 5];
link_index_other = [6, 9];
alllinks = [link_index, link_index_other];
link_index_color = '#FFA500'; % orange
link_index_other_color = 'm';
link_width = 2; 

d.plot
legend('off')
coor = d.getNodeCoordinates;
x = coor{1}(node_index);
y = coor{2}(node_index);

x_other = coor{1}(node_index_other);
y_other = coor{2}(node_index_other);

nodesconnlinks = d.getNodesConnectingLinksID;
nodenameid = d.getNodeNameID;

% plot links
linkHandles = [];
linkHandles_other = [];
for i = 1:length(alllinks)
    FromNode = find(strcmp(nodesconnlinks(alllinks(i), 1), nodenameid), 1);
    ToNode = find(strcmp(nodesconnlinks(alllinks(i), 2), nodenameid), 1);
    if FromNode
        x1 = double(coor{1}(FromNode));
        y1 = double(coor{2}(FromNode));
    end
    if ToNode
        x2 = double(coor{1}(ToNode));
        y2 = double(coor{2}(ToNode));
    end

    if ismember(alllinks(i), link_index_other)
        h1 = line([x1 coor{3}{alllinks(i)} x2], [y1 coor{4}{alllinks(i)} y2], 'LineWidth', link_width, 'Color', link_index_other_color);
        linkHandles_other = [linkHandles_other; h1];
    else
        h2 = line([x1 coor{3}{alllinks(i)} x2], [y1 coor{4}{alllinks(i)} y2], 'LineWidth', link_width, 'Color', link_index_color);
        linkHandles = [linkHandles; h2];
    end
end

% plot nodes
nodeHandles = plot(x, y, 'o', 'LineWidth', 2, 'MarkerEdgeColor', node_index_color, 'MarkerFaceColor', node_index_color, 'MarkerSize', 12);
nodeHandles_other = plot(x_other, y_other, 'o', 'LineWidth', 2, 'MarkerEdgeColor', node_index_other_color, 'MarkerFaceColor', node_index_other_color, 'MarkerSize', 12);

fontweight = 'bold';
fontsize = 11;
for i = 1:length(node_index)
    text(x(i)-5, y(i), d.getNodeNameID{node_index(i)}, 'Color', 'black', 'FontWeight', fontweight, 'Fontsize', fontsize)
end

% Create dummy handles for legend
linkHandles_dummy = line(NaN, NaN, 'Color', link_index_color);
linkHandles_other_dummy = line(NaN, NaN, 'Color', link_index_other_color);

% Combine node and link handles for a single legend
allHandles = [nodeHandles; nodeHandles_other; linkHandles_dummy; linkHandles_other_dummy];
allLabels = {'Nodes Red', 'Nodes Black', 'Links Orange', 'Links Magenta'};

% Adding combined legend
legend(allHandles, allLabels, 'Location', 'northeast'); % southeast, etc..

hold off;
