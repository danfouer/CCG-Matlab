% 这是一个用于计算网格的切割图的 MATLAB 函数。通过切割图，可以使得网格表面变得单连通。该函数有两个版本：
% 如果提供了面和顶点，则调用版本 1；如果只提供了面，则调用版本 2。版本 1 是书籍 [1] 中算法 3 的精确实现。
% 版本 2 是从 David Gu 的 C++ 代码中翻译而来的，比版本 1 快得多。
% 
% 虽然版本 1 考虑了顶点，但两个算法都不能生成最优的切割图（最短的切割图）。事实上，
% 这个问题直到现在似乎仍然是开放的。
% 
% 语法：
% 
% ee = cut_graph(face)
% 
% 参数说明：
% 
% - face：双精度数组，nf x 3，网格的连接性。
% 
% 返回值：
% 
% - ee：双精度数组，n x 2，切割图中的边，每行是网格上的一条边，可能不是连续的顺序。ee 的主要目的是传递给 slice_mesh，它将沿着 ee 中的边将网格切开，形成一个单连通表面。
% 
% 函数实现过程：
% 
% - 首先，通过面的连接性计算邻接矩阵 am 和 amd。
% - 然后，使用数组模拟队列，将第一个面加入队列中。
% - 对于队列中的每个面，遍历其三个边，找到相邻的面，并将其加入队列中。
% - 在遍历过程中，将已经遍历过的面标记为 true。
% - 遍历结束后，得到一个邻接矩阵 G，表示切割图。
% - 最后，对切割图进行修剪，得到最终的切割图 ee。
% 
% 该函数的作者是 Wen Cheng Feng，属于香港中文大学数学系计算几何小组。
%% cut graph 
% Compute a cut graph of mesh, such that surface becomes simply-connected
% if slice mesh along the cut-graph. There are two versions: if both face 
% and vertex provided, invoke version 1; if only face provided, invoke 
% version 2. Version 1 is exact implementation of algorithm 3 in book [1].
% Version 2 is translated from David Gu's C++ code of cut graph, which is
% much faster than version 1.
% 
% Though version 1 takes vertex into consideration, both algorithms do not
% generated optimal cut graph (shortest one). In fact this problem seems
% to be open until now.
%
%% Syntax
%   ee = cut_graph(face)
%
%% Description
%  face  : double array, nf x 3, connectivity of mesh
% 
%  ee: double array, n x 2, edges in the cut graph, each row is an edge on 
%      mesh, may not be in consecutive order. ee's mainly purpose is been 
%      passed to slice_mesh, which will slice the mesh open along edges in
%      ee, to form a simply-connected surface
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/13
%  Revised: 2014/03/13 by Wen, implement another cut graph algorithm
%  Revised: 2014/03/17 by Wen, merge two cut graph algorithm
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui
% 这段代码将一个由三角面片组成的网格图割成多个不相连的部分，并返回连接这些部分的边缘列表
function ee = cut_graph(face)%  1. 该函数接受一个面列表作为输入，并计算该面列表所表示的三角形网格图的邻接矩阵。
[am,amd] = adjacency_matrix(face);% 2. 使用 adjacency_matrix(face) 函数计算邻接矩阵，并返回两个输出：'am' - 邻接矩阵和 'amd' - 有向邻接矩阵。 
nf = size(face,1);% 3. 确定 'face' 的大小，并将一个队列初始化为模拟先进先出列表。 
% use array to emulate queue
queue = zeros(nf,1);
queue(1) = 1;
qs = 1; % point to queue start
qe = 2; % point to queue end

ft = false(nf,1); % 4. 将网格的第一个面添加到队列中，并使用一个布尔数组 'ft' 初始化，将第一个元素设置为 true。
ft(1) = true;
face4 = face(:,[1 2 3 1]);

while qe > qs %5. while 循环运行，直到队列为空。对于队列中的每个面，检查其三个相邻面，以查看是否存在共享的边缘。 
    fi = queue(qs);
    qs = qs+1;
    for i = 1:3
        he = face4(fi,[i i+1]);
        sf = amd(he(2),he(1));
        if sf % 6. 如果找到共享的边缘，并且相邻的面之前未被访问过，则将其添加到队列中，标记为已访问，并从邻接矩阵中删除边缘。      
            if ~ft(sf)
                queue(qe) = sf;
                qe = qe+1;
                ft(sf) = true;
                am(he(1),he(2)) = -1;
            end
        end
    end
end
am((am<0)') = 0; % 7. 一旦队列为空，将 'am' 中的负数条目更改为 0，并通过获取所有 'am' 中正数条目的上三角矩阵创建二进制邻接矩阵 'G'。
G = triu(am>0);

% prune the graph cut
while true % 8. while 循环通过删除只连接到图的其余部分的一条边缘的任何顶点来修剪图。 
    Gs = full(sum(G,2))+full(sum(G,1))';
    ind = (Gs == 1);
    if sum(ind) ==0
        break;
    end
    G(ind,:) = 0;
    G(:,ind) = 0;
end

[I,J,~] = find(G);
ee = [I,J]; % 9. 最终，将连接不相连部分的边缘作为矩阵 'ee' 返回，其中包含两列 - 连接顶点的索引。 
% 注意：该代码使用广度优先搜索算法遍历整个网格并查找连接的组件。