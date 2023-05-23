% 本 MATLAB 代码实现了计算拓扑空间 $M$ 的一维同调群 $H_1(M,\mathbb{Z})$ 的基的算法。
% 该算法基于书籍 [1] 中的第六个算法。
% 
% 该函数的输入参数为一个网格结构体 `mesh`，输出参数为一个基的单元格数组 `hb`。对于亏格为零的曲面，返回空数组。
% 对于每个手柄，返回两个环。如果曲面有边界，则每个边界都是 `hb` 的一个元素。
% 
% 该函数的实现过程如下：
% 
% 1. 通过 `cut_graph` 函数计算网格的割图，得到边的列表 `ee`。
% 
% 2. 将 `ee` 转换为无向图 `G`。
% 
% 3. 通过 `minspantree` 函数计算 `G` 的最小生成树 `tree`，并记录每个节点的父节点 `pred`。
% 
% 4. 从最小生成树的根节点 `v` 开始，对于每条不在最小生成树上的边 `eh(i,:)`，分别从 `eh(i,1)` 
% 和 `eh(i,2)` 开始，沿着 `pred` 数组追踪路径，得到两个环的路径 `p1` 和 `p2`，
% 然后将 `p1` 和 `p2` 拼接成一个环 `loop`。
% 
% 5. 对于每个环 `loop`，通过 `prune_path` 函数将其修剪为最短的环，并将其存储在 `hb` 的相应位置。
% 
% 其中，`trace_path` 函数用于追踪路径，`prune_path` 函数用于修剪环。
% 
% 该函数的作者为 Wen Cheng Feng，最初创建于 2014 年 3 月 13 日，最近更新于 2018 年 9 月 16 日。
% 该函数的版权归香港中文大学数学系计算几何小组所有。
%% homology basis 
% Compute a basis for the homology group H_1(M,Z), based on the algorithm
% 6 in book [1].
%  
% # Gu, Xianfeng David, and Shing-Tung Yau, eds. Computational conformal
%   geometry. Vol. 3. Somerville: International Press, 2008.
%
%% Syntax
%   hb = homology_basis(mesh)
%
%% Description
%  mesh: mesh structure
%
%  hb: cell array, n x 1, a basis of homology group, each cell is a closed 
%      loop based. Return empty for genus zero surface. Two loops on each
%      handle. If there is boundary on surface, each boundary will be an
%      element of hb
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/13
%  Revised: 2014/03/24 by Wen, add doc
%  Revised: 2018/09/16 by Wen, simplify code
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui
function hb = homology_basis(mesh)
face = mesh.face;
ee = cut_graph(face);
G = graph(ee(:,1),ee(:,2));
[tree,pred] = minspantree(G,'Type','forest','Root',ee(1));
v = ee(1);
tree_edges = table2array(tree.Edges);
I = [tree_edges(:,1);tree_edges(:,2)];
J = [tree_edges(:,2);tree_edges(:,1)];
eh = setdiff(ee,[I,J],'rows');
hb = cell(size(eh,1),1);
for i = 1:size(eh,1)
    p1 = trace_path(pred,eh(i,1),v);
    p2 = trace_path(pred,eh(i,2),v);
    loop = [flipud(p1);eh(i,1);eh(i,2);p2];
    hb{i} = prune_path(loop);
end

function path = trace_path(pred,v,root)
path = [];
while true
    path = [path;pred(v)];
    v = pred(v);
    if v == root
        break;
    end
end

function path_new = prune_path(path)
ind = path ~= flipud(path);
i = find(ind,1)-1;
path_new = path(i:end-i+1);
