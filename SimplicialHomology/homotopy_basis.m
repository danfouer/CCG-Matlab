% 本 MATLAB 代码实现了计算高亏格曲面的贪心同伦基的算法。该算法基于论文 [1] 中的算法。
% 
% 该函数的输入参数为一个网格结构体 `mesh`，输出参数为一个基的单元格数组 `hb`。对于亏格为零的曲面，返回空数组。
% 
% 该函数的实现过程如下：
% 
% 1. 通过 `compute_edge` 函数计算网格的边列表 `edge` 和每个面的边列表 `eif`。
% 
% 2. 通过 `compute_adjacency_matrix` 函数计算网格的邻接矩阵 `am` 和对角线矩阵 `amd`。
% 
% 3. 构造邻接矩阵 `G`，其中边的权重为其长度。
% 
% 4. 通过 `graphshortestpath` 函数或 `dijkstra` 函数计算从基点 `bi` 开始的最短路径 `dist` 和路径上的前驱节点 `pred`。
% 
% 5. 通过 `compute_dual_graph` 函数计算网格的对偶图 `amf`。
% 
% 6. 通过 `pred` 和 `amf` 计算出最大生成树 `tree`。
% 
% 7. 构造邻接矩阵 `G2`，其中不包含在 `tree` 中的边和被 `tree` 中的边穿过的边。
% 
% 8. 构造贪心同伦基，其中每个基都是 `G2` 中的一条边。
% 
% 其中，`compute_edge` 函数用于计算网格的边列表和每个面的边列表，`compute_adjacency_matrix` 函数用于计算网格的邻接矩阵和对角线矩阵，`compute_dual_graph` 函数用于计算网格的对偶图，`graphshortestpath` 函数或 `dijkstra` 函数用于计算最短路径，`minimum_spanning_tree` 函数或 `graphminspantree` 函数用于计算最大生成树。
% 
% 需要注意的是，该函数的实现过程中，使用了其他函数的实现，这些函数的实现过程可以参考其他代码的解释。
%% compute greedy homotopy basis 
% Compute a greedy homotopy group basis based on the algorithm in
% paper[1]. Works for closed surface.
% 
% Two graph algorithms are needed: minimum spanning tree and shortest
% path, provided in Matlab's bioinformatics toolbox. We supply 
% alternatives for these two functions, implemented purely in Matlab,
% which are a little slower and will be invoked when Maltab's built-in
% functions are not available.
%  
% # Erickson, Jeff, and Kim Whittlesey. "Greedy optimal homotopy and 
%   homology generators." Proceedings of the sixteenth annual ACM-SIAM 
%   symposium on Discrete algorithms. Society for Industrial and Applied 
%   Mathematics, 2005.
%
%% Syntax
%   hb = compute_greedy_homotopy_basis(face,vertex,bi)
%
%% Description
%  face  : double array, nf x 3, connectivity of mesh
%  vertex: double array, nv x 3, vertex of mesh
%  bi    : integer scaler, base point of homotopy group
%
%  hb: cell array, n x 1, a basis of homotopy group, each cell is a closed 
%      loop based at bi. Return empty for genus zero surface.
% 
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2013/02/23
%  Revised: 2014/03/22 by Wen, optimize code with vectorized operation,
%           remove dependency on Matlab's built-in function.
%  Revised: 2014/03/24 by Wen, add doc
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function hb = compute_greedy_homotopy_basis(face,vertex,bi)
% greedy_homotopy_basis(face,vertex,i) compute a homotopy
% basis of a high genus surface S = (face,vertex), at a basis point bi;
nf = size(face,1);
% nv = size(vertex,1);
[edge,eif] = compute_edge(face);
[am,amd] = compute_adjacency_matrix(face);
nv = size(am,1);

% G is adjacency matrix constructed from the triangle mesh, with weight the
% edge length.
[I,J] = find(am);
el = sqrt(dot(vertex(I,:)-vertex(J,:),vertex(I,:)-vertex(J,:),2));
G = sparse(I,J,el,nv,nv);
G = (G + G'); % G is undirected

% the shortest path in graph G, with source node bi, T = G_pred is a the tree
if exist('graphshortestpath')
    [dist,path,pred] = graphshortestpath(G,bi,'METHOD','Dijkstra');
    % we always use column array
    dist = dist(:);
    pred = pred(:);
else
    [dist,path,pred] = dijkstra(G,bi);
end

% amf is the dual graph of G
amf = compute_dual_graph(face);

% (G\T)* 
% dual graph G* that do not consist edge correspond to edge in T = pred
I = (1:nv)';
I(bi) = [];
J = pred(I);
% (I2,J2) and (J2,I2) are faces corresponding to edge (I,pred(I))
I2 = full(amd(I+(J-1)*nv));
J2 = full(amd(J+(I-1)*nv));
ind = (I2==0 | J2==0);
I2(ind) = [];
J2(ind) = [];
amf(I2+(J2-1)*nf) = 0;
amf(J2+(I2-1)*nf) = 0;

% tree is the maximum spanning tree of (G\T)*, where the weight of any
% edge e* is length(shortest_loop(e))
[I,J] = find(amf);
ind = (eif(:,1)==-1 | eif(:,2)==-1);
eif(ind,:) = [];
edge(ind,:) = [];
F2E = sparse([eif(:,1);eif(:,2)],[eif(:,2);eif(:,1)],[edge(:,1);edge(:,2)],nf,nf);
ei = [F2E(I+(J-1)*nf),F2E(J+(I-1)*nf)];
dvi = vertex(ei(:,1),:)-vertex(ei(:,2),:);
V = -(dist(ei(:,1))+dist(ei(:,2))+sqrt(dot(dvi,dvi,2)));
amf_w = sparse(I,J,V,nf,nf);
if exist('graphminspantree')
    tree = graphminspantree(amf_w,'METHOD','Kruskal');
else
    tree = minimum_spanning_tree(amf_w);
end

% G2 is the graph, with edges neither in T nor are crossed by edges in T*
G2 = G;
I = (1:nv)';
I(bi) = [];
J = pred(I);
% (I,J) are edges in T
G2(I+(J-1)*nv) = 0;
G2(J+(I-1)*nv) = 0;

[I,J] = find(tree);
% ei are edges in T*
ei = [F2E(I+(J-1)*nf),F2E(J+(I-1)*nf)];
index = [ei(:,1)+(ei(:,2)-1)*nv;ei(:,2)+(ei(:,1)-1)*nv];
G2(index) = 0;

% the greedy homotopy basis consists of all loops (e), where e is an edge 
%  of G2.
G2 = tril(G2);
[I,J] = find(G2);
hb = cell(size(I));
for i = 1:length(I)
    pi = path{I(i)};
    pj = path{J(i)};
    hb{i} = [pi(:);flipud(pj(:))];
end
if isempty(I)
    hb = [];
end
