% 这段 MATLAB 代码实现了计算三角网格的对偶图。对偶图中，原网格的每个面对应一个顶点，顶点的位置为原面的重心。
% 具体实现过程如下：
% 
% 首先，从输入的 mesh 结构中获取 eif 矩阵，它是一个大小为 ne x 2 的矩阵，表示每条边所在的两个面的编号。
% 如果一条边只属于一个面，则另一个面的编号为 0。将 eif 矩阵中为 0 的元素替换为 inf，方便后续处理。
% 
% 然后，根据 eif 矩阵构造一个大小为 ne x 2 的矩阵，其中第一列为 eif 矩阵的第一列，第二列为 eif 矩阵的第二列。
% 对于 eif 矩阵中为 inf 的元素，将其替换为 NaN，方便后续处理。
% 
% 接着，根据 eif 矩阵构造一个大小为 ne x 1 的逻辑向量 ind，其中 ind(i) 表示第 i 条边所在的两个面的编号都不为 0。
% 根据 ind 向量筛选出 eif 矩阵中不为 inf 的行，得到一个大小为 n x 2 的矩阵，其中 n 为满足 ind(i) 的 i 的个数。
% 这个矩阵表示对偶图中的边。
% 
% 最后，使用 MATLAB 自带的 graph 函数构造对偶图。将上一步得到的矩阵作为输入的边，构造一个无向图 G。函数的输出
% 即为 G。
%% compute dual graph 
% Dual graph of a triangle mesh, regarded as graph. 
% Each face in original mesh corresponds to a vertex in dual graph, vertex
% position be the centroid of the original face.
%
%% Syntax
%   [amf] = compute_dual_graph(face);
%   [amf,dual_vertex] = compute_dual_graph(face,vertex);
%
%% Description
%  face  : double array, nf x 3, connectivity of mesh
%  vertex: double array, nv x 3, vertex of mesh
% 
%  amf: sparse matrix, nf x nf, connectivity of dual graph
%  dual_vertex: nf x 3, dual vertex in dual graph, if vertex is not
%               supplied, will return []
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/14
%  Revised: 2014/03/18 by Wen, add doc
%  Revised: 2014/03/23 by Wen, revise doc
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function G = dual_graph(mesh)
eif = mesh.eif;
ind1 = eif(:,1)>0;
eif(~ind1,1) = inf;
ind2 = eif(:,2)>0;
eif(~ind2,2) = inf;
ind = ind1&ind2;
G = graph(eif(ind,1),eif(ind,2));
