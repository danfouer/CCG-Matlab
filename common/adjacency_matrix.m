% 这是一个用于计算三角网格邻接矩阵的 MATLAB 函数。该函数只需要连接性信息即可计算出无向和有向邻接矩阵，
% 并以稀疏矩阵的形式存储。对于无向邻接矩阵 am，am(i,j) 的值为 2，表示顶点 i 和顶点 j 之间有一个邻接关系。
% 对于有向邻接矩阵 amd，amd(i,j) 存储一个面的索引，表示半边 (i,j) 所在的面。
% 
% 因此，am 存储了边的信息，而 amd 存储了半边的信息。
% 
% 语法：
% 
% [am,amd] = adjacency_matrix(face);
% 
% 参数说明：
% 
% - face：双精度数组，nf x 3，网格的连接性。
% 
% 返回值：
% 
% - am：稀疏矩阵，nv x nv，无向邻接矩阵。
% - amd：稀疏矩阵，nv x nv，有向邻接矩阵。
% 
% 函数实现过程：
% 
% - 首先，将面的连接性信息转换为三个向量 I、J 和 V，其中 I 和 J 分别表示半边的起点和终点，V 表示半边所在的面的索引。
% - 然后，使用稀疏矩阵函数 sparse 创建有向邻接矩阵 amd。
% - 最后，通过 spones 函数将有向邻接矩阵 amd 转换为无向邻接矩阵 am，并将其与其转置相加，得到最终的无向邻接矩阵 am。
% 
% 该函数的作者是 Wen Cheng Feng，属于香港中文大学数学系计算几何小组。
%% adjacency_matrix
% Compute adjacency matrix of triangle mesh, only connectivity needed.
% 
% Both undirected and directed adjacency matrix computed, stored with 
% sparse matrix. For undirected one (am), am(i,j) has value 2, indicating
% an adjacency between vertex i and vertex j. For directed one (amd),
% amd(i,j) stores a face index, indicating which face the halfedge (i,j) 
% lies in.
% 
% So am stores the information of edge, while amd stores the information
% of half edge.
% 
%% Syntax
%   [am,amd] = adjacency_matrix(face);
% 
%% Description
%  face: double array, nf x 3, connectivity of mesh
% 
%  am : sparse matrix, nv x nv, undirected adjacency matrix
%  amd: sparse matrix, nv x nv, directed adjacency matrix
% 
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/03
%  Revised: 2014/03/14 by Wen, add doc
%  Revised: 2014/03/23 by Wen, remove example
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function [am,amd] = adjacency_matrix(face)
nf = size(face,1);
I = reshape(face',nf*3,1);
J = reshape(face(:,[2 3 1])',nf*3,1);
V = reshape(repmat(1:nf,[3,1]),nf*3,1);
amd = sparse(I,J,V);
am = spones(amd);
am = am+am';
