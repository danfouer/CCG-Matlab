% 这是一个用于计算网格边缘的 MATLAB 函数。该函数可以找到网格的边缘，并且可以指示每条边缘所在的两个面。
% 对于边界边缘，只有一个面与之相连，另一个面用 -1 表示。
% 
% 该函数的输入参数为一个 nf x 3 的双精度数组 face，表示网格的连接关系。输出参数为一个 ne x 2 的双精度数组 edge，
% 表示网格的边缘，以及一个 ne x 2 的双精度数组 eif，表示每条边缘所在的两个面。
% 
% 该函数的实现过程如下：
% 
% 1. 调用 adjacency_matrix 函数计算面的邻接矩阵 am 和度矩阵 amd。
% 
% 2. 调用 find 函数找到邻接矩阵中非零元素的行列下标，然后筛选出行下标小于列下标的元素，得到边缘的起点和终点。
% 
% 3. 调用 xor 函数计算邻接矩阵 amd 和 am 的异或矩阵，然后用 find 函数找到非零元素的行列下标和对应的值，
% 得到每条边缘所在的两个面。
% 
% 4. 将起点、终点和面的下标组成的数组合并成一个 ne x 2 的数组 eif，作为函数的输出之一。
% 
% 5. 将起点和终点组成的数组合并成一个 ne x 2 的数组 edge，作为函数的输出之一。
%% compute edge 
% Find edge of mesh, undirected. eif indicates the faces in which the edge
% lies in. For boundary edge, there is only one face attached, in such
% case, the other one is indicated with -1.
% 
% Use this function to replace edgeAttachments method of
% trianglulation/TriRep class.
%
%% Syntax
%   [edge,eif] = edges(face)
%
%% Description
%  face: double array, nf x 3, connectivity of mesh
% 
%  edge: double array, ne x 2, undirected edge
%  eif : double array, ne x 2, each row indicates two faces in which the
%        edge lies in, -1 indicates a boundary edge
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/
%  Revised: 2014/03/18 by Wen, add doc
%  Revised: 2014/03/23 by Wen, revise doc
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function [edge,eif] = edges(face)
[am,amd] = adjacency_matrix(face);
[I,J,~] = find(am);
ind = I<J;
edge = [I(ind),J(ind)];
[~,~,V] = find(amd-xor(amd,am));
[~,~,V2] = find((amd-xor(amd,am))');
eif = [V(ind),V2(ind)];
