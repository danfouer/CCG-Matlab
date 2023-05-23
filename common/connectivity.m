% 这段 MATLAB 代码实现了半边数据结构，将三角网格的面连接信息转换为易于访问的形式。具体来说，
% 它实现了以下四个稀疏矩阵：
% 
% - vvif：大小为 nv x nv，其中 (i,j) 元素表示边 (i,j) 所在的面的编号。
% - nvif：大小为 nf x nv，其中 (i,j) 元素表示面 i 中顶点 j 的下一个顶点编号。
% - pvif：大小为 nf x nv，其中 (i,j) 元素表示面 i 中顶点 j 的上一个顶点编号。
% 
% 这些矩阵可以方便地用于访问三角网格的拓扑结构。具体实现过程如下：
% 
% 首先，将面的三个顶点编号分别存储在 fi、fj 和 fk 中，将面的编号从 1 到 nf 存储在 ff 中。
% 
% 然后，根据边的两个顶点编号，构造 vvif 矩阵。将 fi、fj 和 fk 按列连接起来，得到大小为 3nf x 1 的列向量，
% 分别表示每个面的三条边的两个顶点编号。将 ff 按列连接三次，得到大小为 3nf x 1 的列向量，分别表示每个面的编号。
% 将这两个列向量作为稀疏矩阵的行、列索引，将每个面的编号作为稀疏矩阵的值，构造出 vvif 矩阵。
% 
% 接着，根据面的编号和顶点的编号，构造 nvif 和 pvif 矩阵。将 ff 按列连接三次，得到大小为 3nf x 1 的列向量，
% 分别表示每个面的编号。将 fi、fj 和 fk 按列连接起来，得到大小为 3nf x 1 的列向量，分别表示每个面的三个顶点的编号。
% 将 fi、fj 和 fk 按列连接起来，得到大小为 3nf x 1 的列向量，分别表示每个面的三个顶点的下一个顶点的编号。
% 将 fi、fj 和 fk 按列连接起来，得到大小为 3nf x 1 的列向量，分别表示每个面的三个顶点的上一个顶点的编号。
% 将这三个列向量作为稀疏矩阵的行、列索引，构造出 nvif 和 pvif 矩阵。
% 
% 最后，将 vvif、nvif 和 pvif 作为函数的输出。
%% compute connectivity
% Transform connectivity face to other form to assistant easy access from
% face to vertex, or vise verse. 
% 
% * From face we can access vertex in each face. 
% * From vvif we can access face given two vertex of an edge. 
% * From nvif we can access the "next" vertex in a face and one vertex of
% the face, "next" in the sense of ccw order.
% * From pvif we can access the "previous" vertex in a face and one vertex of
% the face, "previous" in the sense of ccw order.
% 
% Basically, we implement halfedge structure in sparse matrix form. 
%
%% Syntax
%   [vvif,nvif,pvif] = connectivity(face)
%
%% Description
%  face: double array, nf x 3, connectivity of mesh
% 
%  vvif: sparse matrix, nv x nv, element (i,j) indicates the face in which
%        edge (i,j) lies in
%  nvif: sparse matrix, nf x nv, element (i,j) indicates next vertex of
%        vertex j in face i
%  pvif: sparse matrix, nf x nv, element (i,j) indicates previous vertex of
%        vertex j in face i
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/06
%  Revised: 2014/03/23 by Wen, add doc
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function [vvif,nvif,pvif] = connectivity(face)
fi = face(:,1);
fj = face(:,2);
fk = face(:,3);
ff = (1:size(face,1))';
vvif = sparse([fi;fj;fk],[fj;fk;fi],[ff;ff;ff]);
nvif = sparse([ff;ff;ff],[fi;fj;fk],[fj;fk;fi]);
pvif = sparse([ff;ff;ff],[fj;fk;fi],[fi;fj;fk]);
