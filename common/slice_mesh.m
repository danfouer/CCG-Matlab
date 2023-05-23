% 该函数的作用是将网格沿着一组边切开，通常这些边来自于cut_graph(directly)、compute_greedy_homotopy_basis
% 和compute_homology_basis等函数，这些函数需要从基础的闭合环中形成边。
% 
% 函数的输入参数有三个，分别是原始网格的面片(face)、顶点(vertex)和需要切开的边(ee)。
% 其中，face是一个nf x 3的双精度数组，表示网格的连接关系；
% vertex是一个nv x 3的双精度数组，表示网格的顶点坐标；ee是一个n x 2的双精度数组，表示需要切开的边，
% 每一行表示一条边，可能不是连续的。
% 
% 函数的输出参数有三个，分别是切开后的新网格的面片(face_new)、顶点(vertex_new)和
% 每个新顶点对应的原始网格顶点(father)。其中，face_new是一个nf x 3的双精度数组，
% 表示切开后的新网格的连接关系；vertex_new是一个nv' x 3的双精度数组，表示切开后的新网格的顶点坐标，
% 其中nv'大于原始网格的顶点数，因为切开网格会将每个ee上的顶点分成两个或更多个；
% father是一个nv' x 1的双精度数组，表示每个新顶点对应的原始网格顶点。
% 
% 函数的具体实现过程如下：
% 
% 1. 通过原始网格的面片(face)计算邻接矩阵，并将其存储在amd中。
% 
% 2. 如果输入的ee是一个单元素的cell数组，则将其转换为一个二维数组ees，ees中的每一行表示一个边。
% 
% 3. 根据ee构建一个稀疏矩阵G，其中G(i,j)表示原始网格中是否存在一条从i到j的边。
% 
% 4. 对于ee中的每个顶点ev，找到与之相邻的顶点集合vre，并按照与ev相邻的边的顺序重新排列vre中的顶点。
% 
% 5. 对于每个evr中的相邻顶点对(i,j)，如果它们在原始网格中构成了一个面片，
% 则将该面片中与ev相邻的顶点替换为一个新的顶点nv+k，并将该新顶点的坐标设置为ev的坐标。
% 同时，将该新顶点的father设置为ev。
% 
% 6. 将所有新的顶点和原始网格的顶点合并到一起，得到新的顶点数组vert_new。
% 
% 7. 将face_new中的面片索引替换为vert_new中的新顶点索引。
% 
% 8. 将father中的原始网格顶点索引和新顶点索引合并到一起，得到新的father数组。
% 
% 9. 将face_new和vert_new合并成一个新的网格mesh_new，并将father设置为mesh_new的属性。最后返回mesh_new。
%% slice_mesh 
% Slice mesh open along a collection of edges ee, which usually comes from 
% cut_graph(directly), or compute_greedy_homotopy_basis and
% compute_homology_basis (need to form edges from closed loops in basis).
% ee can form a single closed loops or multiple closed loops. 
% 
%% Syntax
%   [face_new,vertex_new,father] = slice_mesh(face,vertex,ee)
%
%% Description
%  face  : double array, nf x 3, connectivity of mesh
%  vertex: double array, nv x 3, vertex of mesh
%  ee    : double array, n x 2, a collection of edges, each row is an edge on 
%          mesh, may not be in consecutive order. 
% 
%  face_new  : double array, nf x 3, connectivity of new mesh after slice
%  vertex_new: double array, nv' x 3, vertex of new mesh, vertex number is
%              more than original mesh, since slice mesh will separate each
%              vertex on ee to two vertices or more.
%  father    : double array, nv' x 1, father indicates the vertex on original
%              mesh that new vertex comes from.
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/17
%  Revised: 2014/03/24 by Wen, add doc
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function mesh_new = slice_mesh(mesh,ee)
face = mesh.face;
vert = mesh.vert;
nv = mesh.nv;
[~,amd] = adjacency_matrix(face);
if iscell(ee)
    ees = [];
    for i = 1:length(ee)
        ei = ee{i};
        ei(:,2) = ei([2:end,1]);
        ees = [ees;ei];
    end
    ee = ees;
end
G = sparse(ee(:,1),ee(:,2),ones(size(ee,1),1),nv,nv);
G = G+G';

ev = unique(ee(:));
vre = vert_vert_ring(mesh,ev,true);
face_new = face;
vert2 = zeros(size(ee,1)*2,3);
father2 = zeros(size(ee,1)*2,1);
k = 1;
for i = 1:size(ev,1)
    evr = vre{i};
    for i0 = 1:length(evr)
        if G(evr(i0),ev(i))
            break;
        end
    end
    if evr(1) == evr(end) % interior point
        evr = evr([i0:end-1,1:i0]);
    else % boundary point
        evr = evr([i0:end,1:i0-1]);
    end
    for j = 2:length(evr)
        fi = amd(evr(j),ev(i));
        if fi
            fij = face_new(fi,:)==ev(i);
            face_new(fi,fij) = nv+k;
        end
        if G(ev(i),evr(j))
            vert2(k,:) = vert(ev(i),:);
            father2(k) = ev(i);
            k = k+1;
        end
    end
end
vert_new = [vert;vert2];
father = (1:nv)';
father = [father;father2];

fu = unique(face_new);
index = zeros(max(fu),1);
index(fu) = (1:size(fu,1));
face_new = index(face_new);
vert_new = vert_new(fu,:);
father = father(fu);
mesh_new = make_mesh(face_new,vert_new);
mesh_new.father = father;
