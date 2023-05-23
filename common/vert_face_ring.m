% 该函数用于计算三角网格中每个顶点的一环邻域面片，可以选择是否按照逆时针顺序排列。默认情况下，不进行顺序排列，
% 基于函数 vert_vert_ring。
% 
% 函数的输入参数包括三角网格的面片连接信息 face，顶点集合 vc（可选），以及是否需要顺序排列 ordered（可选）。
% 输出参数为每个顶点的一环邻域面片 vfr，是一个 cell 数组。
% 
% 函数的具体实现如下：
% 
% 1. 首先获取顶点数量 nv。
% 
% 2. 如果没有输入顶点集合 vc，则默认为所有顶点。
% 
% 3. 如果没有输入是否需要顺序排列 ordered，则默认为 false。
% 
% 4. 调用函数 vert_vert_ring 计算每个顶点的一环邻域顶点集合 vr。
% 
% 5. 将每个顶点的一环邻域面片存储在稀疏矩阵 eifs 中，
% 其中第 i 行 j 列的元素表示以顶点 i 和 j 为起点的半边对应的面片编号。
% 
% 6. 对于每个顶点 i，通过 eifs(vr{i}+nv*(vc(i)-1)) 获取其一环邻域面片编号，
% 然后使用 full 函数将其转换为完整的数组。最后使用 cellfun 函数将每个数组转换为 cell 数组，并去除其中的零元素。
% 
% 该函数的作者为 Wen Cheng Feng，代码和文档的最后更新时间为 2014 年 3 月。
%% compute_vertex_face_ring 
% Compute one-ring neighbor faces of given vertex or all vertex, with or 
% without ccw order. Default is no order, based on vertex_ring
%
%% Syntax
%   vfr = vert_face_ring(face)
%   vfr = vert_face_ring(face,vc)
%   vfr = vert_face_ring(face,vc,ordered)
%
%% Description
%  face: double array, nf x 3, connectivity of mesh
%  vc  : double array, n x 1 or 1 x n, vertex collection, can be empty, 
%        which equivalent to all vertex.
%  ordered: bool, scaler, indicate if ccw order needed.
% 
%  vfr: cell array, nv x 1, each cell is one ring neighbor face, which is
%      a double array
%
%% Example
%   % compute one ring of all vertex, without order
%   vfr = vert_face_ring(face)
% 
%   % compute one ring of vertex 1:100, without ccw order
%   vfr = vert_face_ring(face,1:100,false)
%
%   % compute one ring of vertex 1:100, with ccw order
%   vfr = vert_face_ring(face,1:100,true)
% 
%   % compute one ring of all vertex, with ccw order (may be slow)
%   vfr = vert_face_ring(face,[],true)
% 
%   % same with last one
%   vfr = vert_face_ring(face,1:nv,true)
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/28
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function vfr = vert_face_ring(mesh,vc,ordered)
nv = max(max(mesh.face));
if ~exist('vc','var') || isempty(vc)
	vc = (1:nv)';
end
if ~exist('ordered','var')
    ordered = false;
end
vr = vert_vert_ring(mesh,vc,ordered);
he = mesh.halfedge;
heif = mesh.heif;
eifs = sparse(he(:,1),he(:,2),heif);
vfr = arrayfun(@(i) full(eifs(vr{i}+nv*(vc(i)-1))),(1:length(vc))','UniformOutput',false);
vfr = cellfun(@(vi) vi(vi>0),vfr,'UniformOutput',false);
