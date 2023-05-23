% 该函数用于计算三角网格在单位圆盘上的调和映射，即将三维网格映射到二维圆盘上，保持角度不变。
% 该函数假设三角网格是单连通的。
% 
% 函数的输入参数包括三角网格的面片连接信息 face 和顶点坐标信息 vertex，以及圆盘边界信息 bd 和目标边界信息
% bd_target（可选）。输出参数为每个顶点在圆盘上的二维坐标 uv。
% 
% 函数的具体实现如下：
% 
% 1. 首先获取顶点数量 nv。
% 
% 2. 如果没有输入圆盘边界信息 bd，则调用 compute_bd 函数计算边界。
% 
% 3. 如果没有输入目标边界信息 bd_target，则根据边界顶点的位置计算圆盘上的坐标 uvbd。
% 
% 4. 初始化 uv 为全零矩阵。
% 
% 5. 将边界顶点的坐标设置为预设的圆盘坐标 uvbd。
% 
% 6. 将内部顶点标记为 in。
% 
% 7. 计算拉普拉斯-贝尔特拉米算子 A。
% 
% 8. 将内部顶点和边界顶点分别提取出来，然后将 A 矩阵和右侧向量 rhs 分别缩小为内部顶点的大小。
% 
% 9. 使用线性求解器求解线性方程组 Ain * uvin = rhs，其中 Ain 是内部顶点对应的拉普拉斯-贝尔特拉米算子，
% rhs 是右侧向量，uvin 是内部顶点的圆盘坐标。
% 
% 10. 将内部顶点的圆盘坐标 uvin 赋值给 uv。
% 
% 该函数的作者为 Wen Cheng Feng，代码和文档的最后更新时间为 2014 年 3 月。
%% disk harmonic map 
% Disk harmonic map of a 3D simply-connected surface.
%
%% Syntax
%   uv = disk_harmonic_map(face,vertex)
%
%% Description
%  face  : double array, nf x 3, connectivity of mesh
%  vertex: double array, nv x 3, vertex of mesh
% 
%  uv: double array, nv x 2, uv coordinates of vertex on 2D circle domain
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/18
%  Revised: 2014/03/24 by Wen, add doc
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function uv = disk_harmonic_map(face,vertex,bd,bd_target)
nv = size(vertex,1);
if ~exist('bd','var') || isempty(bd)
    bd = compute_bd(face);
end
if ~exist('bd_target','var') || isempty(bd_target)
    % bl is boundary edge length
    db = vertex(bd,:) - vertex(bd([2:end,1]),:);
    bl = sqrt(dot(db,db,2));
    t = cumsum(bl)/sum(bl)*2*pi;
    t = t([end,1:end-1]);
    % use edge length to parameterize boundary
    uvbd = [cos(t),sin(t)];
else
    uvbd = bd_target;
end
uv = zeros(nv,2);
uv(bd,:) = uvbd;
in = false(nv,1);
in(face) = true;
in(bd) = false;
A = laplace_beltrami(face,vertex);
Ain = A(in,in);
rhs = -A(in,bd)*uvbd;
uvin = Ain\rhs;
uv(in,:) = uvin;
