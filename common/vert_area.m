% 该函数用于计算三角网格中每个顶点周围的面积，有两种计算方式可选：一环和混合。
% 其中，“一环”是指只考虑与该顶点相邻的三角形的面积，而“混合”则是考虑该顶点周围所有三角形的面积。
% 该函数的计算方法参考了论文 [1]。
% 
% 函数的输入参数包括三角网格的面片连接信息 face 和顶点坐标信息 vert，以及计算方式 type（可选）。
% 输出参数为每个顶点的面积 va。
% 
% 函数的具体实现如下：
% 
% 1. 如果没有输入计算方式 type，则默认为“一环”。
% 
% 2. 调用 face_area 函数计算每个面片的面积 fa。
% 
% 3. 根据计算方式 type 的不同，分别进行计算。
% 
%    a. 当 type 为“一环”时，首先调用 halfedge 函数计算每个半边的面积，然后使用 accumarray 函数
% 将每个顶点周围的面积相加得到 va。
% 
%    b. 当 type 为“混合”时，首先计算每个面片的三个角的余切值 c1、c2 和 c3，然后根据公式计算每个面片的
% 三个顶点的面积 vaf1、vaf2 和 vaf3。最后根据余切值的正负情况，将每个顶点的面积进行调整，
% 最终使用 accumarray 函数将每个顶点周围的面积相加得到 va。
% 
% 4. 如果输入的计算方式 type 不是“一环”或“混合”，则抛出错误信息。
% 
% 5. 函数 vector_cot 用于计算两个向量之间的余切值。
% 
% 该函数的作者为 Wen Cheng Feng，代码和文档的最后更新时间为 2014 年 3 月。
%% vertex area 
% Compute area around vertex, two variants available: "one_ring" and "mixed",
% see paper [1] for more details.
% 
% # Meyer, Mark, et al. "Discrete differential-geometry operators for 
%   triangulated 2-manifolds." Visualization and mathematics III. Springer 
%   Berlin Heidelberg, 2003. 35-57.
%
%% Syntax
%   va = vert_area(face,vertex)
%   va = vert_area(face,vertex,type)
%
%% Description
%  face: double array, nf x 3, connectivity of mesh
%  vert: double array, nv x 3, vertex of mesh
%  type: string, area type, either "one_ring", or "mixed"
% 
%  va: double array, nv x 1, area of all vertex.
% 
%% Example
%   va = vert_area(face,vertex)
%   va = vert_area(face,vertex,'one_ring') % same as last
%   va = vert_area(face,vertex,'mixed')

%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/03
%  Revised: 2014/03/28 by Wen, add code and doc
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function va = vert_area(face,vert,type)
if ~exist('type','var')
    type = 'one_ring';
end
fa = face_area(face,vert);

switch type
    case 'one_ring'
        [he,heif] = halfedge(face);
        va = accumarray(he(:,1),fa(heif));
    case 'mixed'
        dvf12 = vert(face(:,2),:)-vert(face(:,1),:);
        dvf23 = vert(face(:,3),:)-vert(face(:,2),:);
        dvf31 = vert(face(:,1),:)-vert(face(:,3),:);
        c1 = vector_cot(dvf12,-dvf31);
        c2 = vector_cot(dvf23,-dvf12);
        c3 = vector_cot(dvf31,-dvf23);
        vaf1 = (dot(dvf12,dvf12,2).*c3+dot(dvf31,dvf31,2).*c2)/8;
        vaf2 = (dot(dvf23,dvf23,2).*c1+dot(dvf12,dvf12,2).*c3)/8;
        vaf3 = (dot(dvf31,dvf31,2).*c2+dot(dvf23,dvf23,2).*c1)/8;
        ind1 = c1<0;
        vaf1(ind1) = fa(ind1)/2;
        vaf2(ind1) = fa(ind1)/4;
        vaf3(ind1) = fa(ind1)/4;
        ind2 = c2<0;
        vaf1(ind2) = fa(ind2)/4;
        vaf2(ind2) = fa(ind2)/2;
        vaf3(ind2) = fa(ind2)/4;
        ind3 = c3<0;
        vaf1(ind3) = fa(ind3)/4;
        vaf2(ind3) = fa(ind3)/4;
        vaf3(ind3) = fa(ind3)/2;
        va = accumarray(face(:),[vaf1;vaf2;vaf3]);
    otherwise
        error('unknown area type.')
end

function vc = vector_cot(v1,v2)
% cot of angle between two vectors
cs = dot(v1,v2,2);
d1 = dot(v1,v1,2);
d2 = dot(v2,v2,2);
vc = cs./sqrt(d1.*d2-cs.*cs);
