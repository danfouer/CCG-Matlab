% 该函数用于计算给定顶点或所有顶点的一环邻居，可以选择是否按逆时针顺序排序。默认情况下不排序。
% 
% 输入参数：
% - face: double类型，nf x 3，表示网格的连接性。
% - vc: double类型，n x 1或1 x n，表示顶点集合，可以为空，相当于所有顶点。
% - ordered: bool类型，标量，表示是否需要逆时针排序。
% 
% 输出参数：
% - vr: cell类型，nv x 1，每个单元格是一个一环邻居顶点，是一个double类型的数组。
% 
% 函数的使用方法：
% - 计算所有顶点的一环邻居，不排序：vr = vert_vert_ring(face)
% - 计算顶点1:100的一环邻居，不排序：vr = vert_vert_ring(face,1:100,false)
% - 计算顶点1:100的一环邻居，按逆时针排序：vr = vert_vert_ring(face,1:100,true)
% - 计算所有顶点的一环邻居，按逆时针排序（可能很慢）：vr = vert_vert_ring(face,[],ordered)
% - 与上一行代码等价：vr = vert_vert_ring(face,1:nv,ordered)
% 
% 函数的实现过程：
% 1. 计算顶点数nv，假设面的编号从1开始，并按顺序编号。
% 2. 如果vc为空，则vc为1:nv。
% 3. 如果ordered未定义，则ordered为false。
% 4. 初始化vr为一个大小为vc的cell数组。
% 5. 计算边界bds，isbd表示每个顶点是否在边界上。
% 6. 如果不排序，则计算邻接矩阵am，然后找到vc中每个顶点的邻居，将其存储在vr中。
% 7. 如果需要排序，则使用三角剖分DT和vertexAttachments函数计算每个顶点的邻居，然后按逆时针顺序排序。
% 8. 如果需要排序，但是使用vertexAttachments函数计算的排序速度太慢，则使用connectivity函数计算排序。
%% vert_vert_ring 
% Compute one-ring neighbor of given vertex or all vertex, with or 
% without ccw order. Default is no order.
% 
% In some algorithms, ordered one-ring neighbor is necessary. However,
% compute ordered one-ring is significantly slower than unordered one, so
% compute with order only absolutely necessary.
%
%% Syntax
%   vr = vert_vert_ring(face)
%   vr = vert_vert_ring(face,vc)
%   vr = vert_vert_ring(face,vc,ordered)
%
%% Description
%  face: double array, nf x 3, connectivity of mesh
%  vc  : double array, n x 1 or 1 x n, vertex collection, can be empty, 
%        which equivalent to all vertex.
%  ordered: bool, scaler, indicate if ccw order needed.
% 
%  vr: cell array, nv x 1, each cell is one ring neighbor vertex, which is
%      a double array
%
%% Example
%   % compute one ring of all vertex, without order
%   vr = vert_vert_ring(face)
% 
%   % compute one ring of vertex 1:100, without ccw order
%   vr = vert_vert_ring(face,1:100,false)
% 
%   % compute one ring of vertex 1:100, with ccw order
%   vr = vert_vert_ring(face,1:100,true)
% 
%   % compute one ring of all vertex, with ccw order (may be slow)
%   vr = vert_vert_ring(face,[],ordered)
% 
%   % same with last one
%   vr = vert_vert_ring(face,1:nv,ordered)
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/06
%  Revised: 2014/03/23 by Wen, add doc
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function vr = vert_vert_ring(mesh,vc,ordered)
% number of vertex, assume face are numbered from 1, and in consecutive
% order
nv = mesh.nv;
if ~exist('vc','var') || isempty(vc)
    vc = (1:nv)';
end
if ~exist('ordered','var')
    ordered = false;
end
vr = cell(size(vc));
bds = boundary(mesh.face);
isbd = false(nv,1);
if iscell(bds)
    for i = 1:length(bds)
        isbd(bds{i}) = true;
    end
else
    isbd(bds) = true;
end
if ~ordered
    [am,~] = adjacency_matrix(mesh.face);
    [I,J,~] = find(am(:,vc));
    vr = cell(size(vc,1),1);
    for i = 1:length(J)
        vr{J(i)}(end+1) = I(i);
    end
end
if ordered
    DT = triangulation(mesh.face,mesh.vert);
    va = vertexAttachments(DT,vc);
    for i = 1:size(vc,1)                
        vai = va{i};        
        fai = mesh.face(vai,:);
        ind = mod(find(fai'==vc(i)),3)-1;
        if isempty(ind)
            vr{i} = [];
            continue;
        end
        ind1 = ind(1)+2;
        ind(ind<=0) = ind(ind<=0)+3;
        vri = zeros(1,length(ind));
        for j = 1:length(ind)
            vri(j) = fai(j,ind(j));
        end
        %if (length(unique(fai))-1)*2 == length(fai(:))-length(ind)
        if ~isbd(vc(i))
            vri(end+1) = vri(1);
        else
            vri = [fai(1,ind1),vri];
        end
        vr{i} = vri;
    end
end
if ordered && false
    [vvif,nvif,pvif] = connectivity(mesh.face);
    for i = 1:size(vc,1)
        fs = vvif(vc(i),:);
        v1 = full(find(fs,1,'first'));
        if isbd(vc(i))
            while vvif(v1,vc(i))
                f2 = full(vvif(v1,vc(i)));
                v1 = full(pvif(f2,v1));
            end
        end
        vi = v1;
        v0 = v1;
        while vvif(vc(i),v1)
            f1 = full(vvif(vc(i),v1));
            v1 = full(nvif(f1,v1));
            vi = [vi,v1];
            if v0 == v1
                break;
            end
        end
        vr{i} = vi;
    end
end