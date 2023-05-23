% 这是一个用于计算网格面积的 MATLAB 函数。该函数的输入参数包括一个面数组 face 和一个顶点数组 vert，
% 输出参数为一个面积数组 fa。
% 
% 该函数的实现过程如下：
% 
% 1. 从面数组 face 中获取每个面的三个顶点的下标 fi、fj 和 fk。
% 
% 2. 从顶点数组 vert 中获取每个面的三个顶点的坐标 vi、vj 和 vk。
% 
% 3. 计算每个面的三条边的向量 vij、vjk 和 vki，以及每条边的长度 a、b 和 c。
% 
% 4. 计算每个面的半周长 s，然后根据海伦公式计算面积 fa。
% 
% 5. 将每个面的面积 fa 组成一个 nf x 1 的数组，作为函数的输出。
%% face area 
% Compute area of all face
%
%% Syntax
%   fa = face_area(face,vert)
%
%% Description
%  face: double array, nf x 3, connectivity of mesh
%  vert: double array, nv x 3, vertex of mesh
% 
%  fa: double array, nf x 1, area of all faces.
% 
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/03
%  Revised: 2014/03/23 by Wen, add doc
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function fa = face_area(face,vert)
fi = face(:,1);
fj = face(:,2);
fk = face(:,3);
vij = vert(fj,:)-vert(fi,:);
vjk = vert(fk,:)-vert(fj,:);
vki = vert(fi,:)-vert(fk,:);
a = sqrt(dot(vij,vij,2));
b = sqrt(dot(vjk,vjk,2));
c = sqrt(dot(vki,vki,2));
s = (a+b+c)/2.0;
fa = sqrt(s.*(s-a).*(s-b).*(s-c));
