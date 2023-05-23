% 这是一个用于计算半边数据结构的 MATLAB 函数。该函数的输入参数为一个面数组 face，输出参数包括一个半边数组 he
% 和一个面下标数组 heif。
% 
% 该函数的实现过程如下：
% 
% 1. 从面数组 face 中获取每个面的三个顶点的下标，分别为 fi、fj 和 fk。
% 
% 2. 将每个面的三个顶点的下标按照顺序组成三个半边，分别为 (fi,fj)、(fj,fk) 和 (fk,fi)。
% 
% 3. 将所有半边组成一个 nf x 3 的数组 he，其中每行表示一个半边，包括起点和终点的下标。
% 
% 4. 将每个半边所在的面的下标组成一个 nf x 3 的数组 heif，其中每行表示一个半边所在的面的下标。
% 
% 5. 将 he 和 heif 组成一个元组作为函数的输出。
%% compute halfedge 
% Halfedge is simply directed edge, each face has three halfedges.
% This function will return all nf x 3 halfedges, as well as a nf x 3
% vector indicate which face the halfedge belongs to.
%
%% Syntax
%   [he,heif] = halfedges(face)
%
%% Description
%  face: double array, nf x 3, connectivity of mesh
% 
%  he  : double array, (nf x 3) x 2, each row is a halfedge
%  heif: double array, (nf x 3) x 1, face id in which the halfedge lies in
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/06
%  Revised: 2014/03/23 by Wen, add doc
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function [he,heif] = halfedges(face)
nf = size(face,1);
he = [reshape(face',nf*3,1),reshape(face(:,[2 3 1])',nf*3,1)];
heif = reshape(repmat(1:nf,[3,1]),nf*3,1);
