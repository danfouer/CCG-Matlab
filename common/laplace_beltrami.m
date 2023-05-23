% 这是一个计算网格上Laplace Beltrami算子的MATLAB函数。Laplace Beltrami算子是微分几何中的一个概念，
% 用于描述曲面上的函数的变化率。该函数使用余切公式来计算Laplace Beltrami算子，其中有三种变体：'Polthier'、
% 'Meyer'和'Desbrun'。默认使用'Polthier'方法。
% 
% 函数的输入参数为一个网格，包括三个部分：face、edge和vert。face是一个nf x 3的矩阵，
% 表示网格的三角面片的连接关系；edge是一个ne x 2的矩阵，表示网格的边的连接关系；
% vert是一个nv x 3的矩阵，表示网格的顶点坐标。函数的第二个输入参数是一个字符串，
% 表示使用的余切公式的方法，可以是'Polthier'、'Meyer'或'Desbrun'。
% 函数的输出是一个nv x nv的稀疏矩阵，表示Laplace Beltrami算子。
% 
% 函数的具体实现如下：
% 
% 1. 如果没有输入第二个参数，则默认使用'Polthier'方法。
% 
% 2. 根据输入的网格，计算边的权重。如果网格中已经包含了边的权重信息，则直接使用该信息，
% 否则调用edge_weight函数计算边的权重。
% 
% 3. 根据选择的余切公式的方法，计算边的权重的加权平均值，并构造出Laplace Beltrami算子的稀疏矩阵。
% 
% 4. 对于每个顶点，将其对应的行减去该行的和，得到最终的Laplace Beltrami算子。
% 
% 函数中使用了vertex_area函数来计算顶点的面积，该函数根据选择的方法不同，计算的顶点面积也不同。
% 其中，'mixed'方法计算的是顶点周围所有三角形的面积的平均值，
% 'one_ring'方法计算的是顶点周围所有三角形的面积之和。
%% laplace_beltrami 
% Laplace Beltrami operator on the mesh.
% 
% Cotangent formula is used, while there are some variants:
% 
% * 'Polthier', see paper [1]
% * 'Meyer', see paper [2]
% * 'Desbrun', see paper [3]
% 
% For comparison and convergence analysis, see paper [4]
% 
% # K. Polthier. Computational Aspects of Discrete Minimal Surfaces. In 
%   Proc. of the Clay Summer School on Global Theory of Minimal Surfaces, 
%   J. Hass, D. Hoffman, A. Jaffe, H. Rosenberg, R. Schoen, M. Wolf (Eds.), 
%   to appear, 2002.
% # M. Meyer, M. Desbrun, P. Schröder, and A. Barr. Discrete 
%   Differential-Geometry Operator for Triangulated 2-manifolds. In Proc. 
%   VisMath'02, Berlin, Germany, 2002.
% # M. Desbrun, M. Meyer, P. Schröder, and A. H. Barr. Implicit Fairing 
%   of Irregular Meshes using Diffusion and Curvature Flow. SIGGRAPH99, 
%   pages 317-324, 1999.
% # Xu, Guoliang. "Convergent discrete laplace-beltrami operators over 
%   triangular surfaces." Geometric Modeling and Processing, 2004. 
%   Proceedings. IEEE, 2004.
%
%% Syntax
%   A = laplace_beltrami(mesh)
%   A = laplace_beltrami(mesh, method)
%
%% Description
%  face  : double array, nf x 3, connectivity of mesh
%  vertex: double array, nv x 3, vertex of mesh
%  method: string, optional, method of cotangent formula, can be one of
%          three: 'Polthier', 'Meyer', 'Desbrun'. Default is 'Polthier'.
% 
%  A: sparse matrix, nv x nv, Laplace Beltrami operator
%
%% Example
%   A = laplace_beltrami(mesh)
%   A = laplace_beltrami(mesh,'Polthier') % same as last 
%   A = laplace_beltrami(mesh,'Meyer')
%   A = laplace_beltrami(mesh,'Desbrun')
%
%% Contribution
%  Author : Wen Cheng Feng
%  Created: 2014/03/03
%  Revised: 2014/03/03 by Wen, add more cotangent formula variants, not
%           implemented
%  Revised: 2014/03/23 by Wen, add doc
%  Revised: 2014/03/28 by Wen, add code for 'Meyer' and 'Desbrun' methods
% 
%  Copyright 2014 Computational Geometry Group
%  Department of Mathematics, CUHK
%  http://www.math.cuhk.edu.hk/~lmlui

function A = laplace_beltrami(mesh, method)
% default method is 'Polthier'
if nargin == 1
    method = 'Polthier';
end
face = mesh.face;
edge = mesh.edge;
vert = mesh.vert;
if ~isfield(mesh,'ew')
    ew = edge_weight(mesh);
else
    ew = mesh.ew;
end

switch method
    case 'Polthier'
        A = sparse([edge(:,1);edge(:,2)],[edge(:,2);edge(:,1)],[ew;ew]);
        sA = sum(A,2);
        A = A - diag(sA);
    case 'Meyer'
        va = vertex_area(face,vert,'mixed');
        ew = (ew./va(edge(:,1))+ew./va(edge(:,2)))/2;
        A = sparse([edge(:,1);edge(:,2)],[edge(:,2);edge(:,1)],[ew;ew]);
        sA = sum(A,2);
        A = A - diag(sA);
    case 'Desbrun'
        va = vertex_area(face,vert,'one_ring');
        ew = (ew./va(edge(:,1))+ew./va(edge(:,2)))/2*3;
        A = sparse([edge(:,1);edge(:,2)],[edge(:,2);edge(:,1)],[ew;ew]);
        sA = sum(A,2);
        A = A - diag(sA);
    otherwise
        error('Wrong method. Available methods are: Polthier, Meyer, Desbrun')
end
