% 这是一个用于计算网格中心点的 MATLAB 函数。该函数接受一个网格作为输入，返回中心点的索引。
% 
% 函数实现过程：
% 
% - 首先，通过 laplace_beltrami 函数计算网格的拉普拉斯-贝尔特拉米算子 L。
% - 然后，通过 eigs 函数计算 L 的前三个最小特征值和特征向量，其中第二个和第三个特征向量的交点即为网格的中心点。
% - 最后，返回中心点的索引 iv。
% 
% 该函数的作者没有给出详细的注释，但是根据代码可以大致理解其实现过程。
function iv = center(mesh)
% find the (rough) central vertex of input mesh, return index of central
% vertex. input mesh must be simply-connected with one boundary
L = laplace_beltrami(mesh);
% the smallest eigenvalue is 0, eigenvector is constant vector
% second smallest eigenvector (Fiedler vector)'s signs partition mesh into
% two approximately even parts, the cut is approximately the central
% foliation along that direction.
% third smallest eigenvector has similar property (as observed, not 
% confirmed), the cut foliates in orthogonal direction.
% thus the intersection of second and third smallest eigenvector's central
% foliation gives an rough center of mesh.
[V,~] = eigs(L,3,0);
[~,iv] = min(dot(V,V,2));
