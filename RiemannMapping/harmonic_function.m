% 该函数的作用是给定边界上的初始值（由f0指定），计算相应的调和函数。
% 
% 输入参数：
% - mesh: 网格数据结构，包括laplace_beltrami等。
% - f0: 双精度数组，nv x k的数组，表示在k个边界上的初始值。
% 
% 输出参数：
% - f: 双精度数组，nv x k的数组，表示计算得到的调和函数。
% 
% 函数的实现过程：
% 1. 计算拉普拉斯-贝尔特拉米算子L。
% 2. 对于每个边界，找到边界上的顶点并将它们的值设为NaN或Inf，将L中对应的行和列删除得到L2。
% 3. 将L2和f0(~ind,i)代入L2\b中求解得到f(ind,i)，即在边界ind上的调和函数。
% 
% 该函数的作者未知，版权归原作者所有。
function f = harmonic_function(mesh, f0)
% given initial value on boundary (specified by f0), compute corresponding
% harmonic function

L = laplace_beltrami(mesh);
f = f0;
for i = 1:size(f0,2)
    ind = isnan(f0(:,i)) | isinf(f0(:,i));
    L2 = L(ind,ind);
    b = -L(ind,~ind)*f0(~ind,i);
    f(ind,i) = L2\b;
end
