% 这是一个用于计算网格面法向量的 MATLAB 函数。该函数的输入参数包括一个面数组 face 和一个顶点数组 vert，
% 输出参数为一个法向量数组 normal。
% 
% 该函数的实现过程如下：
% 
% 1. 从面数组 face 中获取每个面的三个顶点的下标，分别为 fi、fj 和 fk。
% 
% 2. 从顶点数组 vert 中获取每个面的三个顶点的坐标，分别为 vi、vj 和 vk。
% 
% 3. 计算每个面的两个边的向量 vij 和 vik，然后计算它们的叉积，得到面的法向量 normal。
% 
% 4. 对于每个法向量，使用 bsxfun 函数将其除以其模长，得到单位法向量。
% 
% 5. 将每个面的法向量组成一个 nf x 3 的数组，作为函数的输出。
function normal = face_normal(face, vert)
normal = cross(vert(face(:, 2), :) - vert(face(:, 1), :), ...
    vert(face(:, 3), :) - vert(face(:, 1), :));
normal = bsxfun(@rdivide,normal,sqrt(sum(normal.^2,2)));
end