% 本 MATLAB 代码实现了从一个单连通曲面到单位圆盘的 Riemann 映射的计算。该函数的输入参数为一个网格结构体
% `mesh`，输出参数为一个复平面上的点的数组 `uv`，表示将 `mesh` 映射到单位圆盘上的结果。
% 
% 该函数的实现过程如下：
% 
% 1. 通过 `center` 函数找到一个中心三角形，并将其从 `mesh` 中移除，得到一个拓扑环面。
% 这个中心三角形可以通过找到一个顶点，然后找到以该顶点为公共顶点的三个面，这三个面构成的三角形即为中心三角形。
% 
% 2. 通过 `vert_face_ring` 函数找到中心三角形的面环，将这个面从 `mesh` 中移除，得到一个拓扑环面。
% 这个函数返回一个单元格数组，其中第一个单元格是中心三角形的面环。
% 
% 3. 通过 `make_mesh` 函数将新的面列表 `face2` 和原始的顶点列表 `mesh.vert` 组合成一个新的网格结构体 `mesh2`。
% 
% 4. 通过 `annulus_riemann_map` 函数计算 `mesh2` 到单位圆盘的 Riemann 映射。
% 
% 其中，`center` 函数用于找到中心三角形，`vert_face_ring` 函数用于找到中心三角形的面环，
% `make_mesh` 函数用于构造新的网格结构体，`annulus_riemann_map` 函数用于计算 Riemann 映射。
% 
% 需要注意的是，该函数的实现过程中，使用了其他函数的实现，这些函数的实现过程可以参考其他代码的解释。
function uv = riemann_map(mesh)
% compute riemann map from simply-connected surface to unit disk

% first find a central triangle and remove it, remaining mesh becomes 
% topological annulus, with same vertex set
iv = center(mesh);
vfr = vert_face_ring(mesh,iv);
ivf = vfr{1};
ind = true(mesh.nf,1);
ind(ivf(1)) = false;
face2 = mesh.face(ind,:);
mesh2 = make_mesh(face2,mesh.vert);
uv = annulus_riemann_map(mesh2);
