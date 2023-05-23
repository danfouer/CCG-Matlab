% 函数的作用是将一个三角网格映射到单位球面上。该函数使用了球面谐函数映射的方法，
% 将三角网格映射到单位圆盘上，然后再将圆盘映射到单位球面上。
% 
% 函数的输入参数有两个，分别是面片(face)和顶点(vertex)。其中，face是一个nf x 3的双精度数组，
% 表示网格的连接关系；vert是一个nv x 3的双精度数组，表示网格的顶点坐标。
% 
% 函数的输出参数是一个nv x 3的双精度数组uvw，表示映射后的顶点坐标。其中，uvw的前两列表示映射到单位圆盘上的坐标，
% 第三列表示映射到单位球面上的坐标。
% 
% 函数的具体实现过程如下：
% 
% 1. 计算网格的拉普拉斯-贝尔特拉米算子A。
% 
% 2. 计算A的前两个特征向量，将第二个特征向量fd作为初始映射函数。
% 
% 3. 将fd的正值部分映射到单位圆盘上，并计算该部分的边界bd1。
% 
% 4. 将整个网格映射到单位圆盘上，并使用bd1的映射结果作为边界条件，计算圆盘上的调和映射uv。
% 
% 5. 将uv映射到单位球面上，得到uvw的前两列。
% 
% 6. 计算uvw的第三列，即将前两列的坐标映射到球面上，并计算其长度的平方减1。
% 
% 7. 将uvw的长度归一化，得到单位球面上的坐标。
% 
% 8. 对于不在正值部分的顶点，将其第三列的坐标取相反数。
% 
% 9. 对于每个顶点，计算其邻域内所有顶点的平均值，并将该平均值的长度归一化，得到该顶点的新坐标。
% 
% 10. 重复步骤9，直到收敛或达到最大迭代次数。
% 
% 11. 返回uvw。
function uvw = spherical_harmonic_map(face,vert)
nv = size(vert,1);
A = laplace_beltrami(face,vert);
[V,D] = eigs(A,2,0);
fd = V(:,2);
ind = fd>0;
indf = sum(ind(face),2)>=2;
f1 = face(indf,:);
ind = false(nv,1);
ind(f1) = true;
uv1 = disk_harmonic_map(f1,vert);
bd1 = compute_bd(f1);
uv = disk_harmonic_map(face,vert,bd1,uv1(bd1,:));
uvw = zeros(nv,3);
uvw(:,1:2) = uv*2;
uvw(:,3) = dot(uv,uv,2)-1;
uvw = uvw./(dot(uv,uv,2)+1);
uvw(~ind,3) = -uvw(~ind,3);

vr = compute_vertex_ring(face,uvw,bd1);
for i = 1:length(vr)
    ri = vr{i};
    x = mean(uvw(ri,:));
    uvw(bd1(i),:) = x./norm(x);
end

dt = 0.2;
k = 0;
while k <= 3000
    d = A*uvw;
    dd = d-dot(d,uvw,2).*uvw;
    if mod(k,100) == 0
        fprintf('#%d |df/dt| = %.10f\n',k,norm(dd));
    end
    if norm(dd) < eps
        break;
    end
    while true
        uvw1 = uvw + dt*(dd);
        l = sqrt(dot(uvw1,uvw1,2));
        uvw1 = uvw1./l;
        d = A*uvw1;
        dd1 = d-dot(d,uvw1,2).*uvw1;
        if norm(dd1)>norm(dd)
            dt = dt/2;
        else
            break;
        end
    end
    uvw = uvw + dt*(dd);
    l = sqrt(dot(uvw,uvw,2));
    uvw = uvw./l;
    k = k+1;
end
