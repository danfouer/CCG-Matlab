% 该函数的作用是计算环面到单位圆盘上的Riemann映射，假设最长的边界是外边界。
% 
% 输入参数：
% - mesh: 网格数据结构，包括face、edge、eif、halfedge、heif等。
% 
% 输出参数：
% - uv: 双精度数组，nv x 2的数组，表示环面上每个顶点在单位圆盘上的映射坐标。
% 
% 函数的实现过程：
% 1. 找到环面上的两个边界b1和b2，并将它们的值设为1和0，计算出一个调和函数f，使得f在b1和b2上的值分别为1和0。
% 2. 计算f的外微分df。
% 3. 找到b1和b2之间的最短路径cc，并将路径上的顶点标记为1，其余顶点标记为0。
% 4. 标记每个面的位置，indf等于1表示面在路径cc的一侧，等于-1表示面在路径cc的另一侧。
% 5. 定义一个函数g，使得g在路径cc上的值为0，在路径cc的另一侧的值为1，其余顶点的值为随机值。
% 6. 计算g的外微分w，并找到与w同一同调类的调和1-形式dh。
% 7. 计算df + 1i*dh，并对其进行积分得到eta。
% 8. 将eta映射到单位圆盘上得到z，将z的实部和虚部作为uv的输出。
% 
% 该函数的作者未知，版权归原作者所有。
function uv = annulus_riemann_map(mesh)
% compute riemann map from annulus surface to canonical annulus in unit
% disk, assume the longest boundary to be outer boundary
face = mesh.face;
edge = mesh.edge;
eif = mesh.eif;
halfedge = mesh.halfedge;
heif = mesh.heif;
nf = mesh.nf;
ne = mesh.ne;
nv = mesh.nv;
% compute harmonic function f with fixed value on two boundaries of annulus
% f(b1) = 1, f(b2) = 0
bds = boundary2(mesh);
b1 = bds{1};
b2 = bds{2};
if length(b1)<length(b2)
    b1 = bds{2};
    b2 = bds{1};
end
f = nan(mesh.nv,1);
f(b1) = 1;
f(b2) = 0;
f = harmonic_function(mesh, f);
df = exterior_derivative(mesh,f,0);
% find a shortest path between two boundaries b1,b2
cc = shortest_path(mesh,[b1(1),b2(1)]);
% prune path
ind = zeros(nv,1);
ind(b1) = 1;
ind(b2) = 2;
i = find(ind(cc)==1,1,'last');
j = find(ind(cc)==2,1,'first');
cc = cc(i:j);
% label face
% indf == 1 indicate faces on one side of path
% indf == -1 indicate faces on the other side
indf = zeros(nf,1);
hes = sparse(halfedge(:,1),halfedge(:,2),heif,nv,nv);
f1 = hes(cc(1:end-1)+(cc(2:end)-1)*nv);
f2 = hes(cc(2:end)+(cc(1:end-1)-1)*nv);
indf(f1) = 1;
indf(f2) = -1;
indf2 = zeros(ne,1);
ind1 = eif(:,1)>0;
indf2(ind1) = indf(eif(ind1,1));
ind2 = eif(:,2)>0;
indf2(ind2) = indf2(ind2) + indf(eif(ind2,2));
indf(eif(indf2==1&ind1,1)) = 1;
indf(eif(indf2==1&ind2,2)) = 1;
indf(eif(indf2==-1&ind1,1)) = -1;
indf(eif(indf2==-1&ind2,2)) = -1;
% define g on open surface, w = dg is closed one form and can be defined on
% original closed surface
g = rand(nv,1);
g(cc) = 0;
w = exterior_derivative(mesh,g,0);
% g == 1 on the other side, adjust w
g(cc) = 1;
inde = zeros(ne,1);
inde(ind1) = inde(ind1) + indf(eif(ind1,1));
inde(ind2) = inde(ind2) + indf(eif(ind2,2));
inde = inde>0;
w(inde) = g(edge(inde,2)) - g(edge(inde,1));
% w is closed one form, find harmonic one form dh in the homological class
% of w
dh = harmonic_form(mesh, w);
% integrate df + 1i*dh
des = sparse(edge(:,1),edge(:,2),dh,nv,nv);
des = des - conj(des');
k = sum(full(des(b1+(b1([2:end,1])-1)*nv)));
dh = dh/k;
deta = df + 1i*dh;
% eta has period 1
eta = integration(mesh,deta);
% exponential map
z = exp(2*pi*eta);
uv = [real(z),imag(z)];
