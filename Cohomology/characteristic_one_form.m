% 该函数用于计算特征一形式。
% 
% 输入参数：
% - mesh: 包含网格信息的结构体。
% - hb: cell类型，表示基础同调群的基。如果未定义，则使用homology_basis函数计算基。
% 
% 输出参数：
% - cf: cell类型，表示特征一形式。
% 
% 函数的实现过程：
% 1. 对于每个基hb{i}，使用slice_mesh函数将网格沿着hb{i}切开。
% 2. 构造一个函数f2，其中沿着hb{i}的一侧为0，另一侧为1。
% 3. 使用exterior_derivative函数计算df2，它是沿着切开的网格的一形式。
% 4. 将df2映射回原始网格，得到df。
% 5. 将df存储在cf{i}中。
function cf = characteristic_one_form(mesh, hb)
if ~exist('hb','var') || isempty(hb)
    hb = homology_basis(mesh);
end
cf = cell(size(hb));
for i = 1:length(hb)
    % slice mesh open along basis hb{i}
    bi = hb{i};
    if bi(1) == bi(end)
        ee = [bi,bi([2:end,1])];
    else
        ee = [bi(1:end-1),bi(2:end)];
    end
    mesh2 = slice_mesh(mesh,ee);
    % construct f
    f2 = ones(mesh2.nv,1);
    bd2 = boundary2(mesh2);
    f2(bd2{1}) = 0;
    f2(bd2{2}) = 1;
    % df is closed on original surface
    df2 = exterior_derivative(mesh2,f2,0);
    % map df back to original surface
    edge = mesh2.father(mesh2.edge);
    F = sparse(edge(:,1),edge(:,2),df2,mesh.nv,mesh.nv);
    F = F-F';
    df = F(mesh.edge(:,1) + mesh.nv*(mesh.edge(:,2)-1));
    cf{i} = full(df);
end
