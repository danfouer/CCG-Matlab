function holomorphic_form_basis = holomorphic_one_form(mesh, cf)

if ~exist('cf','var') || isempty(cf)
    cf = characteristic_one_form(mesh);
end

hf = cell(size(cf));
for i = 1:length(cf)
    hf{i} = harmonic_form(mesh, cf{i});
end

dual_hf = cell(size(cf));
for i = 1:length(dual_hf)
    dual_hf{i} = conjugate_harmonic_one_form(mesh, hf{i}, hf);
end

holomorphic_form_basis = cell(size(dual_hf));
for i = 1:length(holomorphic_form_basis)
    holomorphic_form_basis{i} = complex(hf{i}, dual_hf{i});
end
