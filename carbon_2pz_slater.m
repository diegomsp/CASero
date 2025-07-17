function orbital = carbon_2pz_slater(x, y, z, z_eff)
% Carbon 2pz slater orbital
%
% z_eff determines the effective nuclear charge interacting with the pz orbital
% Potential options:
%
% z_eff = 1
%     This corresponds to a hydrogen-like 2pz orbital and in
%     some cases matches well with DFT reference
%
% z_eff = 3.136
%     Value shown in https://en.wikipedia.org/wiki/Effective_nuclear_charge
%
% z_eff = 3.25
%     This is the value calculated by Slater's rules (https://en.wikipedia.org/wiki/Slater%27s_rules)
%     This value is also used in https://doi.org/10.1038/s41557-019-0316-8
%     This is the default.

if nargin < 4
    z_eff = 3.25;
end

r_grid = sqrt(x.^2 + y.^2 + z.^2); % angstrom
a0 = 0.529177; % Bohr radius in angstrom

orbital = z * exp(-z_eff * r_grid / (2 * a0));
end