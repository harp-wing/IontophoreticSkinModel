% Physical Parameters of 4-Layer Model
D_sc = 1.49e-4;     % cm^2/h
D_rs = 1.80e-3;     % cm^2/h
P1 = 4.84;          % SC/Donor partition coefficient
P2 = 0.96;          % RS/SC partition coefficient
P3 = 0.111;         % RC/RS partition coefficient
K3 = 0.043;         % RS-RC interface mass transfer coefficient

% Exact Geometry from Jonsdottir et al.
A = 0.6362;         % Exposure area (cm^2)
L = [1.572, 0.002, 0.098, 18.263]; % D, SC, RS, RC thicknesses (cm)
M_total = 5.0;      % mg

% Diclofenac Properties
pKa = 4.2;
molW = 296.1; % g/mol

% Assume slightly colder than body temp
% (https://pmc.ncbi.nlm.nih.gov/articles/PMC7992731/)
T = 34 + 273.15; % K

% Electrical Parameters
z = -1; % Approximately -1 because pKa (4.2) << pH (7.4)
I = 0.5e-3; % A/cm^2
R_sc = 45000; % ohm
R_derm = 500; % ohm
R_total = R_sc + R_derm;

% Effective skin thickness (SC + RS)
L_skin = L(2) + L(3);
% Voltage from Ohm's law
V_sc = I * R_sc; % V
V_rs = I * R_derm;
% Electric fields
E(1) = V_sc / L(2); % V/cm
E(2) = V_rs / L(3);

v_eo = 0.02; % cm/h

% Numerical Setup
N = [1, 4, 196, 1]; % Nodes per layer
dx = L ./ N;
total_nodes = sum(N);
C0 = zeros(total_nodes, 1);
C0(1:N(1)) = M_total / (A * L(1)); % Initial dose in donor

% Solving
% t_select = [1, 2, 5, 21, 41, 101, 207]; % Chosen time points in the paper
% t_display = (t_select-1)./4;
% Have to rescale because iontophoresis speeds up transport so significantly.
t_select = [0, 500, 1000, 1500, 2000] + 1;
t_display = (t_select-1)./1000;
% t_span = linspace(0, 51.5, 207);
t_span = linspace(0, 2, 2001);
[t, C] = ode15s(@(t, c) odefun(t, c, D_sc, D_rs, P1, P2, P3, K3, N, dx, L, z, T, E, v_eo), t_span, C0);

% Plotting
close all;
figure(1); hold on;
colors = lines(length(t_span));

% x coordinates for SC and RS only (skin layers)
x_sc_start = 0;
x_sc_end   = L(2);
x_rs_end   = L(2) + L(3);

% Node centers for SC and RS
x_sc = x_sc_start + dx(2)*((1:N(2)) - 0.5);
x_rs = x_sc_end   + dx(3)*((1:N(3)) - 0.5);
x_skin = [x_sc, x_rs];  % skin depth > 0

for i = 2:length(t_select)
    col = colors(i-1, :);

    % Concatenate all x and C vectors
    x_all = [-L(1), 0, x_skin, x_rs_end, x_rs_end + L(4)];
    C_all = [C(t_select(i),1), C(t_select(i),1), C(t_select(i), 2:end-1), C(t_select(i),end), C(t_select(i),end)];

    plot(x_all, C_all, 'Color', col, 'LineWidth', 1.5, 'DisplayName', sprintf('t = %g h', t_display(i)));
end

xlim([-0.019, 0.12])
xlabel('Skin Depth [cm]'); ylabel('Average C [mg/cm^3]');
title('Diff + EM + EO (v_{eo} = 0.02 cm/hr)');
legend('Location', 'northwest'); grid on;

figure(2);
plot(t, C(:, end), 'k-', 'LineWidth', 1.5);
% xlim([0, 51.5]);
xlabel('Time (h)');
ylabel('Avg. Conc. in Receiver (mg/cm^3)');
title('Diff + EM + EO (v_{eo} = 0.02 cm/hr)');
grid on;

function dCdt = odefun(~, C, D_sc, D_rs, P1, P2, P3, K3, N, dx, L, z, T, E, v_eo)
    % Constants
    F = 96485.3321; % s*A/mol
    R = 8.314462618; % g*cm^2/(s^2*K*mol)

    N_sc = N(2); % Number of nodes in SC
    N_rs = N(3); % NUmber of nodes in RS

    dx_sc = dx(2); % SC node width
    dx_rs = dx(3); % RS node width

    E_sc = E(1);
    E_rs = E(2);

    L_d   = L(1); % thickness of donor compartment
    L_rc  = L(4); % thickness of receiver

    dCdt = zeros(length(C), 1); % Initialize output var

    % Assign nodes to the compartments
    C_d  = C(1);
    C_sc = C(2:N_sc+1);
    C_rs = C(N_sc+2:end-1);
    C_r  = C(end);

    % Calculate electromigration velocities (cm/h)
    % Should be always positive as the drug is placed at the same charge electrode
    v_em_sc = abs(z) * D_sc * F * E_sc / (R*T);
    v_em_rs = abs(z) * D_rs * F * E_rs / (R*T);
    
    % Total velocity = electromigration + electroosmosis
    % In the case of anion delivery from the cathode, EO will oppose
    % transport, while diffusion and EM will assist transport.
    if z > 0
        v_sc = v_em_sc + v_eo;
        v_rs = v_em_rs + v_eo;
    else
        v_sc = v_em_sc - v_eo;
        v_rs = v_em_rs - v_eo;
    end

    % Pe_sc = v_sc*dx_sc/D_sc
    % Pe_rs = v_rs*dx_rs/D_rs

    % Interface fluxes
    % Advenction-Diffusion: J = -D*dC/dx + v*C
    
    % Donor to SC (x = 0)
    % C_sc(1) = P1 * C_d

    J12_diff = D_sc * (P1 * C_d - C_sc(1)) / (dx_sc / 2);

    % Use upstream node to calculate advective flux
    if v_sc >= 0
        J12_adv = v_sc * P1 * C_d;
    else
        J12_adv = v_sc * C_sc(1);
    end

    J12 = J12_diff + J12_adv;

    
    % SC to RS
    % C_rs(1) = P2 * C_sc(end)

    % Harmonic mean diffusivity between compartments
    D_eff_23 = 2 * D_sc * D_rs / (D_sc * dx_rs + D_rs * dx_sc);
    J23_diff = D_eff_23 * (P2 * C_sc(end) - C_rs(1));

    if v_rs >= 0
        J23_adv = v_rs * P2 * C_sc(end);
    else
        J23_adv = v_rs * C_rs(1);
    end

    J23 = J23_diff + J23_adv;

    % RS to RC: diffusive kinetic resistance + advection
    if v_rs >= 0
        J34 = K3 * (C_rs(end) - C_r / P3) + v_rs * P3 * C_rs(end);
    else
        J34 = K3 * (C_rs(end) - C_r / P3) + v_rs * C_r;
    end

    % Lumped Compartments
    dCdt(1)   = -J12 / L_d;
    dCdt(end) =  J34 / L_rc;
    
    % Debug
    % fprintf('E_sc = %.4f V/cm,  E_rs = %.4f V/cm\n', E_sc, E_rs)
    % fprintf('v_em_sc = %.4f cm/h,  v_em_rs = %.4f cm/h\n', v_em_sc, v_em_rs)
    % fprintf('v_sc = %.4f cm/h,  v_rs = %.4f cm/h\n', v_sc, v_rs)
    % fprintf('Pe_sc = %.2f,  Pe_rs = %.2f\n', v_sc*dx_sc/D_sc, v_rs*dx_rs/D_rs)
    % fprintf('C_d(0) = %.4f\n', C_d)
    % fprintf('J12 = %.6f,  J23 = %.6f,  J34 = %.6f\n', J12, J23, J34)
    % fprintf('dCdt(1) = %.6f,  dCdt(end) = %.6f\n', dCdt(1), dCdt(end))

    %% SC
    % Left boundary node (receives flux from donor)
    J_sc_interior_1 = -D_sc*(C_sc(2)-C_sc(1))/dx_sc + v_sc*upwind(C_sc(1), C_sc(2), v_sc);
    dCdt(2) = (J12 - J_sc_interior_1) / dx_sc;

    % Interior SC nodes
    for k = 2 : N_sc-1
        J_left  = -D_sc*(C_sc(k) - C_sc(k-1))/dx_sc + v_sc*upwind(C_sc(k-1), C_sc(k),   v_sc);
        J_right = -D_sc*(C_sc(k+1) - C_sc(k))  /dx_sc + v_sc*upwind(C_sc(k),   C_sc(k+1), v_sc);
        dCdt(k+1) = (J_left - J_right) / dx_sc;
    end

    % Right boundary node of SC (loses flux to RS)
    J_sc_interior_end = -D_sc*(C_sc(end)-C_sc(end-1))/dx_sc + v_sc*upwind(C_sc(end-1), C_sc(end), v_sc);
    dCdt(N_sc+1) = (J_sc_interior_end - J23) / dx_sc;

    %% RS
    % Left boundary node (receives flux from SC)
    J_rs_interior_1 = -D_rs*(C_rs(2)-C_rs(1))/dx_rs + v_rs*upwind(C_rs(1), C_rs(2), v_rs);
    dCdt(N_sc+2) = (J23 - J_rs_interior_1) / dx_rs;

    % Interior RS nodes
    for k = 2 : N_rs-1
        J_left  = -D_rs*(C_rs(k)   - C_rs(k-1))/dx_rs + v_rs*upwind(C_rs(k-1), C_rs(k),   v_rs);
        J_right = -D_rs*(C_rs(k+1) - C_rs(k))  /dx_rs + v_rs*upwind(C_rs(k),   C_rs(k+1), v_rs);
        dCdt(N_sc+1+k) = (J_left - J_right) / dx_rs;
    end

    % Right boundary node of RS (loses flux to RC)
    J_rs_interior_end = -D_rs*(C_rs(end) - C_rs(end-1))/dx_rs + v_rs*upwind(C_rs(end-1), C_rs(end), v_rs);
    dCdt(N_sc+N_rs+1) = (J_rs_interior_end - J34) / dx_rs;
end

function C_up = upwind(C_left, C_right, v)
% Returns the upwind concentration for the advective flux v*C at an interface.
    if v >= 0
        C_up = C_left;
    else
        C_up = C_right;
    end
end