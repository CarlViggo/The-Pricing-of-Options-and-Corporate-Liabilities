clear; clc; close all;

K = 1; 
Smax = 3; 
N = 200; 
M = 500;

sigma_vec = 0.1:0.05:0.5;
r_vec = 0.01:0.01:0.10;
T_vec = 0.25:0.25:3.0;

plot_EEP(r_vec, sigma_vec, @(r,s) solve_BS(s, r, 1.0, K,Smax,N,M), "r", "\sigma", "T=1.00", 1);
plot_EEP(T_vec, sigma_vec, @(T,s) solve_BS(s, 0.05, T, K,Smax,N,M), "T", "\sigma", "r=0.05", 2);
plot_EEP(T_vec, r_vec, @(T,r) solve_BS(0.2, r, T, K,Smax,N,M), "T", "r", "\sigma=0.2",3);

%%
function plot_EEP(xvec, yvec, f, xl, yl, param_str, num)

    EEP = zeros(length(yvec), length(xvec));

    for yi = 1:length(yvec)
        for xi = 1:length(xvec)
            EEP(yi,xi) = f(xvec(xi), yvec(yi));
        end
    end

    [X, Y] = meshgrid(xvec, yvec);
    figure(num);
    surf(X, Y, EEP, "EdgeColor","none");
    xlabel(xl); 
    ylabel(yl); 
    zlabel("(V_{amer}-V_{eur})/V_{eur}");
    title(sprintf("EEP relativ V_{eur} vid S=0.7K, %s", param_str));
    colorbar; 
    view(45,30);

end

function EEP = solve_BS(sigma, r, T, K, Smax, N, M)

    dS = Smax/N; 
    dt = T/M;
    S = (0:N)' * dS;
    i_K = round(0.7*K/dS) + 1;

    % terminalvillkor
    Ve = max(K-S,0); 
    Va = Ve;

    al = zeros(N+1,1); 
    be = zeros(N+1,1); 
    ga = zeros(N+1,1);

    for i = 2:N
        Si = S(i);
        al(i) = dt/4*(sigma^2*Si^2/dS^2 - r*Si/dS);
        be(i) = dt/2*(sigma^2*Si^2/dS^2 + r);
        ga(i) = dt/4*(sigma^2*Si^2/dS^2 + r*Si/dS);
    end

    A = diag(1+be(2:N)) - diag(ga(2:N-1),1) - diag(al(3:N),-1);
    B = diag(1-be(2:N)) + diag(ga(2:N-1),1) + diag(al(3:N),-1);
    al2 = al(2);

    for j = 1:M
        
        tau = T - j*dt;
        V0 = K*exp(-r*tau);

        % högerled
        rhs_e = B*Ve(2:N); % europeisk
        rhs_e(1) = rhs_e(1) + al2*V0; % europeisk
        rhs_a = B*Va(2:N); % amerikansk
        rhs_a(1) = rhs_a(1) + al2*V0; % amerikansk
        
        Ve(2:N) = A \ rhs_e;

        % amerikansk: optionen får aldrig vara 
        % värd mindre än lösenvärdet
        Va(2:N) = max(A \ rhs_a, K-S(2:N));

        % randpunkter
        Ve(1) = V0; Ve(N+1) = 0; % europeisk
        Va(1) = K; Va(N+1) = 0; % amerikansk 

    end

    EEP = (Va(i_K) - Ve(i_K)) / Ve(i_K);
end