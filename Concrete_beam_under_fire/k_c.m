function k_c=k_c(T,tipo)

% tipo 1 = silicoso
% tipo 2 = calcario

if (tipo == 1)
        k_c = 1.36-0.136*T/100+0.0057*(T/100)^2;
else if (tipo == 2)
        k_c = 2-0.2451*T/100+0.0107*(T/100)^2;
    end
end