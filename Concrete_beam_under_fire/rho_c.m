function rho_c=rho_c(rho0,T)

if (T>=20 && T<=115)
    rho_c=rho0;
else if (T>115 && T<=200)
        rho_c=rho0*(1-0.02*(T-115)/85);
    else if (T>200 && T<=400)
            rho_c=rho0*(0.98-0.03*(T-200)/200);
        else if (T>400)
                rho_c=rho0*(0.95-0.007*(T-400)/800);
            end
        end
    end
end