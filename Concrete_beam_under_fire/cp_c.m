function cp_c=cp_c(T)

if (T>=20 && T<=100)
    cp_c=900;
else if (T>100 && T<=200)
        cp_c=900+(T-100);
    else if (T>200 && T<=400)
            cp_c=1000+(T-200)/2;
        else if (T>400)
                cp_c=1100;
            end
        end
    end
end