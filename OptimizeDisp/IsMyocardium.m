function ItIs = IsMyocardium(Xval, Rendo, Repi, Zbot, Ztop)

Rval = norm(Xval(1:2)); 

ItIs = 0;
tol = 1.0e-7;

if( Rval > Rendo - tol && Rval < Repi + tol )
    if (Xval(3) > Zbot && Xval(3) < Ztop)
        ItIs = 1;
    else
        Xval % This should not happen ... :)
    end
end


end

