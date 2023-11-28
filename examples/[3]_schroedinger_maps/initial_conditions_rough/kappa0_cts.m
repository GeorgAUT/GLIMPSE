function [out] = kappa0_cts(x)
%% Initial curvature

if -pi<=x && x<=-pi/2
    out=2;
elseif -pi/2<x && x<=0
    out=0;    
elseif 0<x && x<=pi/2
    out=2;
elseif pi/2<x && x<=pi
    out=0;        
else 
    fprintf('invalid input for x.')
end


end

