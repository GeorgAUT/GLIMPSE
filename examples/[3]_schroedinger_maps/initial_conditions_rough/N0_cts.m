function [out] = N0_cts(x)
%% Initial normal vector field

if -pi<=x && x<=-pi/2
    out=[-sin(2*x+pi),cos(2*x+pi),0];
elseif -pi/2<x && x<=0
    out=[0,1,0];    
elseif 0<x && x<=pi/2
    out=[-sin(2*x),cos(2*x),0];
elseif pi/2<x && x<=pi
    out=[0,-1,0];        
else 
    fprintf('invalid input for x.')
end

end

