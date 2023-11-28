function [out] = T0_cts(x)
%% Initial tangent vector field

if -pi<=x && x<=-pi/2
    out=[cos(2*x+pi),sin(2*x+pi),0];
elseif -pi/2<x && x<=0
    out=[1,0,0];    
elseif 0<x && x<=pi/2
    out=[cos(2*x),sin(2*x),0];
elseif pi/2<x && x<=pi
    out=[-1,0,0];        
else 
    fprintf('invalid input for x.')
end

end

