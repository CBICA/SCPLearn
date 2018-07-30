%> @file  ProjectionOnUnitBoxSimplex.m
%> @brief Function that projects on the feasible set Boxed-Sparsity
%======================================================================
%> @brief It takes as input a vector value and projects it into a "box"
%> (box feasible set intersected with l-1 norm inequality)
%> @param x value of variable
%> @param lambda l-1 norm sparsity value
%> @retval xp projected value of variable
%> @b Author: 
%> Kayhan Batmanghelich
%>
%> @b Link: 
%> https://www.cbica.upenn.edu/sbia/software/
%> 
%> @b Contact: 
%> sbia-software@uphs.upenn.edu
%======================================================================
% 
function xp = ProjectionOnUnitBoxSimplex(x,lambda)
    %> first project on the box constraints
    x_org = x ;
    x(x < 0) = 0 ;
    x(x > 1) = 1 ;
    xp = x ;
    
    if (sum(xp) > lambda)
        xp = ProjectionOnConstrainedSimplex(double(x_org),double(lambda)) ;  %> it is important to keep "double" because of precision
    end
end

function xp = ProjectionOnConstrainedSimplex(x0,lambda)
    x = x0 ;
    z = min(1,max(0,x)) ;
    l_max = 2*max(x) ;
    l_min = min(x) ;
    
    xx_old = ones(size(x)) ;
    l = l_min + (l_max - l_min)/2 ;
    xx_min = softThresholdTruncate(x,l_min) ;
    xx_max = softThresholdTruncate(x,l_max) ;
    while(1)
        xx = softThresholdTruncate(x,l) ;
        %if (isequal((xx > 0),(xx_old>0)))     %if index does not change break and try to find the accurate solution
         
        %xx_old = xx ;
        if (sum(xx)>lambda) %> threshold is too low
            l_min = l ;
            l = l_min + (l_max - l_min)/2 ;
            xx_min = xx ;
        elseif (sum(xx)<lambda) %> threshold is too high
            l_max = l ;
            l = l_min + (l_max - l_min)/2 ;
            xx_max = xx ;
        else     %> sum(xx)==lambda, this event happens with very-very low probability :(
            xp = xx ;
            break ;
        end
        if (isequal((xx_min == 0),(xx_max == 0)))
            if (isequal((xx_min == 1),(xx_max == 1)))
                ind = ((xx>0) & (xx < 1)) ;
                l = -(lambda - sum(xx==1) - sum(x(ind)) )/sum(ind)*2 ;
                xp = softThresholdTruncate(x,l) ;
                break ;
            end
        end
        if ((l_max - l_min)/2 < 1e-3 )
%            warning('very slow convergence on projection') ;
            ind = ((xx>0) & (xx < 1)) ;
            l = -(lambda - sum(xx==1) - sum(x(ind)) )/sum(ind)*2 ;
            xp = softThresholdTruncate(x,l) ;
            break ;
        end
        %if (abs(sum(xx) - lambda)<0.0001)
        %    xp = xx ;
        %    break ;
        %end
    end
    
end


function xx = softThresholdTruncate(x,l)
    xx = (x - l/2).*(x>l/2) ;
    xx = min(xx,1) ;
end

