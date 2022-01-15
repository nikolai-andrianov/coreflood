function fn = coreyCapPressureAD(pd, lambda, Swr, Sor)
% Corey capillary pressure for AD style fluids
    fn = @(sw, varargin) coreyCapPressure(sw, pd, lambda, Swr, Sor);
end

function pc = coreyCapPressure(sw, pd, lambda, Swr, Sor)
    % Ensure sat > 0 to avoid infinite pc
    %eps = 1e-8;
    eps = 0;
    sat = (sw - (Swr - eps))./(1 - Swr - Sor);
    sat = max(min(sat, 1), 0);
    
    pc = pd*sat.^(-1./lambda);   
   
    if isa(sat, 'ADI')
        if any(~isfinite(pc.val)) 
          ind = ~isfinite(pc.val);
          eps = 1e-8;
          pc_max = pd*eps.^(-1./lambda); 
          pc.val(ind) = pc_max;
        end
       
        dpcdsw = pc.jac{2};
        [i, j, jac] = find(dpcdsw);
        ind = ~isfinite(jac);
        if any(ind) 
            if all(sat.val(ind) == 1) || all(sat.val(ind) == 0)
                dpcdsw(i(ind), j(ind)) = 0;                
                pc.jac{2} = dpcdsw;
            else
                a = 1;
            end
        end   
    end    
  
end
