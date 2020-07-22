function [kx,Ex] = chris_exp_eng_spec_hs(C,eta,dpxdx)

[ny,nx,nt] = size(C);

for i = 1:nt
     Ci = squeeze(C(:,:,i));
     
    Cm = mean(Ci');
   
    C_p = Ci;
    
    for j = 1:nx
        C_p(:,j) = Ci(:,j) - Cm';
    end
    
    itop = find(eta(:)>-0.5);%,i
    range_cp = find(eta(itop)<0.5);%,i
    if ~isempty(range_cp)
        n=0;
        for j = range_cp(1):range_cp(end)
            n = n+1;
            [Eax(n,:),kax] = pspectrum(C_p(j,:),[1:nx]/dpxdx);
        end
        Ex(:,i) = mean(Eax);
        kx(:,i) = kax;
    end
    
end

end

