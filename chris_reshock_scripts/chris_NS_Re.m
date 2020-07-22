function [Ru,Rv] = chris_NS_Re(u,v,rho,mu)

    [~,dru2] = chris_gradient(rho.*u.*u,1,1);
    [drv2,~] = chris_gradient(rho.*v.*v,1,1);
    [druvdy,druvdx] = chris_gradient(rho.*u.*v,1,1);

    d2u = chris_laplacian_2D(u,1,1);
    d2v = chris_laplacian_2D(v,1,1);

    Ru = sqrt(mean(((dru2+druvdy).^2)')./mean(((mu.*d2u).^2)'));
    Rv = sqrt(mean(((drv2+druvdx).^2)')./mean(((mu.*d2v).^2)'));

end

