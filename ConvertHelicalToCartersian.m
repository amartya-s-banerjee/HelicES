function pos_cart = ConvertHelicalToCartersian(pos_helical,alph,tau)

pos_cart = zeros(size(pos_helical)) ;

pos_cart(:,1) = pos_helical(:,1).*cos(2*pi*(alph*pos_helical(:,2)+pos_helical(:,3))) ;
pos_cart(:,2) = pos_helical(:,1).*sin(2*pi*(alph*pos_helical(:,2)+pos_helical(:,3))) ;
pos_cart(:,3) = pos_helical(:,2)*tau;

end
