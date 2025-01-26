function pos_helical = ConvertCartersiantoHelical(pos_cart,alph,tau)

pos_helical = zeros(size(pos_cart)) ;

pos_helical(:,1) = sqrt(pos_cart(:,1).^2+pos_cart(:,2).^2) ;
pos_helical(:,2) = pos_cart(:,3)/tau ;
pos_helical(:,3) = (1/(2*pi))*atan2(pos_cart(:,2),pos_cart(:,1)) - (alph*pos_cart(:,3))/tau ;

end
