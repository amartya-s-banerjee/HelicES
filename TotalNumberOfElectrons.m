function f = TotalNumberOfElectrons(lambda_f_g, EigVal, Nelectron, bet, nkpt, wkpt)

f = 0 ;
for kpt=1:nkpt
   f = f + wkpt(kpt)*2*sum(1./(1+exp(bet*(EigVal(:,kpt)-lambda_f_g)))) ;
end
f = f - Nelectron;

% f = (2*sum(sum(1./(1+exp(bet*(EigVal-lambda_f_g))))))/nkpt-Nelectron; 

end