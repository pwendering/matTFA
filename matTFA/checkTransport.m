function model = checkTransport(model)
% checks if the reaction is a transport reaction
% a transport reaction is defined as having the same metabolite(s) in different
% compartments
model.isTrans = zeros(numel(model.rxns),1);
for i=1:length(model.rxns)
   
   reactants = model.mets(find(model.S(:,i) < 0));
   products = model.mets(find(model.S(:,i) > 0));
   
   model.isTrans(i) = checkIfRxnIsTransport(reactants,products);
end

end
