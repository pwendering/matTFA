function isTrans = checkIfRxnIsTransport(reactants,products)
% This function simply check if there is one or more same metabolite(s)
% on both sides of the reaction. If there is at least one common
% metabolite, then the reaction is classified as transport.
% INPUTS
%  - reactants:
%  - products:
% 
% OUTPUTS
% - isTrans: is a flag. 1 if the reaction is classified as transport
%   reaction, and 0 otherwise.
 
isTrans = 0;

reactants=columnVector(reactants);
products=columnVector(products);

react=repmat({''},numel(reactants),1);
prod=repmat({''},numel(products),1);

for i=1:length(reactants)
    react{i}=parseMet(reactants{i});
end
 
for i=1:length(products)
    prod{i}=parseMet(products{i});
end
for i=1:length(reactants)
    [~, ba]=ismember(prod, react(i));
    if length(find(ba))==1
        isTrans=1;
    elseif find(ba)>1
        isTrans=0;
        return
    end
end
 
for i=1:length(products)
    [~, ba]=ismember(react, prod(i));
    if length(find(ba))==1
        isTrans=1;
    elseif find(ba)>1
        isTrans=0;
        return
    end
end
end