function [tree, minp, tracker] = nni2(tree, parsimony, alignment)

%{
    Function nni that searches for the best possible phylogenetic tree, using
    the nearest-neighbor interchange algorithm
         Input: phylogenetic tree, parsimony score, multiple sequence alignment
        Output: phylogenetic tree, parsimony score, value to see if score changed
%}

%Create a row vector of the leaf names in alphabetical order.
[dummy,order] = sort(get(tree,'LeafNames'));
numBranches = get(tree,'NumBranches');

%Initialize values
matOrder(1,:) = order;
matTrees{1} = tree;
minp = parsimony;
tracker = 1;

%Reorder the phylogenetic tree to match as closely as possible the row vector 
%of alphabetically ordered leaf names, without dividing the clades or having 
%crossing branches.
for i=2:4
    [tree2, order] = reorder(tree,order,'approximate',true);
    display('Finished a reorder');
    matTrees{i} = tree2;
    matOrder(i,:) = order;
    view(tree2);
    display(i);
    display('Starting parsimony calculation');
    parsimony(i) = treeParseAmino(matTrees{i}, alignment);
    display('Finished a calculation');
    display(parsimony);
    
    %check to see if the current tree has the minimum parsimony score
    if (parsimony(i) < minp)
        minp = parsimony(i);
        tracker = i;
    end
end
%display(matOrder);
%display(minp);

%checks to see if the parsimony of was ever lower than the original score
if (tracker > 1)
    tree = matTrees{tracker};
end
