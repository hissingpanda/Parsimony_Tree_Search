function [tree, parsimony] = parsimonyTreeSearch2(alignment)

%{

    Todo:
        recreate 6, 7 seed trees
        get parsimony scores, txt document
        compute ncbi tree
        write essay

    Function parsimonyTreeSearch that performs three tasks:
        1) Builds initial phylogenic tree using pairwise distances of the
            input multiple sequence alignment
        2) Computes the tree's parsimony score
        3) Searches for a better tree (and generates it) using the parsimony score

         Input: multiple sequence alignment
        Output: phylogenetic tree, parsimony score for the input tree

    Called like:
       parsimonyTreeSearch('PF02171/PF02171_seed.fasta')
       tree = parsimonyTreeSearch('PF02171/PF02171_full_dashes.fasta')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1: Build Initial Tree  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Starting Sequence Alignment');
seqs = fastaread(alignment);

% Computes pairwise distances with pam250 matrix as the scoring matrix
%{
distances = seqpdist(seqs,'Method','alignment-score',...
                'Indels','pairwise-delete',...
                'ScoringMatrix','pam250');
%}

% Need to use this one for full sequence (pairwise-delete errors out otherwise)

display('Computing Distances');
distances = seqpdist(seqs);
% Generates and view phylogenetic tree
display('Creating initial tree')
tree = seqneighjoin(distances, 'equivar', seqs);
%view(tree);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 2: Compute Tree's Parsimony Score  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Starting initial parse');
parsimony = treeParseAmino(tree, alignment);
display('Finished initial parse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 3: Implement Tree Search Strategy  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start searching for a better tree, using the nni function
display('Starting search for better tree');
view(tree);
[tree, parsimony, tracker] = nni2(tree, parsimony, alignment);
count = 1;
view(tree);
display(parsimony);

% If better tree was found, recursivelly apply the nni function to
% subsequent trees
%{
while (tracker > 1)
    count = count + 1;
    [tree, parsimony, tracker] = nni(tree, parsimony, alignment);
    %display(tracker);
    view(tree);
    display(parsimony);
    
    % Breaks at standard value to prevent long run-times
    if (count >= 4)
        break;
    end
end
%}
%{
 [tree3, order] = reorder(tree,order,'approximate',true);
    mat(3,:) = order;
    [tree4, order] = reorder(tree,order,'approximate',true);
    mat(4,:) = order;
    [tree5, order] = reorder(tree,order,'approximate',true);
    mat(5,:) = order;
    [tree6, order] = reorder(tree,order,'approximate',true);
    mat(6,:) = order;
    display(mat);
%}
%View the original and the reordered phylogenetic trees.
%view(tree)
%view(tree2)
%view(tree3)
%view(tree4)
%view(tree5)
%view(tree6)

