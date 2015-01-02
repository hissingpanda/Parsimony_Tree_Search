function out = sankoffAmino(node_l, node_r)

%{
    Function sankoffAmino that performs Sankoff's algorithm on two phylogenetic
     nodes
         Input: Left and right node of a branch
        Output: The parsimony values for each amino acid at the node 
%}

% Initializes output array as 24 elements of infinity
out = Inf(1, 24);

%{
    Loops through 24 possible amino acids (includes 'B', 'Z', 'X', and '-') 
     at each node keeping track of the minimum parsimony at each position
     by using the scoring matrix and current position.
%}

% Uses pam250 scoring matrix
scor = pam(250);

len = length(node_l);
for i=1:len
    
    left = inf; 
    for j=1:len      
        cur = scor(i,j) + node_l(j);
        left = min(left, cur);
    end
    
    right = inf;
    for j=1:len
        cur = scor(i,j) + node_r(j);
        right = min(right, cur);
    end
    out(i) = left + right; 
    
end
%display(out);
