function parsimony = treeParseAmino(tree, alignment)

%{
    Function treeParseAmino that performs Sankoff's algorithm on a
    phylogenetic tree.
         Input: phylogenetic tree, multiple sequence alignment
        Output: parsimony score for the input tree
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializes individual amino acids 
%(should put all of these in an array/for loop)-- very gross right now
A = Inf(1, 24); A(1) = 0;
R = Inf(1, 24); R(2) = 0;
N = Inf(1, 24); N(3) = 0;
D = Inf(1, 24); D(4) = 0;
C = Inf(1, 24); C(5) = 0;
E = Inf(1, 24); E(6) = 0;
Q = Inf(1, 24); Q(7) = 0;
G = Inf(1, 24); G(8) = 0;
H = Inf(1, 24); H(9) = 0;
I = Inf(1, 24); I(10) = 0;
L = Inf(1, 24); L(11) = 0;
K = Inf(1, 24); K(12) = 0;
M = Inf(1, 24); M(13) = 0;
F = Inf(1, 24); F(14) = 0;
P = Inf(1, 24); P(15) = 0;
S = Inf(1, 24); S(16) = 0;
T = Inf(1, 24); T(17) = 0;
W = Inf(1, 24); W(18) = 0;
Y = Inf(1, 24); Y(19) = 0;
V = Inf(1, 24); V(20) = 0;
B = Inf(1, 24); B(21) = 0; % Asparagine or aspartic acid
Z = Inf(1, 24); Z(22) = 0; % Glutamine or glutamic acid
X = Inf(1, 24); X(23) = 0; % Unspecified or unknown amino acid
dash = Inf(1, 24); dash(24) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parsimony = 0;

names    = get(tree,'LeafNames');
pointers = get(tree,'Pointers');
numLeaves   = get(tree,'NumLeaves');
numBranches   = get(tree,'NumBranches');
branchNames = get (tree, 'BranchNames');

len = numBranches;
seqs = fastaread(alignment);

% Gets the sequences at each leaf
for i=1:numLeaves
    aligned_seqs{i} = seqs(i).Sequence;
end

%{
    Loops through entire sequence length, translating the lowest parsimony
    value to an amino acid to be used as each internal nodes sequence
    alignment
%}
seqlen = length(aligned_seqs{1});
for z=1:seqlen

    for i=1:len
        % Pointers are the branch/leaf connectivity array
        left = pointers(i,1);
        right = pointers(i,2);
        
        % At a leaf node
        if (left <= numLeaves)
            if (right <= numLeaves)
                leafSeq(1) = aligned_seqs{left}(z);
                leafSeq(2) = aligned_seqs{right}(z);
                
                %bool, true if at leaf, false if at internal node
                atLeaf(1) = 1;
                atLeaf(2) = 1;
        % Second node is at ancestor node
            else
                leafSeq(1) = aligned_seqs{left}(z);
                internalSeq{2} = scores((right - numLeaves),:);
                atLeaf(1) = 1;
                atLeaf(2) = 0;
            end
        % First node is ancestor node
        elseif (left > numLeaves)
            if (right <= numLeaves)
                internalSeq{1} = scores((left - numLeaves),:);
                leafSeq(2) = aligned_seqs{right}(z);
                atLeaf(1) = 0;
                atLeaf(2) = 1;
            % Both nodes are ancestors
            else
                internalSeq{1} = scores((left - numLeaves),:);
                internalSeq{2} = scores((right - numLeaves),:);
                atLeaf(1) = 0;
                atLeaf(2) = 0;
            end
        end
        
        for x=1:2
            %converts to ascii value using giant switch case (should make less sloppy),
            % 65 = 'A' for leaf sequences, this is to give the default parsimony score for
            %each amino acid A = [0, inf, inf, inf, inf, .....]
            
            if (atLeaf(x))
                switch uint8(leafSeq(x))
                    % A, U, G, C, - in order
                    case {uint8('A'), uint8('a')}
                            node{x} = A;
                    case {uint8('R'), uint8('r')}
                            node{x} = R;
                    case {uint8('N'), uint8('n')}
                            node{x} = N;
                    case {uint8('D'), uint8('d')}
                            node{x} = D;
                    case {uint8('C'), uint8('c')}
                            node{x} = C;
                    case {uint8('E'), uint8('e')}
                            node{x} = E;
                    case {uint8('Q'), uint8('q')}
                            node{x} = Q;
                    case {uint8('G'), uint8('g')}
                            node{x} = G;
                    case {uint8('H'), uint8('h')}
                            node{x} = H;
                    case {uint8('I'), uint8('i')}
                            node{x} = I;
                    case {uint8('L'), uint8('l')}
                            node{x} = L;
                    case {uint8('K'), uint8('k')}
                            node{x} = K;
                    case {uint8('M'), uint8('m')}
                            node{x} = M;
                    case {uint8('F'), uint8('f')}
                            node{x} = F;
                    case {uint8('P'), uint8('p')}
                            node{x} = P;
                    case {uint8('S'), uint8('s')}
                            node{x} = S;
                    case {uint8('T'), uint8('t')}
                            node{x} = T;
                    case {uint8('W'), uint8('w')}
                            node{x} = W;
                    case {uint8('Y'), uint8('y')}
                            node{x} = Y;
                    case {uint8('V'), uint8('v')}
                            node{x} = V;
                    case {uint8('B'), uint8('b')}
                            node{x} = B;
                    case {uint8('Z'), uint8('z')}
                            node{x} = Z;
                    case {uint8('X'), uint8('x')}
                            node{x} = X;
                        
                    otherwise
                            node{x} = dash;
                end            
            end
        end
        
        % Performs Sankoff's algorithm on each node by determining if it is
        % an internal node or leaf node
        if (atLeaf(1) && atLeaf(2))
            scores(i,:) = sankoffAmino(node{1}, node{2});
        elseif (atLeaf(1) && atLeaf(2) == 0)
            scores(i,:) = sankoffAmino(node{1}, internalSeq{2});
        elseif (atLeaf(1) == 0 && atLeaf(2))
            scores(i,:) = sankoffAmino(internalSeq{1}, node{2});
        else
            scores(i,:) = sankoffAmino(internalSeq{1}, internalSeq{2});
        end
    end
    
    % Gets the minimum parsimony for each branch at the sequence position
    minimum_val = min(scores,[],2);
    % numLeaves = the root node, what we will use for the parsimony score
    tmp_parsimony = minimum_val(numBranches);
    % adds up parsimony scores at each amino acid in the sequence
    parsimony = parsimony + tmp_parsimony;
    %display(scores);
    
%{    
    % Switch case used to generate a nucleotide sequence at each position
    % for the internal nodes
    for i=1:len
        for k=1:5
            if (scores(i,k) == minimum_val(i))
                switch k
                    % A, U, G, C, -
                    case 1
                        seqs(i,z) = 'A';
                    case 2
                        seqs(i,z) = 'U';
                    case 3
                        seqs(i,z) = 'G';
                    case 4
                        seqs(i,z) = 'C';
                    otherwise
                        seqs(i,z) = '-';
                end            
            end
        end
    end
    %}
    %display(seqs);
    %display(branchNames);
end
%display(seqs);
%display(minimum_val);
