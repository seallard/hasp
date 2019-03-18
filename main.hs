import Needleman 
import ProgressiveAlignment
import Distances
import UPGMA
import qualified Data.Map as Map

{- treeBuilder sequences
   Create phylogenetic tree in Newick format. 
   PRE: sequences contains at least one sequence longer than 3 symbols.
   RETURN: Newick format of phylogenetic tree from sequences.
 -}
treeBuilder :: [String] -> String
treeBuilder seqs = newickFormat (phylo_tree!!0) where
    
    -- Create guide tree by clustering sequences based on k-mer distance matrix.
    guide_tree = upgmaClustering seqs False
    
    -- Traverse guide tree and apply NWA to generate MSA. 
    msa = treeTraversal (guide_tree!!0)

    -- Create phylogenetic tree by clustering sequences in MSA with Jukes-Cantor distance measure. 
    phylo_tree = upgmaClustering msa True