module UPGMA where

import Test.HUnit
import Data.Array
import Distances


{- data UPGMATree
   UPGMATree represents an UPGMA tree.
   It is a full binary tree, leafs contain a sequence and 0 (height).
   The nodes contain a string ("") and height from leaf.
   The tree is ultrametric - height from node to each leaf below it is the same.
 -}
data UPGMATree = Leaf String Double | Node (UPGMATree) [String] Double (UPGMATree) deriving Show


{- upgmaClustering seqs
   Performs UPGMA clustering on sequences (recursively merging of clusters in cluster list
   and recalculation of the distance matrix from these clusters).
   RETURNS: UPGMATree with each sequence present in seqs stored in leafs.
-}
upgmaClustering :: [String] -> Bool -> UPGMATree
upgmaClustering sequences final = upgmaClusteringAux clusters trees final where
    clusters = initialClusters sequences
    trees = initialTrees sequences


{- upgmaClusteringAux clusters trees
   Update clusters and trees based on new distance matrix.
   RETURNS: List of all new UPGMA trees.
   VARIANT: (length (tail clusters))
 -}
upgmaClusteringAux :: [[String]] -> [UPGMATree] -> Bool -> UPGMATree
upgmaClusteringAux (cluster:[]) trees final = trees!!0
upgmaClusteringAux clusters trees final = upgmaClusteringAux new_clusters new_trees final where
    distance_matrix = clusterDistanceMatrix clusters final -- Calculate distance matrix.
    index = fst (minDistance distance_matrix)              -- Find index of closest clusters.
    distance = snd (minDistance distance_matrix)           -- Get distance between closest clusters.
    new_clusters = updateClusters clusters index           -- Replace clusters with merged.
    new_trees = updateTrees trees index distance           -- Replace trees with merged.

{- initialClusters seqs
   Create initital list of clusters from sequences.
   RETURNS: list of clusters consisting of single sequences from seqs.
 -}
initialClusters :: [String] -> [[String]]
initialClusters sequences = [[x] | x <- sequences]


{- initialTrees seqs
   Create initital list of trees from sequences.
   RETURNS: list of Leafs containing each sequence from seqs and 0.0.
   EXAMPLES:
 -}
initialTrees :: [String] -> [UPGMATree]
initialTrees sequences = [Leaf x 0.0 | x <- sequences]


{- clusterDistanceMatrix clusters
   Create distance matrix from clusters.
   PRE: clusters is non-empty and contain no empty string.
   RETURNS: a distance matrix of the clusters.
 -}
clusterDistanceMatrix :: [[String]] -> Bool -> Array (Int, Int) Double
clusterDistanceMatrix clusters final = matrix where
    m = (length clusters) - 1
    matrix = array ((0,0),(m,m)) [((x,y), (clusterDistance (clusters!!x) (clusters!!y)) final) | x<-[0..m], y<-[0..m]]


{- clusterDistance clusterA clusterB final
   Calculates the mean distance between elements of each cluster.
   Uses Jukes as distance measure instead of kmer if final is true.
   PRE: clusterA or clusterB is non-empty.
   RETURNS: mean k-mer distance between clusterA and clusterB.
   EXAMPLES: clusterDistance ["tes"] ["tis"] = 1.0
             clusterDistance ["tes","tis"] ["tas"] = (1.0 + 1.0)/2 = 1.0
             clusterDistance ["tes","tis"] ["tas","tos"] = (1.0 + 1.0 + 1.0 + 1.0)/4 = 1.0
 -}
clusterDistance :: [String] -> [String] -> Bool -> Double
clusterDistance clusterA clusterB False = pair_distances/num_distances where
    pair_distances = sum [kmerDistance a b | a <- clusterA, b <- clusterB]
    num_distances = fromIntegral ((length clusterA)*(length clusterB))

-- Use Jukes-Cantor distance if final is True (building final tree from msa).
clusterDistance clusterA clusterB True = pair_distances/num_distances where
    pair_distances = sum [jukesCantorDistance a b | a <- clusterA, b <- clusterB]
    num_distances = fromIntegral ((length clusterA)*(length clusterB))


{- minDistance distancematrix
   Find closest sequences in distance matrix.
   PRE: distance matrix is non-empty.
   RETURNS: index of smallest value in distancematrix.
 -}
minDistance :: Array (Int, Int) Double -> ((Int, Int), Double)
minDistance m = foldl1 cmp (assocs m)


{- cmp coord_and_val_1 coord_and_val_2
   Determine which of two places of matrix contain smaller value.
   RETURNS: the coordinates of coord_and_val_1 or coord_and_val_2 containing the smallest value.
-}
cmp :: ((Int, Int), Double) -> ((Int, Int), Double) -> ((Int, Int), Double)
cmp (c1, v1) (c2, v2)
    -- check if one coord-pair is in the main diagonal, if so, return the other pair
    | x1 == y1 = (c2, v2)
    | x2 == y2 = (c1, v1)
    -- otherwise, compare the two values and return the coords for the smallest one
    | v1 <= v2 = (c1, v1)
    | otherwise = (c2, v2) where
        (x1, y1) = c1
        (x2, y2) = c2


{- updateClusters clusters min_coordinates
   Removes closest clusters and adds them as a single cluster.
   RETURNS: clusters with with the two clusters at min_coordinates merged.
 -}
updateClusters :: [[String]] -> (Int, Int) -> [[String]]
updateClusters clusters (i,j) = merged_clusters ++ filtered_clusters where
    merged_clusters = [clusters!!i ++ clusters!!j]
    filtered_clusters = [clusters!!x| x <- [0..((length clusters) -1)], x /= i && x /=j]


{- updateTrees trees (i,j) h
   Merge trees corresponding to closest clusters and remove the previously unmerged trees.
   PRE: (i,j) contains two valid indices of trees.
   RETURNS: updated trees with trees previously at position i and j merged and
            added to end of trees with height calculated from h.
 -}
updateTrees :: [UPGMATree] -> (Int, Int) -> Double -> [UPGMATree]
updateTrees trees (i,j) distance = merged_trees ++ filtered_trees where
    height = (distance + (extractHeight (trees!!i)) + (extractHeight (trees!!j)))/2 --arclength ground->ground / 2
    merged_trees = [Node (trees!!i) [""] height (trees!!j)]
    filtered_trees = [trees!!x| x <- [0..((length trees) -1)], x /= i && x /=j]


{- extractHeight tree
   Get height stored in root of tree.
   RETURNS: height of tree.
 -}
extractHeight :: UPGMATree -> Double
extractHeight (Node _ _ height _) = height
extractHeight (Leaf _ height) = height