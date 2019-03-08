module ProgressiveAlignment where

import Test.HUnit
import Data.Array
import qualified Data.Set as Set
import qualified Data.Map as Map
import UPGMA
import Distances
import Needleman

{- treeTraversal guidetree
   Creates multiple sequence alignment (MSA) by traversing the guide tree and applying NWA.
   RETURNS: MSA of sequences in leaves of guidetree. 
   VARIANT: number of nodes not visited.
 -}
treeTraversal :: UPGMATree -> [String]
treeTraversal (Node l str v r) = traceBack (traceBack (treeTraversal l) (treeTraversal r)) str
treeTraversal (Leaf seq v) = [seq]

{- newickFormat tree
   Create Newick format of UPGMA tree.
   PRE: tree is non-empty. 
   RETURN: tree in Newick format.
   VARIANT: tree will eventually be just a leaf.
 -}
newickFormat :: UPGMATree -> String
newickFormat (Node l str v r) = "(" ++ (newickFormat l) ++ ":" ++ (show (edgeLength (Node l str v r) l)) ++ "," ++ (newickFormat r)++ ":" ++ (show (edgeLength (Node l str v r) r)) ++")"
newickFormat (Leaf seq v) = seq

sequenceToLabel :: Map.Map [Char] [Char] -> String -> Maybe String
sequenceToLabel dictionary seq = Map.lookup seq dictionary

{- edgeLength node1 node2
   Calculate length of edge between nodes.
   RETURNS: length of edge between node1 and node2.
 -}
edgeLength :: UPGMATree -> UPGMATree -> Double
edgeLength nodeA nodeB = (extractHeight nodeA) - (extractHeight nodeB)