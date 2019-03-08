module Distances where

import Test.HUnit
import Data.Array
import qualified Data.Set as Set


{- naiveDistance a b
   Calculates proportion of sites at which two aligned sequences differ.
   PRE: a and b are aligned (and therefore of the same length).
   RETURNS: the naive distance between a and b.
   EXAMPLES: naiveDistance "tes-" "test" = 0.25
             naiveDistance "test" "test" = 0.0
-}
naiveDistance :: String -> String -> Double
naiveDistance x y = mismatches/aligned where
    mismatches = fromIntegral (misMatches x y 0)
    aligned    = fromIntegral (length x)

{- misMatches a b diff
   Find number of differing sites in two strings. 
   RETURNS: the number of mismatches between a and b, returned in diff.
   EXAMPLES: misMatches "abc" "def" 0 = 3
             misMatches "abc-" "abcd" 0 = 1
   VARIANT: (length a) + (length b)
-}
misMatches :: String -> String -> Int -> Int
misMatches [] [] diff = diff
misMatches (x:xs) (y:ys) diff
    | x == y = misMatches xs ys diff
    | otherwise = misMatches xs ys (diff+1)


{- jukesCantorDistance a b
   Improved distance measure, accounts for hidden mutations.
   PRE: a and b are aligned (and therefore of the same length).
   RETURNS: Jukes-Cantor estimated distance for a and b.
 -}
jukesCantorDistance :: String -> String -> Double
jukesCantorDistance a b = (-0.75)*(log (1-((4/3))*(naiveDistance a b)))


{- msaDistance msa
   Calculate distance matrix for MSA.
   PRE: msa contains non-empty strings.
   RETURNS: a matrix of all naive distances between sequences in the MSA.
-}
msaDistance :: [String] -> Array (Int, Int) Double
msaDistance msa = matrix where
    m = (length msa) - 1
    matrix = array ((0,0),(m,m)) [((x,y), naiveDistance (msa!!x) (msa!!y)) | x<-[0..m], y<-[0..m]]


{- kmerDistance sequenceA sequenceB
Calculates distance as the fraction of words that are unique to either sequence.
PRE: Neither sequence is shorter than k (default k=3).
RETURNS: k-mer distance between sequenceA and sequenceB. 
EXAMPLES: kmerDistance "GATTACA" "TAGAT" = 6/7
-}
kmerDistance :: Fractional a => String -> String -> a
kmerDistance a b = unique/total where
    x = wordComposition a 3 True
    y = wordComposition b 3 True

    uniqueX = (map f x) where
        f k = (elem k y)

    uniqueY = (map f y) where
        f k = (elem k x)
    
    unique = fromIntegral (length (filter (==False) (uniqueX++uniqueY)))
    total = fromIntegral (length (Set.fromList (x++y)))
        

{- wordComposition sequence k overlap 
   Find all sequential k-mer words of a sequence, option for overlapping or not. 
   PRE: k > 0.
   RETURNS: all words of length k from sequence, overlapping if overlap is True. 
   EXAMPLES: wordComposition "EXAMPLE" 3 False = ["EXA","MPL"]
             wordComposition "EXAMPLE" 3 True = ["EXA","XAM","AMP","MPL","PLE"]
   VARIANT: (length seq) - k + 1
 -}
wordComposition :: String -> Int -> Bool -> [String] 
wordComposition seq k overlap
    | length seq < k = []
    | overlap == True = (take k seq):(wordComposition (drop 1 seq) k overlap)
    | otherwise = (take k seq):(wordComposition (drop k seq) k overlap)
