import Test.HUnit
import Data.Array

import Needleman 
import ProgressiveAlignment
import Distances
import UPGMA

-- Test Cases

-- UPGMA.hs
testU1 = TestCase $ assertEqual "wordComposition TEST overlapping" ["TE", "ES", "ST"] (wordComposition "TEST" 2 True)
testU2 = TestCase $ assertEqual "wordComposition TEST without overlapping" ["TE", "ST"] (wordComposition "TEST" 2 False)
testU3 = TestCase $ assertEqual "wordComposition [] overlapping" [] (wordComposition "" 2 True)
testU4 = TestCase $ assertEqual "wordComposition [] without overlapping" [] (wordComposition "" 2 False)
testU5 = TestCase $ assertEqual "wordComposition with k > sequence, overlapping" [] (wordComposition "TEST" 5 True)
testU6 = TestCase $ assertEqual "wordComposition with k > sequence, without overlapping" [] (wordComposition "TEST" 5 False)
testU7 = TestCase $ assertEqual "wordComposition with k == sequence, overlapping" ["TEST"] (wordComposition "TEST" 4 True)
testU8 = TestCase $ assertEqual "wordComposition with k == sequence, without overlapping" ["TEST"] (wordComposition "TEST" 4 False)
testU10 = TestCase $ assertEqual "kmerDistance with seq = seq" 0.0 (kmerDistance "ABC" "ABC")
testU11 = TestCase $ assertEqual "kmerDistance seq1 < seq2" 0.5 (kmerDistance "ABC" "ABCD")
testU12 = TestCase $ assertEqual "kmerDistance seq1 > seq2" 0.5 (kmerDistance "ABCD" "ABC")
testU14 = TestCase $ assertEqual "clusterDistance [tes,tis] [tas,tos] False" 1.0 (clusterDistance ["tes","tis"] ["tas","tos"] False)
testU15 = TestCase $ assertEqual "minDistance (array ((0,0),(1,1)) [((0,0),0.0),((0,1),0.888),((1,0),0.888),((1,1),0.0)])" ((0,1),0.888) (minDistance (array ((0,0),(1,1)) [((0,0),0.0),((0,1),0.888),((1,0),0.888),((1,1),0.0)]))

-- Needleman.hs
testN1 = TestCase $ assertEqual "sumPairScore [TEST] 0 [TEST] 0 = 1" 1.0 (sumPairScore ["TEST"] 0 ["TEST"] 0)
testN2 = TestCase $ assertEqual "sumPairScore [-EST] 0 [-EST] 0 = 0" 0.0 (sumPairScore ["-EST"] 0 ["-EST"] 0)
testN3 = TestCase $ assertEqual "sumPairScore [TEST,TEST] 0 [TEST] 0 = 2" 1.0 (sumPairScore ["TEST","TEST"] 0 ["TEST"] 0)
testN4 = TestCase $ assertEqual "spAux test 0 test 0 = 1" 1 (spAux "test" 0 "test" 0)
testN5 = TestCase $ assertEqual "traceBack [test] [test]" ["test","test"] (traceBack ["test"] ["test"])
testN6 = TestCase $ assertEqual "traceBack [testing] [test]" ["testing","test---"] (traceBack ["testing"] ["test"])
testN7 = TestCase $ assertEqual "traceBack [IAMAPEPTIDE] [IAMPEPTIDE]" ["IAMAPEPTIDE","IAM-PEPTIDE"] (traceBack ["IAMAPEPTIDE"] ["IAMPEPTIDE"])
testN8 = TestCase $ assertEqual "traceBack [IAMAPEPTIDE, IAM-PEPTIDE] [IAMPEPPED]" ["IAMAPEPTIDE","IAM-PEPTIDE","IAM-PEPPED-"] (traceBack ["IAMAPEPTIDE", "IAM-PEPTIDE"] ["IAMPEPPED"])
testN9 = TestCase $ assertEqual "traceBack [IAMPEPPED] [IAMAPEPTIDE, IAM-PEPTIDE] " ["IAM-PEPPED-","IAMAPEPTIDE","IAM-PEPTIDE"] (traceBack ["IAMPEPPED"] ["IAMAPEPTIDE", "IAM-PEPTIDE"])

-- Distances.hs
testD1 = TestCase $ assertEqual "naiveDistance test test" 0.0 (naiveDistance "test" "test")
testD2 = TestCase $ assertEqual "naiveDistance test tes-" 0.25 (naiveDistance "test" "tes-")
testD3 = TestCase $ assertEqual "misMatches test test 0" 0 (misMatches "test" "test" 0)
testD4 = TestCase $ assertEqual "misMatches test tes- 0" 1 (misMatches "test" "tes-" 0)
testD5 = TestCase $ assertEqual "msaDistance test test" (array ((0,0),(1,1)) [((0,0),0.0),((0,1),0.0),((1,0),0.0),((1,1),0.0)]) (msaDistance ["test","test"])
testD6 = TestCase $ assertEqual "jukesCantorDistance ATGCTCGTAGCTGCTACGTC ATCCTCGAAGCAGCTACGAC" 0.2326161962278796 (jukesCantorDistance "ATGCTCGTAGCTGCTACGTC" "ATCCTCGAAGCAGCTACGAC")

-- ProgressiveAlignment.hs
testPA1 = TestCase $ assertEqual "treeTraversal upgmaClustering" ["IAMAPEPTIDE","IAM-PEPTIDE","IAM-PEPPED-"] (treeTraversal ((upgmaClustering ["IAMAPEPTIDE", "IAMPEPTIDE","IAMPEPPED"] False)!!0))

runtests = runTestTT (TestList [testU1, testU2, testU3, testU4, testU5, testU6, testU7, testU8, testU10, testU11, testU12, testU14, testU15, 
                                testN1, testN2, testN3, testN4, testN5, testN6, testN7, testN8, testN9, testPA1, testD1, testD2, testD3, testD4, testD5, testD6])
