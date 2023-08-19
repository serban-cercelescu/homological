module Simplicial (
  Simplex(..),
  SimplicialComplex(..),
  newSimplicialComplex,
  nSimplices,
  nthBoundary,
  nthHomology,
  nthBettiNumber,
  nthTorsionNumbers,
  homologySequence,
  eight,
  rp2,
  torus
) where

import qualified Data.Map.Strict as Map
import Data.List ( intercalate, sort, sortBy)
import Data.Array ( array )
import Matrix ( Matrix, transpose )
import Homology ( ZModule, homologyModule )


newtype Simplex = Simplex [Int]
data SimplicialComplex = SimplicialComplex {
  numVertices :: Int,
  maximalSimplices :: [Simplex]
}

instance Show Simplex where
  show (Simplex xs) = "(" ++ intercalate "," (show <$> xs) ++ ")"

instance Show SimplicialComplex where
  show (SimplicialComplex _ xs) = show xs

instance Eq Simplex where
  (Simplex xs) == (Simplex ys) = sort xs == sort ys

instance Ord Simplex where
  (<=) (Simplex xs) (Simplex ys) = sort xs <= sort ys

instance Eq SimplicialComplex where
  (==) (SimplicialComplex nx xs) (SimplicialComplex ny ys) = sort xs == sort ys && nx == ny

{-
  the arguments are (n) the number of vertices and (xs) the list of simplices
  it simply stores the maximal simplices and puts them into a normal form
  (each simplex has its vertices sorted in ascending order, and the list of maximal simplices is sorted lexicographically)
-}
newSimplicialComplex :: Int -> [[Int]] -> SimplicialComplex
newSimplicialComplex n xs0 = SimplicialComplex n xs where
  xs1 = sortBy (\xs ys -> (length xs, xs) `compare` (length ys, ys)) $ sort <$> xs0
  xs = Simplex <$> antichain xs1

  antichain [] = []
  antichain (x:xs) = if any (isContained x) xs then antichain xs else x : antichain xs

nSimplices :: SimplicialComplex -> Int -> [Simplex]
-- nSimplices (-1) _ = [Simplex []] -- reduced homology
nSimplices (SimplicialComplex numVertices xs) n' = ans where
  n = n' + 1
  nSpx _ _ [] = []
  nSpx spx lspx maximals = if lspx == n then [spx] else do
    let avb = if null spx then numVertices else head spx - 1
    let needed = n - lspx
    newHead <- [needed .. avb]
    let newSimplex = newHead:spx
    nSpx newSimplex (lspx + 1) (filter (isContained newSimplex) maximals)

  ans
    | n == 0    = []
    | otherwise = Simplex <$> nSpx [] 0 ((\(Simplex x) -> x) <$> xs)  

nthBoundary :: SimplicialComplex -> Int -> Matrix Integer
nthBoundary cpx n = boundaryMap where
  simplices0 = nSimplices cpx n
  faceSimplices = nSimplices cpx (n - 1)

  s0 = length simplices0
  sf = length faceSimplices

  faceToInt = Map.fromList $ zip faceSimplices [1 .. sf]

  (%) = (Map.!)

  faces :: Simplex -> [Simplex]
  faces (Simplex xs) = Simplex <$> [take i xs ++ drop (i + 1) xs | i <- [0 .. length xs - 1]]

  boundaryMap = transpose $ array (1, s0) $ do
    (i, spx) <- zip [1 .. s0] simplices0
    let iFaces = Map.fromList $ (\(face, idx) -> (face, if odd idx then -1 else 1)) <$> zip ((faceToInt %) <$> faces spx) [0..]
    let matRow = array (1, sf) [(j, Map.findWithDefault 0 j iFaces) | j <- [1 .. sf]]
    return (i, matRow)


nthHomology :: SimplicialComplex -> Int -> ZModule
nthHomology spx n = homologyModule (nthBoundary spx (n + 1)) (nthBoundary spx n)

nthBettiNumber :: SimplicialComplex -> Int -> Int
nthBettiNumber cpx n = length $ filter (==0) $ nthHomology cpx n

nthTorsionNumbers :: SimplicialComplex -> Int -> [Integer]
nthTorsionNumbers cpx n = filter (/=0) $ nthHomology cpx n

homologySequence :: SimplicialComplex -> [[Integer]]
homologySequence cpx = ans where
  cpx' = let SimplicialComplex _ cpx' = cpx in cpx'
  dimension = maximum [let Simplex s = s' in length s | s' <- cpx']

  boundaryMaps = [nthBoundary cpx i | i <- [0 .. dimension]]
  ans = uncurry (flip homologyModule) <$> zip boundaryMaps (tail boundaryMaps)



----------------------------- Utils ------------------------------

-- given two *ORDERED* lists, check if the first is contained in the second
isContained :: Ord a => [a] -> [a] -> Bool
isContained [] _ = True
isContained (_:_) [] = False
isContained (x:xs) (y:ys)
  | x == y = isContained xs ys
  | x > y = isContained (x:xs) ys
  | otherwise = False


--------------------------- Examples -----------------------------

eight :: SimplicialComplex
eight = newSimplicialComplex 5 [[1, 2], [2, 3], [3, 1], [1, 4], [4, 5], [5, 1]]

rp2 :: SimplicialComplex
rp2 = newSimplicialComplex 6 [
  [1, 3, 4],
  [1, 4, 5],
  [1, 2, 5],
  [2, 3, 4],
  [2, 4, 6],
  [4, 5, 6],
  [3, 5, 6],
  [3, 2, 5],
  [1, 2, 6],
  [1, 3, 6]
  ]

torus :: SimplicialComplex
torus = newSimplicialComplex 7 [
  [1, 3, 7],
  [1, 2, 3],
  [2, 3, 6],
  [3, 4, 6],
  [4, 6, 3],
  [4, 6, 7],
  [1, 6, 7],
  [1, 5, 6],
  [2, 5, 6],
  [2, 5, 7],
  [3, 5, 7],
  [3, 5, 4],
  [1, 4, 5],
  [1, 2, 4],
  [2, 4, 7]
  ]
