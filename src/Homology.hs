module Homology where

import SNF (smithNormalForm)
import Matrix ( Matrix, mulMatrix, mBounds )
import Data.Array ( (!) )
import Text.Printf ( printf )


type ZModule = [Integer]
type ZChain = [Integer]

{-
    ℤ^n -A-> ℤ^m -B-> ℤ^k is part of a chain complex (i.e. BA = 0)
    the homology module is the direct sum of the torsion and free parts of Ker[B]/Im[A]
    this is represented as a list of integers, where there the positive integers are the orders of the torsion parts
    and the zero integers each represent a free part, so that for example [2, 3, 0, 0] represents
    ℤ/2ℤ ⊕ ℤ/3ℤ ⊕ ℤ ⊕ ℤ
-}
homologyModule :: Matrix Integer -> Matrix Integer -> ZModule
homologyModule a b = ans where
    (na, ma) = mBounds a
    (nb ,mb) = mBounds b

    (_, smithA, _) = smithNormalForm a
    (_, smithB, _) = smithNormalForm b
    torsion = filter (>= 1) [smithA ! i ! i | i <- [1 .. min na ma]]
    bNullity = mb - length (filter (/= 0) [smithB ! i ! i | i <- [1 .. min nb mb]])
    free = replicate (bNullity - length torsion) 0

    ans
      | na == 0 || ma == 0                                                       = replicate bNullity 0
      | nb == 0 || mb == 0                                                       = filter (/= 1) $ torsion ++ replicate (na - length torsion) 0
      | mb /= na                                                                 = error $ printf "Incompatible matrix dimensions of A (%d, %d) and B (%d, %d) to compute the product BA \n A : %s \n B : %s" na ma nb mb (show a) (show b)
      | not $ all (==0) [mulMatrix b a ! i ! j | i <- [1 .. nb], j <- [1 .. ma]] = error $ printf "Matrix product is not zero"
      | otherwise                                                                = filter (/= 1) $ torsion ++ free

homologyTorsion :: Matrix Integer -> Matrix Integer -> [Integer]
homologyTorsion a b = filter (/= 0) $ homologyModule a b

homologyFreeRank :: Matrix Integer -> Matrix Integer -> Int
homologyFreeRank a b = length $ filter (== 0) $ homologyModule a b
