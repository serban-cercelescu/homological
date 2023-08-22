module Matrix (
    Matrix,
    mbounds,
    mulMatrix,
    identity,
    transpose,
    matToList,
    matToString
) where

import Data.Array ( Array, array, (!), bounds )
import Data.List (intercalate)

type Matrix a = Array (Int, Int) a

mbounds :: Matrix a -> (Int, Int)
mbounds matrix = (n, m)
    where
        (n, m) = snd $ bounds matrix

mulMatrix :: (Num a) => Matrix a -> Matrix a -> Matrix a
mulMatrix a b = if m /= m' then error "Matrix dimensions do not match" else array ((1, 1), (n, k)) [((i, j), sum [a ! (i, l) * b ! (l, j) | l <- [1 .. m]]) | i <- [1 .. n], j <- [1 .. k]]
    where
        (n, m) = mbounds a
        (m', k) = mbounds b

identity :: (Num a) => Int -> Matrix a
identity n = array ((1, 1), (n, n)) [((i, j), if i == j then 1 else 0) | i <- [1 .. n], j <- [1 .. n]]

matToList :: Matrix a -> [[a]]
matToList matrix = [[matrix ! (i, j) | j <- [1 .. m]] | i <- [1 .. n]]
    where
        (n, m) = mbounds matrix

transpose :: Matrix a -> Matrix a
transpose matrix = array ((1, 1), (m, n)) [((j, i), matrix ! (i, j)) | j <- [1 .. m], i <- [1 .. n]]
    where
        (n, m) = mbounds matrix

matToString :: (Show a) => Matrix a -> String
matToString matrix = intercalate " " $ fmap show $ concat $ matToList matrix
