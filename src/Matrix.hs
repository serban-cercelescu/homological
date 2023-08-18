module Matrix (
    Matrix,
    ElemOp(..),
    mBounds,
    mulMatrix,
    rowMul,
    rowAdd,
    rowSwap,
    colMul,
    colAdd,
    colSwap,
    identity,
    transpose,
    matToList
) where

import Data.Array ( Array, array, (!), (//), bounds )
import Control.Monad.State.Lazy ()

type Matrix a = Array Int (Array Int a)

data ElemOp =
  Swap Int Int |
  Add Int Int Integer |
  Mul Int Integer
  deriving (Show, Eq)

mBounds :: Matrix a -> (Int, Int)
mBounds matrix = (n, m) where
    n = snd $ bounds matrix
    m = if n == 0 then 0 else snd $ bounds $ matrix ! 1

transpose :: Matrix a -> Matrix a
transpose matrix = array (1, m) [(i, array (1, n) [(j, matrix ! j ! i) | j <- [1 .. n]]) | i <- [1 .. m]]
  where
    (n, m) = mBounds matrix

mulMatrix :: (Num a) => Matrix a -> Matrix a -> Matrix a
mulMatrix a b = array (1, n) [(i, array (1, k) [(j, sum [a ! i ! l * b ! l ! j | l <- [1 .. m]]) | j <- [1 .. k]]) | i <- [1 .. n]]
    where
        (_, n) = bounds a
        (_, m) = bounds (a ! 1)
        (_, k) = bounds (b ! 1)

rowMul :: (Num a) => Int -> a -> Matrix a -> Matrix a
rowMul i a matrix = matrix // [(i, (matrix ! i) // [(j, a * (matrix ! i ! j)) | j <- [1 .. m]])]
    where (_, m) = mBounds matrix

rowAdd :: (Num a) => Int -> Int -> a -> Matrix a -> Matrix a
rowAdd i j a matrix = matrix // [(j, (matrix ! j) // [(k, (matrix ! j ! k) + a * (matrix ! i ! k)) | k <- [1 .. m]])]
    where (_, m) = mBounds matrix

rowSwap :: Int -> Int -> Matrix a -> Matrix a
rowSwap i j matrix = matrix // [(i, matrix ! j), (j, matrix ! i)]

colMul :: (Num a) => Int -> a -> Matrix a -> Matrix a
colMul i a matrix = matrix // [(j, (matrix ! j) // [(i, a * (matrix ! j ! i))]) | j <- [1 .. n]]
    where (n, _) = mBounds matrix

colAdd :: (Num a) => Int -> Int -> a -> Matrix a -> Matrix a
colAdd i j a matrix = matrix // [(k, (matrix ! k) // [(j, (matrix ! k ! i) * a + (matrix ! k ! j))]) | k <- [1 .. n]]
    where (n, _) = mBounds matrix

colSwap :: Int -> Int -> Matrix a -> Matrix a
colSwap i j matrix = matrix // [(k, (matrix ! k) // [(i, matrix ! k ! j), (j, matrix ! k ! i)]) | k <- [1 .. n]]
    where
        (n, _) = mBounds matrix


identity :: (Num a) => Int -> Matrix a
identity n = array (1, n) [(i, array (1, n) [(j, if i == j then 1 else 0) | j <- [1 .. n]]) | i <- [1 .. n]]

matToList :: Matrix a -> [[a]]
matToList matrix = [[matrix ! i ! j | j <- [1 .. m]] | i <- [1 .. n]]
    where
        (n, m) = mBounds matrix
