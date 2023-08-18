module SNF (smithNormalForm) where

import Data.List ( find )
import Data.Array ( (!), bounds )
import Control.Monad.State.Lazy ( when, forM_, modify, execState, MonadState(get), State )
import Data.Maybe ( isJust )

import Matrix (Matrix, ElemOp(..))
import qualified Matrix

data SmithState = SmithState {
  smithMatrix :: Matrix Integer,
  rowOps :: [ElemOp],
  colOps :: [ElemOp]
}
  deriving (Show, Eq)

toSmithState :: Matrix Integer -> SmithState
toSmithState matrix = SmithState matrix [] []

getMat :: State SmithState (Matrix Integer)
getMat = smithMatrix <$> get

euclid :: Integer -> Integer -> (Integer, Integer, Integer)
euclid a b = (x, y, d)
  where
    (x, y, d) = go a b 1 0 0 1
    go a 0 x y _ _ = (x, y, a)
    go a b x y x' y' = go b (a `mod` b) x' y' (x - q * x') (y - q * y')
      where
        q = a `div` b

rowSwap :: Int -> Int -> SmithState -> SmithState
rowSwap i j (SmithState matrix rops cops) = SmithState (Matrix.rowSwap i j matrix) (Swap i j : rops) cops

rowAdd :: Int -> Int -> Integer -> SmithState -> SmithState
rowAdd i j a (SmithState matrix rops cops) = SmithState (Matrix.rowAdd i j a matrix) (Add i j a : rops) cops

rowMul :: Int -> Integer -> SmithState -> SmithState
rowMul i a (SmithState matrix rops cops) = SmithState (Matrix.rowMul i a matrix) (Mul i a : rops) cops

colSwap :: Int -> Int -> SmithState -> SmithState
colSwap i j (SmithState matrix rops cops) = SmithState (Matrix.colSwap i j matrix) rops (Swap i j : cops)

colAdd :: Int -> Int -> Integer -> SmithState -> SmithState
colAdd i j a (SmithState matrix rops cops) = SmithState (Matrix.colAdd i j a matrix) rops (Add i j a : cops)

colMul :: Int -> Integer -> SmithState -> SmithState
colMul i a (SmithState matrix rops cops) = SmithState (Matrix.colMul i a matrix) rops (Mul i a : cops)

  -- 1. find pivot
  -- 2. swap rows
  -- 3. find the rest of the non-zero elements in the column
  -- 4. add rows according to bezout
  -- 5. wipe out the rest of the column
diagonalForm :: Matrix Integer -> SmithState
diagonalForm matrix = execState (go 1) $ toSmithState matrix where
  (1, n) = bounds matrix
  (1, m) = bounds (matrix ! 1)

  go :: Int -> State SmithState ()
  go t = do
    matrix <- getMat
    case fmap (\(Just x, i) -> (x, i)) $ filter (isJust . fst) [(findPivot matrix j t, j) | j <- [t..m]] of
      [] -> return ()
      ((i, j) : _) -> do
        modify $ rowSwap t i -- ensure that the pivot is on (t, t)
        modify $ colSwap t j

        -- make the pivot equal to the gcd of the column, by performing elementary operations
        forM_ [t+1 .. n] $ \i -> euclidCol i t

        -- make the pivot equal to the gcd of the row, by performing elementary operations
        forM_ [t+1 .. m] $ \j -> euclidRow j t


        -- wipe out the rest of the column
        forM_ [t+1 .. n] $ \i -> do
          matrix <- getMat
          let a = matrix ! i ! t `div` matrix ! t ! t
          modify $ rowAdd t i (-a)

        -- wipe out the rest of the row
        forM_ [t+1 .. m] $ \j -> do
          matrix <- getMat
          let a = matrix ! t ! j `div` matrix ! t ! t
          modify $ colAdd t j (-a)

        go (t + 1)

  euclidCol :: Int -> Int -> State SmithState ()
  euclidCol i t = do
    matrix <- getMat
    let (a, b) = (matrix ! t ! t, matrix ! i ! t)
    when (a < 0) $ modify $ rowMul t (-1)
    when (b < 0) $ modify $ rowMul i (-1)
    when (b /= 0) $ do
      modify $ rowAdd i t (-a `div` b)
      modify $ rowSwap t i
      euclidCol i t

  euclidRow :: Int -> Int -> State SmithState ()
  euclidRow i t = do
    matrix <- getMat
    let (a, b) = (matrix ! t ! t, matrix ! t ! i)
    when (a < 0) $ modify $ colMul t (-1)
    when (b < 0) $ modify $ colMul i (-1)
    when (b /= 0) $ do
      modify $ colAdd i t (-a `div` b)
      modify $ colSwap t i
      euclidRow i t

  findPivot :: Matrix Integer -> Int -> Int -> Maybe Int
  findPivot matrix column t = find (\i -> matrix ! i ! column /= 0) [t .. n]

smithNormalForm :: Matrix Integer -> (Matrix Integer, Matrix Integer, Matrix Integer)
smithNormalForm matrix = (leftMat, diagMat, rightMat) where
  df = diagonalForm matrix
  (_, n) = bounds matrix
  (_, m) = bounds (matrix ! 1)

  go :: Int -> State SmithState ()
  go t = when (t <= min n m) $ do
    matrix <- getMat
    let nonzero = filter (\i -> matrix ! i ! i /= 0) [t .. min n m]
    case nonzero of
      [] -> return ()
      (i : _) -> do
        modify $ rowSwap t i -- ensure that the pivot is on (t, t)
        modify $ colSwap t i

    when (matrix ! t ! t /= 0) $ do
      -- make the (t, t) cell equal to the gcd of everything that follows it on the diagonal
      forM_ [t+1 .. min n m] $ \i -> do
        matrix <- getMat
        let a = matrix ! t ! t
        let b = matrix ! i ! i
        let (x, y, d) = euclid a b

        when (b `mod` a /= 0) $ do
          modify $ rowMul t x
          modify $ rowAdd i t 1
          modify $ colAdd i t 1
          modify $ colAdd t i $ negate (b `div` d)
          modify $ rowAdd t i $ negate (b * y `div` d)

    -- ensure that the entry is non-negative
    when (matrix ! t ! t < 0) $ modify $ rowMul t (-1)
    go (t + 1)

  resState = execState (go 1) df
  (leftMat, diagMat, rightMat) = stateToDecomp resState

stateToDecomp :: SmithState -> (Matrix Integer, Matrix Integer, Matrix Integer)
stateToDecomp resState = (lMat, dMat, rMat) where
  diagMat = smithMatrix resState
  lMat = foldr (\op acc -> case op of
    Swap i j -> Matrix.rowSwap i j acc
    Add i j a -> Matrix.rowAdd i j a acc
    Mul i a -> Matrix.rowMul i a acc
    ) (Matrix.identity n) (rowOps resState)
  rMat = foldr (\op acc -> case op of
    Swap i j -> Matrix.colSwap i j acc
    Add i j a -> Matrix.colAdd i j a acc
    Mul i a -> Matrix.colMul i a acc
    ) (Matrix.identity m) (colOps resState)
  dMat = smithMatrix resState

  (_, n) = bounds diagMat
  (_, m) = bounds (diagMat ! 1)

