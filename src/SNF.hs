{-# LANGUAGE ForeignFunctionInterface #-}
module SNF (smithNormalForm) where

import Data.Array ((!))
import System.IO.Unsafe (unsafePerformIO)
import Foreign.C
import Foreign.Ptr
import Foreign.Marshal.Array

import Matrix
import Foreign.Storable

foreign import ccall "csmith_normal_form" cSmithNormalForm :: Ptr CLLong -> IO (Ptr CLLong)

smithNormalForm :: Matrix Integer -> [Integer]
smithNormalForm matrix = unsafePerformIO $ do
    let (n, m) = mbounds matrix
    cinput <- mallocArray (8 * (n * m + 2)) :: IO (Ptr CLLong)
    pokeArray cinput (fromIntegral n : fromIntegral m : [fromIntegral (matrix ! (i, j)) | i <- [1 .. n], j <- [1 .. m]])
    coutput <- cSmithNormalForm cinput
    k <- peek coutput

    map fromIntegral <$> (peekArray (fromIntegral k) (plusPtr coutput 8) :: IO [CLLong])
