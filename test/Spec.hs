import SNF
import Matrix
import Control.Monad

examples :: [Matrix Integer]
examples = [
    example
    ]


main :: IO ()
main = do
    forM_ examples $ \ex -> do
        let (a, b, c) = smithNormalForm ex
        putStrLn $ "Decomposition equality is valid: " ++ show (mulMatrix a (mulMatrix ex c) == b)
