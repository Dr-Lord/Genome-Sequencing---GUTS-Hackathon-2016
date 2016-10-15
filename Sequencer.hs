-- INTRODUCTION ---------------------------------------------------------------
--
--   Author:
--       Dr-Lord
--   Version:
--       0.1 14/10/2016
--
--   Repository:
--       https://github.com/Dr-Lord/
--
--   Description:
--      This program takes in multiple randomised partitions of genomes generated
--      by splitting it on different substrings.
--      It recombines the whole genome by comparing each set of substrings looking
--      for alignments.
--
--   Sections:
--       1 - Imports and Type declarations
--       2 - Testing Stuff
--       3 - To Do
--       4 - Main Functions
--       5 - Other Functions

---- 1 - IMPORTS AND TYPE DECLARATIONS -----------------------------------------

import Data.List.Split as S (split, endsWith)
import Data.List ((\\), union, delete)



---- 2 - TESTING STUFF ---------------------------------------------------------



---- 3 - TO DO -----------------------------------------------------------------



---- 4 - MAIN FUNCTIONS --------------------------------------------------------

    -- Main program
main = do
    -- Read in files and extract lists of strings
    let prefix = ".\\AlienDNAChallenge\\genomePieces\\1k_digest_"
    let tags = ["BC", "DE", "EDA", "DFAD"]
    partitionsStrs <- mapM readFile $ map (prefix++) tags
    let partitions = map (read :: String -> [String]) partitionsStrs
    
    -- Check all of same length
    let allSame = map (sum . map length) partitions
    putStrLn $ show allSame
    
    -- Look for one partitions' parts inside all other partitions' parts
        -- Split at all tags (one at a time) and compare/expand largest matches
        -- Perhaps all the pairings are not necessary; just incrementally work on the confirmed substrings
    let pairings = [(x,y) | x <- [0..3], y <- [0..3], x /= y]
    putStrLn $ show pairings
    
    
    
   
    -- Split partition-Y on tag-X then try and merge partitions X and Y
--merge :: String -> [String] -> [String] -> [String]
--merge xTag xPart yPart = 
    

    -- Join a partitioned-by-tag-X partition-Y element with adjacent partition-Y elements
    -- by looking for split partition-X elements on bondaries
        -- NOTE: Ignore and redo at end instances where there is more than one possible split
        -- partition-X element on a boundary (the set of remaining partition-X elements decreases every step)
        -- NOTE: Cater for possibility of current incompletability (all possibility sets are greater than one for a full pass)
        -- PERHAPS it would be useful to sort the yPartXed list by decreasing total original length of chunks
joinAdjacents :: [[String]] -> [String] -> [String]
joinAdjacents yPartXed []    = yPartXed
joinAdjacents yPartXed xPart = joinAdjacents' yPartXed xPart [] $ length yPartXed
    where joinAdjacents' []       xp res _ = joinAdjacents res xp
          joinAdjacents' (e:ypxd) xp res n = joinAdjacents' ypxd' (xp \\ (middleVals e)) res' n'
            | length lAdjs == 0 && length rAdjs == 0 = 
            | 
            where ypxd' = ypxd \\ (le `union` re)
                  res'  = [[lAdj] ++ e ++ [rAdj]]:(delete e res)
                  n' = length ypxd' + length res'
                  lAdjs = filter () allChunkss
                  rAdjs = filter () allChunkss
                  (l,r) = (head e, last e)
                  allChunkss = ypxd ++ res
                
                
    
    pool = xPart \\ middles
    (extrema, middles) = map extsAndMids yPartXed
        -- And then just combine the exts in order to get pool items
          


---- 5 - OTHER FUNCTIONS -------------------------------------------------------

    -- True if the second string begins with the first string
beginsWith :: String -> String -> Bool
beginsWith xs = (==) xs . take (length xs)


    -- True if the second string ends with the first string
endsWith :: String -> String -> Bool
endsWith xs ys = beginsWith ws zs
    where ws = reverse xs
          zs = take (length xs) $ reverse ys


    -- Splits the second string into chunks ending with the first string
splitOn :: String -> String -> [String]
splitOn xs = S.split (S.endsWith xs)


    -- Returns all elements of a list except edge ones
middleVals :: [a] -> [a]
--middleVals [] = [] -- Should never happen
middleVals [_] = []
middleVals [_,_] = []
middleVals (x:xs) = init xs


    -- Separates edge and middle elements of a list
extsAndMids :: [a] -> ([a],[a])
--extsAndMids [] = ([],[]) -- Should never happen
extsAndMids l@[x]   = (l,[])
extsAndMids l@[x,y] = (l,[])
extsAndMids (x:xs)  = ([x, last xs], init xs)