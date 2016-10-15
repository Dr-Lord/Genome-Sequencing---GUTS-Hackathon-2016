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