-- INTRODUCTION ---------------------------------------------------------------
--
--   Author:
--       Dr-Lord
--   Version:
--       0.2 15/10/2016
--
--   Repository:
--       https://github.com/Dr-Lord/Genome-Sequencing---GUTS-Hackathon-2016
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
--       6 - Unused Functions

---- 1 - IMPORTS AND TYPE DECLARATIONS -----------------------------------------

import Data.List.Split as S (split, endsWith)
import Data.List ((\\), union, delete, intersect)


type Edge = ((Int,Int),String)
data Path = Path {eP1 :: Int, eP2 :: Int, eList :: [Edge]} deriving (Eq, Show, Read)

    -- Make a Path out of an Edge
pathify :: Edge -> Path
pathify (e@((a,b),_)) = Path a b [e]

    -- Prepend or append a (guaranteed to be) linkable edge to a list of ordered edges
correctlyConcat :: Edge -> Path -> Path
correctlyConcat (e@((a,b),s)) el
    | a == eP2 el = Path (eP1 el) b (eList el ++ [e])
    | b == eP1 el = Path a (eP2 el) ((:) e $ eList el)
    | otherwise   = error "This should never happen!!!"
    


---- 2 - TESTING STUFF ---------------------------------------------------------



---- 3 - TO DO -----------------------------------------------------------------

-- FINAL PASS TO REMOVE UNUSED FUNCTIONS


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
--joinAdjacents :: [[String]] -> [String] -> [String]   
--joinAdjacents yPartXed []    = yPartXed
--joinAdjacents yPartXed xPart =
    --(pastedUpBits, remainingIndexedChunks) = foldl step ([], zip [0..] yPartXed) $ unzip linkingSnips 
    --    where step :: ([String],[(Int,[String])]) -> (((Int,Int),String), String) -> ([String],[(Int,[String])])
    --          step (acc,(i,strs)) ((a,b),snip) = (():acc,)
              
              
--    (linkingSnips, remainingPool) = foldl step ([], pool) possSnips
--        where step :: ([((Int,Int),String)], [String]) -> ((Int, Int), String) -> ([((Int,Int),String)], [String]) 
--              step (curRes, curPool) absp@((a,b),snip)
--                | snip `elem` curPool = (absp:curRes, delete snip curPool)
--                | otherwise           = (curRes, curPool)
--    possSnips = [((a,b),snip) | a <- inds, b <- inds, b /= a, let snip = (last $ extrema!!a) ++ (head $ extrema!!b)]
--    inds = [0..(length extrema - 1)]
--    pool = xPart \\ middles
--    (extrema, middles) = map extsAndMids yPartXed
          


---- 5 - OTHER FUNCTIONS -------------------------------------------------------

    -- Link up all possible substrings (using graph analogy: 'chunk-link-chunk's are edges, with chunks being points)
        -- Adapted from my incomplete Travelling Salesman implementation from last year's hackathon, XD
linkUpSubstrings :: [Edge] -> [Path]
linkUpSubstrings = getNextPid []
    where getNextPid :: [Path] -> [Edge] -> [Path]
          getNextPid [] (e:es) = getNextPid [pathify e] es
          getNextPid acc []    = acc
          getNextPid accAcc@(acc:accs) ees@(e:es) =
            case filter (pidsInCommon . getPids) ees of
                []   -> getNextPid ((pathify e):accAcc) es
                [ne] -> getNextPid ((correctlyConcat ne acc):accs) (delete ne es)
                nes  -> error "Is this even possible?"
                    -- NOTE: Is it possible to get two edges with the same point on the same side?
            where getPids ((a,b),_) = [a, b]
                  pidsInCommon = not. null . intersect [eP1 acc, eP2 acc]


    -- Splits the second string into chunks ending with the first string
splitOn :: String -> String -> [String]
splitOn xs = S.split (S.endsWith xs)


    -- Separates edge and middle elements of a list
extsAndMids :: [a] -> ([a],[a])
--extsAndMids [] = ([],[]) -- Should never happen
extsAndMids l@[x]   = (l,[])
extsAndMids l@[x,y] = (l,[])
extsAndMids (x:xs)  = ([x, last xs], init xs)



---- 6 - UNUSED FUNCTIONS ------------------------------------------------------

    -- True if the second string begins with the first string
beginsWith :: String -> String -> Bool
beginsWith xs = (==) xs . take (length xs)


    -- True if the second string ends with the first string
endsWith :: String -> String -> Bool
endsWith xs ys = beginsWith ws zs
    where ws = reverse xs
          zs = take (length xs) $ reverse ys


    -- Delete nth lement of a list
deleteNth :: Int -> [a] -> [a]
deleteNth n xs = [ x | (ix, x) <- zip [0..] xs, ix /= n ]


    -- Returns all elements of a list except edge ones
middleVals :: [a] -> [a]
--middleVals [] = [] -- Should never happen
middleVals [_] = []
middleVals [_,_] = []
middleVals (x:xs) = init xs