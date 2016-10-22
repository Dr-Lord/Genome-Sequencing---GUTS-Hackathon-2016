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

import Data.List.Split (split, endsWith, endsWithOneOf)
import Data.List (isSuffixOf, inits, (\\), union, delete, intersect, nubBy, sortBy, tails)
import Data.Function (on)



---- 2 - TESTING STUFF ---------------------------------------------------------



---- 3 - TO DO -----------------------------------------------------------------

-- FINAL PASS TO REMOVE UNUSED FUNCTIONS AND IMPORTS


---- 4 - MAIN FUNCTIONS --------------------------------------------------------

    -- Main program
main = do
    -- Read in files and extract lists of strings
    let pathPrefix = ".\\AlienDNAChallenge\\genomePieces\\"
    let pathInfix = "_digest_"
    let genomes = ["1k", "10k", "100k", "1m"]
    let genomeLengths = [1000, 10000, 100000, 1000000]
    let tags = ["BC", "DE", "EDA", "DFAD"]

        -- SET THIS VALUE TO CHANGE GENOME
    let genomeIndex = 1

    partitionsStrs <- mapM readFile $ [pathPrefix ++ x ++ pathInfix ++ y | let x = genomes!!genomeIndex, y <- tags]
        -- Interpret the files as lists of strings
    let partitions = map (read :: String -> [String]) partitionsStrs
    
    putStrLn $ "Chosen genome length: " ++ genomes!!genomeIndex
    putStrLn . ("Total lengths of partitions: " ++) . show $ map (sum . map length) partitions
    
    
        -- Remove the elements containing the genome end from each partition and isolate the longest endstring
    let (newPartitions,endString,n,endStrings) = removeGenomeEnd partitions tags
    
    putStrLn . ("Total lengths of endStrings: " ++) . show $ map length endStrings
    putStrLn . ("Total lengths of newPartitions: " ++) . show $ map (sum . map length) newPartitions

    putStrLn . ("Index of longest endstring: " ++) . show $ n
    putStrLn . ("Initial endstring length: " ++) . show $ length endString
    
    
        -- Partition each partition further by splitting on all tags
    let allFullySplit = [ [splitOnAll (delete (tags!!x) tags) y | y <- newPartitions!!x] | x <- [0..3] ]
        -- Length one sublists cannot contribute to the sequencing
    let joinedPools = filter ((/= 1) . length) $ concat allFullySplit 
        -- Split the endstring on all tags too
    let startList = splitOnAll tags endString
    
    --putStrLn . ("Initial endstring: " ++) . show $ startList

        -- Sequence the genome and return it as a list of substring (ending in any tag); also return the remaining pool
    --let (res, newPool) = alignAndJoinWithOverlapFast 1 startList joinedPools
    let (res, newPool) = alignAndJoinDeterministic (genomeLengths!!genomeIndex) startList joinedPools
    
    putStrLn "\n\nRESULTS:\n"
    putStrLn . ("Length of final result: " ++) . show $ sum $ map length res
    --putStrLn . ("Final result: " ++) . show $ res
    putStrLn . ("Pool elements containing first element of result: " ++) . show $ filter ("FFFBDCEABDBFDACFCDFCCBACEFEADCBAEACAAFBADCCADCAAEFEEFACDE" `elem`) newPool

    putStrLn . ("\n\nThe result ends with the endstring: " ++ ) . show $ isSuffixOf startList res
    
    putStrLn . ("\nRemaining pool size: " ++) . show $ length newPool
    putStrLn . ("Remaining pool total length: " ++) . show $ sum $ map (sum . map length) newPool




    -- Progressively align and merge lists of strings to a given end sequence
    -- NOTE: This is a deterministic version, making sure that there is always only one option to prepend to the accumulator
alignAndJoinDeterministic :: Int -> [String] -> [[String]] -> ([String],[[String]])
alignAndJoinDeterministic targetLength startLs joinedPools = alAndJoin (startLs, joinedPools) $ overlapsAndInits maxSubLen startLs
    where maxSubLen = maximum $ map length joinedPools
          alAndJoin :: ([String],[[String]]) -> [(Int,[String])]-> ([String],[[String]])
          alAndJoin out               []                     = out
          alAndJoin out@(acc,curPool) ((overlap,curInit):xs) = case filter (curInit `isSuffixOf`) curPool of
                [match] -> let newAcc = match ++ drop overlap acc in if length newAcc >= targetLength
                            then (newAcc, delete match curPool)
                            else alAndJoin (newAcc, delete match curPool) $ overlapsAndInits maxSubLen newAcc
                _       -> alAndJoin out xs


    -- Return all inits shorter than maxLen and their lengths of a list, except for the first one, which would be empty
overlapsAndInits :: Int -> [a] -> [(Int,[a])]
overlapsAndInits maxLen = take maxLen . tail . zip [0..] . inits







    -- Progressively align and merge lists of strings to a given end sequence
    -- NOTE: This is a fast version relying on the low statistical chance of no false positives
alignAndJoinWithOverlapFast :: Int -> [String] -> [[String]] -> ([String],[[String]])
alignAndJoinWithOverlapFast overlap startLs joinedPools = alAndJoin (startLs, joinedPools) 0 joinedPools
    where alAndJoin :: ([String],[[String]]) -> Int -> [[String]]-> ([String],[[String]])
          alAndJoin out               0        []     = out
          alAndJoin out@(_,  curPool) _        []     = alAndJoin out 0 curPool
          alAndJoin out@(acc,curPool) joinsNum (x:xs) = case areAlignableWithOverlapFast overlap x acc of
            Just joined -> alAndJoin (joined, delete x curPool) (joinsNum+1) xs
            Nothing     -> alAndJoin out joinsNum xs
            

    -- Check whether 2 lists of strings are alignable (by considering overlap elements) and return them overlapped if they are
    -- NOTE: This is a fast version relying on the low statistical chance of no false positives
areAlignableWithOverlapFast :: Int -> [String] -> [String] -> Maybe [String]
areAlignableWithOverlapFast overlap x rest = foldl step Nothing noExtraTails
    where noExtraTails = take (length x) noEndTails
          noEndTails = take (length rest + 1 - overlap) $ tails rest -- Note that length rest - overlap < 0 short-circuits the fold 
          step :: Maybe [String] -> [String] -> Maybe [String]
          step acc aTail
            | aligns && better = Just $ x ++ drop overlap aTail
            | otherwise        = acc
                where aligns = isSuffixOf (take overlap aTail) x
                      better = case acc of
                          Just l  -> length x > length l
                          Nothing -> True



---- 5 - OTHER FUNCTIONS -------------------------------------------------------

    -- Remove the end of the genome from all partitions and preserve the longest one
removeGenomeEnd :: [[String]] -> [String] -> ([[String]], String, Int, [String])
removeGenomeEnd allPartitions allTags = foldr step ([],"",0,[]) [0..3]
    where step n (aPs,curEndStr,curN,endStrs)
            | not $ null end = if length s > length curEndStr then ((delete s aPart):aPs, s, n, s:endStrs) else ((delete s aPart):aPs, curEndStr, curN, s:endStrs)
            | otherwise      = (aPart:aPs, curEndStr, curN, endStrs)
            where [s] = end
                  end = filter (not . isSuffixOf aTag) aPart
                  aTag = allTags!!n
                  aPart = allPartitions!!n
                                          

    -- Splits the second string into chunks ending with the first string
splitOn :: String -> String -> [String]
splitOn xs = split (endsWith xs)


    -- Splits the second string into chunks ending with any of the strings in the first list
splitOnAll :: [String] -> String -> [String]
splitOnAll []       toSplit = [toSplit]
splitOnAll (xs:xss) toSplit = foldl step (splitOn xs toSplit) xss
    where step :: [String] -> String -> [String]
          step acc xs' = concatMap (splitOn xs') acc



---- 6 - UNUSED FUNCTIONS ------------------------------------------------------

{-
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
-}

{-
    let xPart = partitions!!0
    let yPart = partitions!!1
    let yPartXed = map (splitOn (tags!!0)) yPart
    
    let (joinedPaths, remainingPool) = joinAdjacents (tags!!0, tags!!1) yPartXed xPart
    let joinedStrings = map (stringify yPart) joinedPaths
    
    putStrLn "\n\nRESULTS:\n"
    putStrLn . show $ res
    
    putStrLn . show $ joinedStrings
    putStrLn $ show remainingPool
    putStrLn . show $ sum . map length $ joinedStrings ++ remainingPool
    putStrLn . show $ filter (\x-> length x == 97) remainingPool
    
     if remainingPool == xPart then putStrLn "No progress was made" else putStrLn "SOME PROGRESS WAS MADE!!"
    
    
    
    let orderedStrs = sortBy (compare `on` length) joinedStrings
    putStrLn $ "Longest substring length: " ++ (show . length $ last orderedStrs)
    putStrLn $ "Shortest substring length: " ++ (show . length $ head orderedStrs)

    putStrLn $ "Longest substring: " ++ (show $ last orderedStrs)
-}

{-
    -- True if the second list begins with the first list
beginsWithSub :: Eq a => [a] -> [a] -> Bool
beginsWithSub xs = (==) xs . take (length xs)


    -- True if the second list ends with the first list
endsWithSub :: Eq a => [a] -> [a] -> Bool
endsWithSub xs ys = beginsWithSub ws zs
    where ws = reverse xs
          zs = take (length xs) $ reverse y


    -- Short circuit nub for Edge lists
edgeNub :: [Edge] -> [Edge]
edgeNub = nubBy ((==) `on` fst)


    -- Separates edge and middle elements of a list
extsAndMids :: [a] -> ([a],[a])
--extsAndMids [] = ([],[]) -- Should never happen
extsAndMids l@[x]   = (l,[])
extsAndMids l@[x,y] = (l,[])
extsAndMids (x:xs)  = ([x, last xs], init xs)

    -- Split partition-Y on tag-X then try and merge partitions X and Y
--merge :: String -> [String] -> [String] -> [String]
--merge xTag xPart yPart = 
    

    -- Join a partitioned-by-tag-X partition-Y element with adjacent partition-Y elements
    -- by looking for split partition-X elements on bondaries
    -- Returns a list of paths
        -- NOTE: Ignore and redo at end instances where there is more than one possible split
        -- partition-X element on a boundary (the set of remaining partition-X elements decreases every step)
        -- NOTE: Cater for possibility of current incompletability (all possibility sets are greater than one for a full pass)
        -- PERHAPS it would be useful to sort the yPartXed list by decreasing total original length of chunks
joinAdjacents :: (String,String) -> [[String]] -> [String] -> ([Path], [String])
joinAdjacents (tA,tB) yPartXed xPart = (linkedBits, remainingPool) 
    where linkedBits = linkUpSubstrings linkingSnips 
          (linkingSnips, remainingPool) = foldl step ([], pool) possSnips
            where step :: ([Edge], [String]) -> Edge -> ([Edge], [String]) 
                  step (curRes, curPool) absp@(_,snip)
                    | snip `elem` curPool = (absp:curRes, delete snip curPool)
                    | otherwise           = (curRes, curPool)
          possSnips = [((a,b),snip) | a <- inds, b <- inds, b /= a, let snip = (last $ extrema!!a) ++ (head $ extrema!!b)]
          --possSnips = edgeNub [((a,b),snip) | a <- inds, b <- inds, b /= a, let snip = (last $ extrema!!a) ++ (head $ extrema!!b)]
          inds = [0..length extrema - 1]
          pool = (xPart \\) . concat $ extrema ++ middles
          --possBeginnings = filter (\x-> endsWithSub tA x || endsWithSub tB x) $ xPart `intersect` (concat extrema)
          (extrema, middles) = unzip $ map extsAndMids yPartXed


    -- Turn a Path into a String
stringify :: [String] -> Path -> String
stringify yPart (Path _ _ (((w,z),str):es)) = foldl joinUp (yPart!!w ++ str ++ yPart!!z) es
    where joinUp acc ((_,b),s) = acc ++ s ++ yPart!!b


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




    -- Delete nth lement of a list
deleteNth :: Int -> [a] -> [a]
deleteNth n xs = [ x | (ix, x) <- zip [0..] xs, ix /= n ]


    -- Returns all elements of a list except edge ones
middleVals :: [a] -> [a]
--middleVals [] = [] -- Should never happen
middleVals [_] = []
middleVals [_,_] = []
middleVals (x:xs) = init xs

-}