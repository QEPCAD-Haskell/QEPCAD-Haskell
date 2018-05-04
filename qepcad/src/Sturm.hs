module Sturm where
import UPolynomial
import Data.Ord
import Data.Ratio (Rational, (%), numerator, denominator)
import qualified Data.Set as S
import Data.Set (Set)
import Math.NumberTheory.ArithmeticFunctions (divisors)

import Numeric (fromRat)


type Sturm r = [UPolynomial r]
type Interval r = (r, r)

width :: Num r => Interval r -> r
inside :: Ord r => r -> Interval r -> Bool
midpoint :: Fractional r => Interval r -> r
width (a, b) = b - a
inside x (a,b) = a < x && x < b
midpoint (a, b) = (a + b)/2


data Root r = Mere r | Algebraic (UPolynomial r) (Interval r)

instance (Num r, Eq r, Show r) => Show (Root r) where
  show (Mere x)   = show x
  show (Algebraic p _) = "⟨root of " ++ show p ++ "⟩"



bisect :: Fractional r => Interval r -> (Interval r, Interval r, r)
bisect i@(a, b) = ((a, c), (c, b), c) where c = midpoint i


var' :: (Eq r, Num r) => UPolynomial r -> Int
var :: (Num r, Eq r) => UPolynomial r -> Interval r -> Int
var' p = length . filter (==(-1)) $ zipWith (*) (init cs) (tail cs)
  where cs = map signum $ coeffs p
var p (a,b) = var' q where
  a' = constant a
  b' = constant b
  d = deg p
  x = uvar
  q = sum [constant c * (a'*x + b')^e * (1+x)^(d-e) | (c,e) <- enumerate p]



sturm :: (Eq r, Fractional r) => UPolynomial r -> Sturm r
sturm p = stseq p (sep p) where
  stseq p 0 = [p]
  stseq p q = p:stseq q (-r) where
    r = remainder p q

vp :: (Eq r, Num r) => Sturm r -> r -> Maybe Int
vp ss x = if 0 `elem` sgs
            then Nothing
            else Just . length . filter id $ zipWith (/=) (tail sgs) (init sgs)
  where sgs = map (signum . eval x) ss


nroots :: (Eq r, Num r) => Sturm r -> Interval r -> Maybe Int
nroots ss (a, b) =
  case (vp ss a,  vp ss b) of
    (Just na, Just nb) -> Just $ na  - nb
    _                  -> Nothing


-- | [SAG, Prop. 1.3]
bound :: UPolynomial Rational -> Int
bound p = ceiling $ 0.1 + maximum [(fromIntegral d * abs (c'/a0))**(1.0/fromIntegral i) |
  (c,e) <- enumerate p,
  let i = d - e,  0 < i, let c' = fromRat c] where
 d = deg p
 a0 = fromRat $ lc p

bound' :: UPolynomial Rational -> Float
bound' p = 1 + maximum [fabs c/fabs c0 | c <- coeffs (trunc p)]
  where c0 = lc p
        fabs = abs . fromRational


refine :: Root Rational -> Root Rational
refine r@(Algebraic p iv) =
  case S.toList $ qroots p iv of
    [rq] -> Mere rq
    _       -> r
  -- | from https://en.wikipedia.org/wiki/Rational_root_theorem
refine r      = r

ptoZ :: UPolynomial Rational -> UPolynomial Integer
ptoZ p =  fmap (\q -> let n = numerator q
                          d = denominator q
                      in n * m `div` d ) p
   where
     m :: Integer
     m = foldl1 lcm . map denominator $ coeffs p
qroots :: UPolynomial Rational -> Interval Rational -> Set Rational
qroots p iv  =
  if a0 == 0
  then S.insert 0 $ qroots (p `qquot` t) iv
  else S.fromList . filter ((==0) . flip eval p) $ rs
    where p' :: UPolynomial Integer
          p' = ptoZ p
          t = uvar
          an = head (coeffs p')
          a0 = last (coeffs p')
          factors = S.toList . divisors
          rs = [r | p <- factors a0, q <- factors an,
                    s <- [1, -1],
                    let r = s * p % q, r `inside` iv]


roots :: UPolynomial Rational -> [Root Rational]
roots p = wrap (-b, b) where
  b = fromIntegral $ bound p
  nr = nroots (sturm p)
  wiggle w x = if w < 1%10000
    then error $ "too small at " ++ show x
    else let iv = (x - w, x + w)
         in  case nr iv of
            Just 1 -> iv
            _      -> wiggle (w/2) x
  wrap iv@(a, b) =
    case nr iv of
      Nothing | eval a p == 0 ->
                let (_,a') = wiggle (width iv/3) a in
                  if eval b p == 0 then
                    let (b',_) = wiggle (width iv/3) b in
                      Mere a : wrap (a',b') ++ [Mere b]
                    else Mere a : wrap (a', b)
              | eval b p == 0 ->
                let (b',_) = wiggle (width iv/3) b in wrap (a, b') ++ [Mere b]
              | otherwise -> error "cannot be"
      Just 0  -> []
      Just 1  -> [Algebraic p (a,b)]
      Just _  -> let (i,j, c) = bisect (a,b)
            in if eval c p == 0
               then wrap i ++ Mere c :  wrap j
               else wrap i ++ wrap j

abs1 :: Num r =>  UPolynomial Rational -> Float
abs1 = fromRat . sum . map abs . coeffs

type Approx2 = Int -> UPolynomial Rational

-- | Assume all roots lie in a disk of radius 1/2 centered at 0


-- newtype Cell r = Cell (Either r (Int, UPolynomial r)) deriving Show


p, x :: UPolynomial Rational
x = uvar
p = x^3 - 3*x +1
