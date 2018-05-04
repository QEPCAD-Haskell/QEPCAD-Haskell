module Sturm where
import UPolynomial
import Data.Ord
import Data.Ratio (Rational, (%), numerator, denominator)
import Math.NumberTheory.ArithmeticFunctions (divisors)

import Numeric (fromRat)


type Sturm r = [UPolynomial r]
type Interval r = (r, r)

inside :: r -> Interval r -> Bool
inside x (a,b) = a < x && x < b

data Root r = Mere r | Algebraic (UPolynomial r) (Interval r)

instance (Num r, Eq r, Show r) => Show (Root r) where
  show (Mere x)   = show x
  show (Algebraic p _) = "<root of " ++ show p ++ ">"



midpoint :: Fractional r => Interval r -> r
bisect :: Fractional r => Interval r -> (Interval r, Interval r)
midpoint (a, b) = (a + b)/2
bisect i@(a, b) = ((a, c), (c, b)) where c = midpoint i


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

vp :: (Eq r, Num r) => Sturm r -> r -> Int
vp ss x = length $ filter id $ zipWith (/=) (tail sgs) (init sgs) where
      sgs = filter (/=0) . map (signum . eval x) $ ss

nroots :: (Eq r, Num r) => Sturm r -> Interval r -> Int
nroots ss (a, b) = vp ss a - vp ss b

-- | [SAG, Prop. 1.3]
bound :: UPolynomial Rational -> Int
bound p = ceiling $ 0.1 + maximum [(fromIntegral d * abs (c'/a0))**(1.0/fromIntegral i) |
  (c,e) <- enumerate p,
  let i = d - e,  0 < i, let c' = fromRat c] where
 d = deg p
 a0 = fromRat $ lc p

bound' :: UPolynomial Rational -> Float
bound' p = 1 + maximum [fabs(c)/fabs(c0) | c <- coeffs (trunc p)]
  where c0 = lc p
        fabs = abs . fromRational

isolate' :: (Eq r, Fractional r)  => Sturm r -> Interval r -> [Interval r]
isolate' st i@(a,b) =
  if nroots st i <= 1
  then [i]
  else let (l,r) = bisect i
       in isolate' st l ++ isolate' st r
isolate'' :: (Eq r, Fractional r) => Interval r -> [Sturm r] -> [Interval r]
isolate'' i = foldr (\s -> concat . map (isolate' s)) [i]


refine :: Root Rational -> Root Rational
refine (Algebraic p (a,b)) = undefined where
  -- | from https://en.wikipedia.org/wiki/Rational_root_theorem
refine r      = r

ptoZ :: UPolynomial Rational -> UPolynomial Integer
ptoZ p =  fmap (\q -> let n = numerator q
                         d = denominator q
                     in n * m `div` d ) p
   where
     m :: Integer
     m = foldl1 lcm . map denominator $ coeffs p


qrootIn :: UPolynomial Integer -> Interval Rational -> Maybe Rational
qrootIn p iv | last (coeffs p) == 0 =
  if 0 `inside` iv
  then Just 0
  else qrrotIn (p `division` uvar) iv
qrootIn p iv                     =
    case filter ((==0) . flip eval p) rs of
      [] -> Nothing
      [r] -> Just r
      _   -> error "Not isolating"
    where
      an = head (coeffs p)
      a0 = tail (coeffs p)
      rs = [r | p <- divisors a0, q <- divisors an, s <- [1,(-1)], let r = (s*p%q), r `inside` iv]


roots :: UPolynomial Rational -> [Root Rational]
roots p = wrap (-b, b) where
  b = fromIntegral $ bound p
  nr a b = nroots (sturm p) (a,b)
  wrap (a, b) = case nr a b of
    0 -> []
    1 -> [Algebraic p (a,b)]
    _ -> let (i,j) = bisect (a,b)
         in wrap i ++ wrap j

abs1 :: Num r =>  UPolynomial Rational -> Float
abs1 = fromRat . sum . map abs . coeffs

type Approx2 = Int -> UPolynomial Rational

-- | Assume all roots lie in a disk of radius 1/2 centered at 0


-- newtype Cell r = Cell (Either r (Int, UPolynomial r)) deriving Show


p, x :: UPolynomial Rational
x = uvar
p = x^3 - 3*x +1
