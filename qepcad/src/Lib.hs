module Lib where
import CAD
import qualified Data.Set as S
import Polynomial

examples =
    let [x,y,z] = map single "xyz"
        p1 = [x^2 + y^2 + z^2 -2]
        p2 = [x * y * z - x^2 - y^2]
        p3 = [p, q, r] where
          p = z
          q = x^2 + y^2 - z
          r = x^2 + y^2 - 2 * x
    in map S.fromList [p1, p2, p3]
