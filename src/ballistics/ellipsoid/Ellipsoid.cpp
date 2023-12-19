#include "Ellipsoid.hpp"

namespace Ballistics::Ellipsoid
{
    [[nodiscard]] scalar Ellipsoid::calcHeight(const Vector3d &position) const noexcept
    {
        scalar R = std::hypot(position(0), position(1));
        scalar height = std::hypot(R, position(2));  // Distance to center of earth
        if (height > maxRad_) {
            return height;
        } else if (e4a_ == 0) {
            height -= a_;
            return height;
        } else {
            // Treat prolate spheroids by swapping R and Z here and by switching
            // the arguments to phi = atan2(...) at the end.
            scalar p = R / a_;
            p *= p;
            scalar q = position(2) / a_;
            q *= q;
            q *= e2m_;
            scalar r = (p + q - e4a_) / 6;
            if (f_ < 0) std::swap(p, q);
            if (!(e4a_ * q == 0 && r <= 0)) {
                // Avoid possible division by zero when r = 0 by multiplying
                // equations for s and t by r^3 and r, resp.
                scalar S = e4a_ * p * q / 4;  // S = r^3 * s
                scalar r2 = r * r, r3 = r * r2, disc = S * (2 * r3 + S);
                scalar u = r;
                if (disc >= 0) {
                    scalar T3 = S + r3;
                    // Pick the sign on the sqrt to maximize abs(T3).  This minimizes
                    // loss of precision due to cancellation.  The result is unchanged
                    // because of the way the T is used in definition of u.
                    T3 += T3 < 0 ? -std::sqrt(disc) : std::sqrt(disc);  // T3 = (r * t)^3
                    // N.B. cbrt always returns the real root.  cbrt(-8) = -2.
                    scalar T = std::cbrt(T3);  // T = r * t
                    // T can be zero; but then r2 / T -> 0.
                    u += T + (T != 0 ? r2 / T : 0);
                } else {
                    // T is complex, but the way u is defined the result is real.
                    scalar ang = std::atan2(std::sqrt(-disc), -(S + r3));
                    // There are three possible cube roots.  We choose the root which
                    // avoids cancellation.  Note that disc < 0 implies that r < 0.
                    u += 2 * r * std::cos(ang / 3);
                }
                scalar v = std::sqrt(u * u + e4a_ * q),  // guaranteed positive
                // Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
                // e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
                uv = u < 0 ? e4a_ * q / (v - u) : u + v,  // u+v, guaranteed positive
                // Need to guard against w going negative due to roundoff in uv - q.
                w = std::max(static_cast<scalar>(0.), e2a_ * (uv - q) / (2 * v)),
                // Rearrange expression for k to avoid loss of accuracy due to
                // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
                k = uv / (std::sqrt(uv + w * w) + w), k1 = f_ >= 0 ? k : k - e2_, k2 = f_ >= 0 ? k + e2_ : k,
                        d = k1 * R / k2;
                height = (1 - e2m_ / k1) * std::hypot(d, position(2));
            } else {  // e4 * q == 0 && r <= 0
                // This leads to k = 0 (oblate, equatorial plane) and k + e^2 = 0
                // (prolate, rotation axis) and the generation of 0/0 in the general
                // formulas for phi and h.  using the general formula and division by 0
                // in formula for h.  So handle this case by taking the limits:
                // f > 0: z -> 0, k      ->   e2_ * sqrt(q)/sqrt(e4 - p)
                // f < 0: R -> 0, k + e2_ -> - e2_ * sqrt(q)/sqrt(e4 - p)
                scalar zz = std::sqrt((f_ >= 0 ? e4a_ - p : p) / e2m_), xx = std::sqrt(f_ < 0 ? e4a_ - p : p),
                        H = std::hypot(zz, xx);
                height = -a_ * (f_ >= 0 ? e2m_ : 1) * H / e2a_;
            }
        }
        return height;
    }
}
