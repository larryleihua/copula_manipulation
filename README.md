# copula manipulation

Based on some bivariate copulas, one may need to rotate, reflect, mix, and distort them to construct a new bivariate copula. This project is for building C++ functions to handle these model constructions. The license is GPL-3

- Rotation: Rotation can be achieved by C'(u, v) = 0.5*Cbar(1-u, 1-v) + 0.5*C(u,v), which is a rotation of C(u,v) around the point (0.5, 0.5). The basic functions for C'(u,v) are given for some cases that are based on commonly-used bivariate copulas, such as, student t, Clayton, Joe, Gumbel, etc.

- Non-exchangeability: The Khoudraji's transformation can be used.
