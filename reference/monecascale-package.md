# monecascale: Scalable Backends for MONECA Mobility Clustering

Sibling package to moneca providing scalable clustering backends whose
output is consumed unchanged by moneca's analysis and plotting stack.
The load-bearing theoretical observation is that moneca's relative-risk
matrix `RR = O / E` is the degree-corrected stochastic block model
(DC-SBM) residual under identity block interaction; SBM inference is
therefore a principled drop-in for moneca's clique-enumeration step, and
it scales to settings (10^6+ nodes) where clique enumeration is not
viable.

## See also

Useful links:

- <https://github.com/gmontaletti/monecascale>

- Report bugs at <https://github.com/gmontaletti/monecascale/issues>

## Author

**Maintainer**: Giampaolo Montaletti <giampaolo.montaletti@gmail.com>
([ORCID](https://orcid.org/0009-0002-5327-1122))
