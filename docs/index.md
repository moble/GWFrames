# GWFrames

GWFrames manipulates gravitational waveforms: it changes the frames in which
their modes are decomposed, aligns and compares them, hybridizes them, and
constructs post-Newtonian waveforms for precessing binaries.  The code
implements the methods described in the paper "[Angular velocity of
gravitational radiation from precessing binaries and the corotating
frame](https://arxiv.org/abs/1302.2919)," whose first commit it shares.

:::{note}
**For post-Newtonian waveforms, use
[PostNewtonian.jl](https://github.com/moble/PostNewtonian.jl).**  It is far
more complete, accurate, and flexible than the post-Newtonian code in this
package.  Most of the rest of the functionality here has moved to
[`scri`](https://github.com/moble/scri) and
[`sxs`](https://github.com/sxs-collaboration/sxs).  GWFrames is kept working
mainly as a reference implementation to compare other codes against.
:::

The central objects are the {py:class}`~GWFrames.Waveform` class, which holds
the spin-weighted spherical-harmonic modes of a waveform together with the
frame in which they are decomposed, and its subclass
{py:class}`~GWFrames.PNWaveform`, which computes a post-Newtonian waveform.
Quaternions, which represent the rotations between frames, are in the
{py:mod}`GWFrames.Quaternions` submodule.

```python
import GWFrames

# A post-Newtonian waveform for a precessing binary
W = GWFrames.PNWaveform("TaylorT4", 0.0, [0.0, 0.0, 0.9], [0.5, 0.0, 0.0], 0.02)

# Transform it to the corotating frame, in which it is nearly constant
W.TransformToCorotatingFrame()
```

The [tutorial](tutorial) works through the paper with examples like this one,
and the [API reference](api/index) documents every class and function.

```{toctree}
:hidden:
:maxdepth: 2

installation
tutorial
api/index
changes
```

## Citing GWFrames

The license for using this software is basically open (see the
[LICENSE](https://github.com/moble/GWFrames/blob/main/LICENSE) file), though
citations to [the original paper](https://arxiv.org/abs/1302.2919) are
appreciated, where relevant.  Also, if your work depends on features of this
module that have not been described in a separate publication of mine, I would
appreciate the opportunity to be a coauthor.

This code is not optimized.  In many cases, obvious optimizations were
rejected in favor of clearer or simpler code, or code that simply reflected the
discussion in the paper more directly.  The code is generally more than fast
enough for interactive use on individual waveforms, but rotations may be
painfully slow for large data sets.

## Contributions

The code is primarily written by Mike Boyle.  Dan Hemberger helped by porting
some older code from [Triton](https://github.com/moble/Triton) to perform
noise-weighted overlap calculations, and with numerous bug reports and helpful
suggestions.  Serguei Ossokine also helped substantially by cross-checking the
post-Newtonian formulas and results.  Dante Iozzo helped to keep things working
as this code entered its dotage.  Claude (Opus 5.5) resuscitated the package for
Python 3.
