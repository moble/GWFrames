"""
This submodule specializes the usual pyplot.plot function to be used
as methods with Waveform objects.

"""


def _ModeLabel(LM):
    """Return a label like '(2, -2)', even when LM holds numpy integers"""
    return str(tuple(int(x) for x in LM))


def plotWaveform(this, WaveformPart='Abs', Modes=(), t_fid=None, *pyplot_args, **pyplot_kwargs) :
    """Plot some part of the modes of a Waveform against time.

    Importing `GWFrames.plot` attaches this function to `Waveform` and
    `PNWaveform` as the method ``plot``, so that it is called as, for example,

    >>> import GWFrames.plot
    >>> W.plot('Abs', Modes=[[2,2], [2,-2]], c='r')

    which plots the magnitudes of the (2, 2) and (2, -2) modes of ``W`` in red.

    Parameters
    ----------
    WaveformPart : str, default 'Abs'
        The part of the data to plot: one of ``'Abs'``, ``'LogAbs'``,
        ``'LogLogAbs'``, ``'Arg'``, ``'ArgUnwrapped'``, ``'Re'``, or ``'Im'``.
        The ``Log`` forms use logarithmic axes.
    Modes : int or sequence, optional
        The modes to plot, given as a single mode index, as one ``[ell, m]``
        pair, or as a sequence of such pairs.  If empty (the default), all
        modes are plotted.
    t_fid : float, optional
        For ``'ArgUnwrapped'``, the phases are shifted by multiples of
        :math:`2\\pi` so that they lie between 0 and :math:`2\\pi` at the
        first time step after ``t_fid``.
    *pyplot_args, **pyplot_kwargs
        Passed on to `matplotlib.pyplot.plot` (or ``semilogy`` or
        ``loglog``).  The legend label of each line is set to its
        :math:`(\\ell, m)` automatically, unless a single line is plotted and a
        ``label`` is given.

    Returns
    -------
    list of matplotlib.lines.Line2D
        The lines plotted.
    """

    from matplotlib.pyplot import plot, semilogx, semilogy, loglog, xlabel, ylabel, setp, isinteractive, ioff, ion, draw, show, gcf
    try : # If this matplotlib is too old, just ignore this
        from matplotlib.pyplot import tight_layout
    except :
        pass
    import matplotlib
    # # We probably shouldn't be overwriting people's color choices
    # try:
    #     matplotlib.rcParams['axes.color_cycle'] = ['#000000', '#cc79a7', '#d55e00', '#0072b2', '#f0e442', '#56b4e9', '#e69f00', '#2b9f78']
    # except KeyError :
    #     matplotlib.axes.set_default_color_cycle(['#000000', '#cc79a7', '#d55e00', '#0072b2', '#f0e442', '#56b4e9', '#e69f00', '#2b9f78'])
    from warnings import warn
    from numpy import array, empty, transpose, sin, cos, modf, pi

    def NotAbs(Anything) :
        return Anything

    def _ArgUnwrapped(Mode=None) :
        if (t_fid) :
            i_fid = int(array(range(this.NTimes()))[this.T()>t_fid][0])
            if (Mode is None) :
                Arg = this.ArgUnwrapped()
                x = modf(Arg[:,i_fid]/(2*pi))
                x[1][x[0]<0.0] -= 1.0
                return (Arg.transpose() - x[1]*2*pi).transpose()
            else :
                Arg = this.ArgUnwrapped(Mode)
                w,x = modf(Arg[i_fid]/(2*pi))
                if (w<0.0) : x -= 1.0
                return Arg - x*2*pi
        else :
            if (Mode is None) :
                return this.ArgUnwrapped()
            else :
                return this.ArgUnwrapped(Mode)

    XLabel = r'$(t-r_\ast)/M$'
    YLabel = ''

    Labels = []
    Lines = []

    WasInteractive = isinteractive()
    ioff()

    ## This decides which quantity to plot, what type of plot, and what the labels should be
    if (WaveformPart.lower()=='abs') :
        YLabel =r'$\mathrm{abs} \left( '+this.GetLaTeXDataDescription()+r' \right) $'
        styledplot = plot
        quantity = this.Abs
        AbsOrNot = NotAbs
    elif (WaveformPart.lower()=='logabs') :
        YLabel =r'$\mathrm{abs} \left( '+this.GetLaTeXDataDescription()+r' \right) $'
        styledplot = semilogy
        quantity = this.Abs
        AbsOrNot = abs
    elif (WaveformPart.lower()=='loglogabs') :
        YLabel =r'$\mathrm{abs} \left( '+this.GetLaTeXDataDescription()+r' \right) $'
        styledplot = loglog
        quantity = this.Abs
        AbsOrNot = abs
    elif (WaveformPart.lower()=='arg') :
        YLabel =r'$\mathrm{arg} \left( '+this.GetLaTeXDataDescription()+r' \right) $'
        styledplot = plot
        quantity = this.Arg
        AbsOrNot = NotAbs
    elif (WaveformPart.lower()=='argunwrapped' or WaveformPart.lower()=='uarg' or WaveformPart.lower()=='argu') :
        YLabel =r'$\mathrm{uarg} \left( '+this.GetLaTeXDataDescription()+r' \right) $'
        styledplot = plot
        quantity = _ArgUnwrapped
        AbsOrNot = NotAbs
    elif (WaveformPart.lower()=='real' or WaveformPart.lower()=='re') :
        YLabel =r'$\mathrm{Re} \left( '+this.GetLaTeXDataDescription()+r' \right) $'
        styledplot = plot
        quantity = this.Re
        AbsOrNot = NotAbs
    elif (WaveformPart.lower()=='imaginary' or WaveformPart.lower()=='imag' or WaveformPart.lower()=='im') :
        YLabel =r'$\mathrm{Im} \left( '+this.GetLaTeXDataDescription()+r' \right) $'
        styledplot = plot
        quantity = this.Im
        AbsOrNot = NotAbs
    else :
        print("Unsupported plot type '{0}'".format(WaveformPart))
        return []

    ## This does the actual work of plotting, depending on what Modes are needed
    if (type(Modes)==int) : # The requested mode is plotted
        Lines = styledplot(this.T(), AbsOrNot(quantity(Modes)), *pyplot_args, **pyplot_kwargs)
        Labels = [_ModeLabel(this.LM(Modes))]
    elif (len(Modes)==0) : # All modes are plotted
        Lines = styledplot(this.T(), AbsOrNot(quantity()).transpose(), *pyplot_args, **pyplot_kwargs)
        if ((this.NModes()==1) and ('label' in pyplot_kwargs)) :
            Labels = [pyplot_kwargs['label']]
        else :
            Labels = [_ModeLabel(this.LM(mode)) for mode in range(this.NModes())]
    elif ( (len(Modes)==2) and (type(Modes[0])==int and type(Modes[1])==int) ) :
        ModeIndex = this.FindModeIndex(Modes[0], Modes[1])
        Lines = styledplot(this.T(), transpose(AbsOrNot(quantity(ModeIndex))), *pyplot_args, **pyplot_kwargs)
        if ('label' in pyplot_kwargs) :
            Labels = [pyplot_kwargs['label']]
        else :
            Labels = [_ModeLabel(this.LM(ModeIndex))]
    elif ( (len(Modes)==1) and (type(Modes[0])==list) and (type(Modes[0][0])==int and type(Modes[0][1])==int) ) :
        ModeIndex = this.FindModeIndex(Modes[0][0], Modes[0][1])
        Lines = styledplot(this.T(), transpose(AbsOrNot(quantity(ModeIndex))), *pyplot_args, **pyplot_kwargs)
        if ('label' in pyplot_kwargs) :
            Labels = [pyplot_kwargs['label']]
        else :
            Labels = [_ModeLabel(this.LM(ModeIndex))]
    else :
        Modes = array(Modes, dtype=int)
        for i in range(Modes.shape[0]) :
            ModeIndex = this.FindModeIndex(int(Modes[i][0]), int(Modes[i][1]))
            Labels.append(_ModeLabel(Modes[i]))
            Lines.append(styledplot(this.T(), AbsOrNot(quantity(ModeIndex)).transpose(), *pyplot_args, **pyplot_kwargs))

    xlabel(XLabel)
    ylabel(YLabel)

    # draw()

    if(WasInteractive) :
        ion()
        # gcf().show()

    # try :
    #     tight_layout(pad=0.1)
    # except :
    #     pass

    for i in range(len(Lines)) :
        setp(Lines[i], label=Labels[i])

    return Lines



### The following allows us to write things like
###   W = GWFrames.Waveform()
###   W.plot('Abs', ((2,2), (2,-2)))
import GWFrames
GWFrames.Waveform.plot = plotWaveform
GWFrames.PNWaveform.plot = plotWaveform
