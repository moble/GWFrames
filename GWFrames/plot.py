"""
This submodule specializes the usual pyplot.plot function to be used
as methods with Waveform objects.

"""

def plotWaveform(this, WaveformPart='Abs', Modes=(), *pyplot_args, **pyplot_kwargs) :
    """
    This function should be called as a method of the Waveform class,
    e.g., as
    
    >>> W.plot('Abs', Modes=[[2,2], [2,-2]], c='r')
    
    where W is a Waveform object.  The first parameter should be a
    string --- one of ['Abs', 'LogAbs', 'LogLogAbs', 'Arg',
    'ArgUnwrapped', 'Re', Im'].  The second (optional) parameter is a
    list (using square brackets), where each element is some [l,m].
    Only modes included in this list will be plotted, unless the list
    is empty (the default), in which case all modes are plotted.
    
    All following arguments are passed to the usual pyplot.plot (or
    semilogx, semilogy, or loglog) function; in the example above, the
    "c='r'" argument is passed, making the line color red.  The only
    exception to this rule is the "label" keyword, which is handled
    internally when there is just one line.  Otherwise, the legend
    labels are set automatically to contain the [l,m] mode (though the
    legend is not shown by default).
    
    """
    
    from matplotlib.pyplot import plot, semilogx, semilogy, loglog, xlabel, ylabel, setp, isinteractive, ioff, ion, draw, show, gcf
    try : # If this matplotlib is too old, just ignore this
        from matplotlib.pyplot import tight_layout
    except :
        pass
    import matplotlib
    try:
        matplotlib.rcParams['axes.color_cycle'] = ['#000000', '#cc79a7', '#d55e00', '#0072b2', '#f0e442', '#56b4e9', '#e69f00', '#2b9f78']
    except KeyError :
        matplotlib.axes.set_default_color_cycle(['#000000', '#cc79a7', '#d55e00', '#0072b2', '#f0e442', '#56b4e9', '#e69f00', '#2b9f78'])
    from warnings import warn
    from numpy import array, empty, transpose, sin, cos
    
    def NotAbs(Anything) :
        return Anything
    
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
        quantity = this.ArgUnwrapped
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
        Labels = [str(tuple(this.LM(Modes)))]
    elif (len(Modes)==0) : # All modes are plotted
        Lines = styledplot(this.T(), AbsOrNot(quantity()).transpose(), *pyplot_args, **pyplot_kwargs)
        if ((this.NModes()==1) and ('label' in pyplot_kwargs)) :
            Labels = [pyplot_kwargs['label']]
        else :
            Labels = [str(tuple(this.LM(mode))) for mode in range(this.NModes())]
    elif ( (len(Modes)==2) and (type(Modes[0])==int and type(Modes[1])==int) ) :
        ModeIndex = this.FindModeIndex(Modes[0], Modes[1])
        Lines = styledplot(this.T(), transpose(AbsOrNot(quantity(ModeIndex))), *pyplot_args, **pyplot_kwargs)
        if ('label' in pyplot_kwargs) :
            Labels = [pyplot_kwargs['label']]
        else :
            Labels = [str(tuple(this.LM(ModeIndex)))]
    elif ( (len(Modes)==1) and (type(Modes[0])==list) and (type(Modes[0][0])==int and type(Modes[0][1])==int) ) :
        ModeIndex = this.FindModeIndex(Modes[0][0], Modes[0][1])
        Lines = styledplot(this.T(), transpose(AbsOrNot(quantity(ModeIndex))), *pyplot_args, **pyplot_kwargs)
        if ('label' in pyplot_kwargs) :
            Labels = [pyplot_kwargs['label']]
        else :
            Labels = [str(tuple(this.LM(ModeIndex)))]
    else :
        Modes = array(Modes, dtype=int)
        for i in range(Modes.shape[0]) :
            ModeIndex = this.FindModeIndex(int(Modes[i][0]), int(Modes[i][1]))
            Labels.append(str(tuple(Modes[i])))
            Lines.append(styledplot(this.T(), AbsOrNot(quantity(ModeIndex)).transpose(), *pyplot_args, **pyplot_kwargs))
    
    xlabel(XLabel)
    ylabel(YLabel)
    
    draw()
    
    if(WasInteractive) :
        ion()
        gcf().show()
    
    try :
        tight_layout(pad=0.1)
    except :
        pass
    
    for i in range(len(Lines)) :
        setp(Lines[i], label=Labels[i])
    
    return Lines



### The following allows us to write things like
###   W = GWFrames.Waveform()
###   W.plot('Abs', ((2,2), (2,-2)))
import GWFrames
GWFrames.Waveform.plot = plotWaveform
