import numpy as np
from math import ceil
from matplotlib import gridspec
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

from calc.core_functionality import calc_broad

# Graph setup
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage[cm]{sfmath} \usepackage{amsmath}'
mpl.rcParams['font.size'] = 18
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 2.5

class BarsSpectrum():
    """
    Creates bars to be plotted
    Usage:
     BarsSpectrum(start,end,numpts,peaks,width)
    where
     peaks -- ( [List of peaks],[list of heights] )
    """
    def __init__(self, start, end, numpts, peaks, width):
        self.start = start
        self.end = end
        self.numpts = numpts
        self.peaks = peaks[0]
        self.heights = peaks[1]
        self.width = width


        # make heights a local variable as it's accessed in the inner loop
        heights = self.heights

        self.xvalues = np.arange(self.numpts)*float(self.end-self.start)/(self.numpts-1) + self.start
        step_ind = ceil(width/(self.xvalues[1] - self.xvalues[0]))
        indices = np.arange(-step_ind, step_ind+1, 1)
        spec = np.zeros(self.xvalues.shape)
        indix = np.abs(self.peaks[:, np.newaxis] - self.xvalues[np.newaxis, :]).argmin(axis=1)
        for ith in range(heights.shape[0]):
            pos = indices + indix[ith]
            pos = pos[np.logical_or(pos >= 0, pos < self.xvalues.shape[0])].tolist()
            spec[pos] += heights[ith]

        self.spectrum = spec


def plot_colorbar(data, start, end, hwhm=10, save=False):
    imgext = 'pdf'
    fig_geom = (12, 4)
    smin = start
    smax = end
    hwhmo = hwhm
    numpts = int((smax-smin)*5)
    xvalues = np.arange(numpts)*float(end - start)/(numpts-1) + start
    #fig = plt.figure(dpi=300, figsize=fig_geom) #height_ratios=fig_geom)
    nfig = len(data)
    ratios = [1 for _ in range(nfig)]
    ratios[0] = 8
    figid = 1
    fig, axes = plt.subplots(nfig, 1, dpi=300, num=figid, figsize=fig_geom,  # sharex=True,
                             gridspec_kw={'height_ratios': ratios})
    # gs = gridspec.GridSpec(nrows=nfig, ncols=1, hspace=0, wspace=0.01, height_ratios=ratios)
    extent = smin, smax, -10, 10
    result = []
    # axes = []
    for i, subdata in enumerate(data):
        values = calc_broad(subdata, smin, smax, hwhmo)
        result.append(values.spectrum[0, :])
        # axes.append(fig.add_subplot(gs[i, 0]))
        if not i:
            axes[i].plot(xvalues, result[-1])
            axes[i].set_xlim(start, end)
            axes[i].set_xticklabels([])
            axes[i].set_xticks([])
            # axes[-1].set_facecolor('black')
        else:
            axes[i].imshow(result[-1].reshape((1, -1)), extent=extent,
                           cmap='bwr_r', vmax=1, vmin=-1,
                           aspect='auto')
            axes[i].set_yticklabels([])
            axes[i].set_yticks([])
            if not i == len(data)-1:
                axes[i].set_xticklabels([])
                axes[i].set_xticks([])
            axes[i].set_ylabel('Frg: {}'.format(i), rotation=0, labelpad=20)
            # axes[-1].set_facecolor('black')

    axes[-1].set_xlabel(r'Frequency ($cm^{-1}$)')
    fig.subplots_adjust(hspace=0.05)
    if save:
        fmt_o = 'ecd_spect.'+imgext
        outname = fmt_o
        plt.savefig(outname, num=figid, bbox_inches='tight')
        plt.close(figid)
    else:
        plt.show()

def plot_colorbar2(data, start=None, end=None):
    imgext = 'pdf'
    fig_geom = (12, 4)
    if not start:
        smin = data.xvals.min()
    else:
        smin = start
    if not end:
        smax = data.xvals.max()
    else:
        smax = end
    nfig = len(data.yvals)
    ratios = [1 for _ in range(nfig)]
    ratios[0] = 8
    fig, axes = plt.subplots(nfig,1, dpi=300, figsize=fig_geom, # sharex=True,
                             gridspec_kw = {'height_ratios':ratios})
    extent = smin, smax, -10, 10

    maxim = data.get_ymax()


    for i in range(nfig):
        if not i:
            axes[i].plot(data.xvals, data.yvals[i])
            axes[i].set_xlim(smin,smax)
            axes[i].set_xticklabels([])
            axes[i].set_xticks([])
            # axes[-1].set_facecolor('black')
        else:
            axes[i].imshow(data.yvals[i].reshape((1,-1)), extent=extent,
                           cmap='bwr_r', vmax=maxim, vmin=-maxim,
                           aspect='auto')
            axes[i].set_yticklabels([])
            axes[i].set_yticks([])
            if not i == nfig-1:
                axes[i].set_xticklabels([])
                axes[i].set_xticks([])
            axes[i].set_ylabel('Frg: {}'.format(i), rotation=0, labelpad=20)

    axes[-1].set_xlabel(r'Frequency ($cm^{-1}$)')

    fig.subplots_adjust(hspace=0.05)
    plt.show()



def plot_colorbar2(data, start=None, end=None):
    fig_geom = (12, 4)
    if not start:
        smin = data.xvals.min()
    else:
        smin = start
    if not end:
        smax = data.xvals.max()
    else:
        smax = end
    nfig = len(data.yvals)
    ratios = [1 for _ in range(nfig)]
    ratios[0] = 8
    fig, axes = plt.subplots(nfig,1, dpi=300, figsize=fig_geom, # sharex=True,
                             gridspec_kw = {'height_ratios':ratios})
    extent = smin, smax, -10, 10

    maxim = data.get_ymax()

    for i in range(nfig):
        if not i:
            axes[i].plot(data.xvals, data.yvals[i])
            axes[i].set_xlim(smax,smin)
            axes[i].set_xticklabels([])
            axes[i].set_xticks([])
            # axes[-1].set_facecolor('black')
        else:
            axes[i].imshow(data.yvals[i].reshape((1,-1)), extent=extent,
                           cmap='bwr_r', vmax=maxim, vmin=-maxim,
                           aspect='auto')
            axes[i].set_yticklabels([])
            axes[i].set_yticks([])
            if not i == nfig-1:
                axes[i].set_xticklabels([])
                axes[i].set_xticks([])
            axes[i].set_ylabel('Frg: {}'.format(i), rotation=0, labelpad=20)

    axes[-1].set_xlabel(r'Wavenumber ($cm^{-1}$)')

    fig.subplots_adjust(hspace=0.05)
    return fig

def align_yaxis(ax1, ax2):
    y_lims = np.array([ax.get_ylim() for ax in [ax1, ax2]])

    # force 0 to appear on both axes, comment if don't need
    y_lims[:, 0] = y_lims[:, 0].clip(None, 0)
    y_lims[:, 1] = y_lims[:, 1].clip(0, None)

    # normalize both axes
    y_mags = (y_lims[:,1] - y_lims[:,0]).reshape(len(y_lims),1)
    y_lims_normalized = y_lims / y_mags

    # find combined range
    y_new_lims_normalized = np.array([np.min(y_lims_normalized), np.max(y_lims_normalized)])

    # denormalize combined range to get new axes
    new_lim1, new_lim2 = y_new_lims_normalized * y_mags
    ax1.set_ylim(new_lim1)
    ax2.set_ylim(new_lim2)

def add_second_axis(axes_org, xvals, yvals):
    new_axes = axes_org.twinx()
    new_axes.vlines(xvals, 0, yvals, picker=5)
    align_yaxis(axes_org, new_axes)


def hextorgb(color_hex):
    hexi = color_hex[1:]
    return tuple(int(hexi[i:i+2], 16) for i in (0, 2, 4))

#see https://stackoverflow.com/a/3943023
def text_color(background):
    res = []
    for col in hextorgb(background):
        col = col / 255.0
        if col <= 0.03928:
            col = col/12.92
        else:
            col = ((col+0.055)/1.055)** 2.4
        res.append(col)
    lumin = 0.2126 * res[0] + 0.7152 * res[1] + 0.0722 * res[2]
    if lumin > np.sqrt(1.05 * 0.05) - 0.05:
        return '#000000'
    else:
        return '#ffffff'


def plot_colorbar3(fig, data, colors, start=None, end=None):
    if not start:
        smin = data.xvals.min()
    else:
        smin = start
    if not end:
        smax = data.xvals.max()
    else:
        smax = end
    nfig = len(data.yvals)
    ratios = [1 for _ in range(nfig)]
    ratios[0] = 12
    gs = fig.add_gridspec(nfig, 1, height_ratios=ratios)
    axes = []
    axes.append(fig.add_subplot(gs[0, 0]))
    for i in range(1, nfig):
        axes.append(fig.add_subplot(gs[i, 0], sharex=axes[0]))
    #axes = [fig.add_subplot(gs[x,0]) for x in range(nfig)]
    extent = smin, smax, -10, 10

    maxim = data.get_ymax()

    for i in range(nfig):
        if not i:
            axes[i].plot(data.xvals, data.yvals[i])
            axes[i].set_xlim(smax,smin)
            #axes[i].set_xticklabels([])
            axes[i].label_outer()
            # axes[-1].set_facecolor('black')
        else:
            axes[i].imshow(data.yvals[i].reshape((1,-1)), extent=extent,
                           cmap='bwr_r', vmax=maxim, vmin=-maxim,
                           aspect='auto')
            axes[i].set_yticklabels([])
            axes[i].set_yticks([])
            if not i == nfig-1:
                axes[i].label_outer()
                #axes[i].set_xticks([])
            textcolor = text_color(colors[i-1])
            axes[i].set_ylabel('Frg: {}'.format(i), color=textcolor,
                               backgroundcolor=colors[i-1],rotation=0,fontsize='small',
                               labelpad=30)
            #axes[i].yaxis.set_label_coords(-0.03,0.25)

    axes[-1].set_xlabel(r'Wavenumber ($cm^{-1}$)')

    fig.subplots_adjust(top=0.98, bottom=0.15, hspace=0.05)
    #fig.tight_layout()


def getcolorshexa(vector, cmap=cm.bwr_r):
    lim = max(vector.max(), np.abs(vector.min()))
    norm = mpl.colors.Normalize(vmin=-lim, vmax=lim)
    # Fix colormap usage for modern matplotlib
    rgba_values = cmap(norm(vector))
    collist = []
    for rgba in rgba_values:
        # Convert RGBA to hex, ignoring alpha channel
        r, g, b = rgba[:3]
        collist.append("#{0:02x}{1:02x}{2:02x}".format(
            int(r * 255), int(g * 255), int(b * 255)
        ))
    return collist
