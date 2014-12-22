#!/usr/bin/env python
# encoding: utf-8
"""
Plots tide based sea level. We will maintain a separate GeoJSON file
containing station names and positions.
"""
# Standard library imports
import datetime
import cStringIO
import glob
import os.path
import collections
import math

# External imports
import cherrypy
import numpy as np
import pytz
import matplotlib as mpl
mpl.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# Setup Matplotlib parameters
mpl.rcParams['font.size'] = 8
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['axes.titlesize'] = 8

FIGSIZE = 6 # Inches
STD_DPI = 100
STD_XSIZE = 600

def nice_axis_limits(v):
    """Sets limits according to data."""
    vmin = np.amin(v)
    vmax = np.amax(v)
    exponent, remainder = divmod(math.log10(vmax - vmin), 1)
    if remainder < 0.5:
        exponent -= 1
    scale = 10**(-exponent)
    vmin = math.floor(scale*vmin)/scale
    vmax = math.ceil(scale*vmax)/scale
    vmin, vmax = mpl.transforms.nonsingular(vmin, vmax, expander=0.5)
    return vmin, vmax

def plot_timeseries(tseries, plotopts={}):
    """\
    Plots a time series.

    Arguments:
    plotbuf -- output png buffer
    tseries -- timeseries dict
    plotopts -- Plot options (not yet used)
    """
    timezone = plotopts.get('tz', pytz.timezone('UTC'))
    xsize = float(plotopts.get('xsize', STD_XSIZE))
    std_ysize = int(round(0.85*xsize))
    ysize = float(plotopts.get('ysize', std_ysize))
    autotitle = plotopts.get('autotitle', False)
    vmin_range = plotopts.get('vmin_range', None)

    plotbuf = cStringIO.StringIO()
    data = tseries['values']
    dates = [d.astimezone(timezone) for d in tseries['datetimes']]
    mdates = mpl.dates.date2num(dates)

    xfig = xsize/STD_DPI
    yfig = ysize/STD_DPI
    fig = mpl.figure.Figure(figsize=(xfig, yfig), dpi=STD_DPI)
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(1, 1, 1)
    ax.grid(True)
    plines1 = ax.plot_date(mdates, data, 'b-', tz=timezone)

    # Set y axis limits
    vmin, vmax = nice_axis_limits(data)
    if vmin_range is not None:
        dv = vmin_range - (vmax - vmin)
        if dv > 0:
            vmax += 0.5*vmin_range
            vmin -= 0.5*vmin_range
    ax.set_ylim(vmin, vmax)

    # Set x axis limits
    vmin = mdates[0]
    vmax = mdates[-1]
    ax.set_xlim(vmin, vmax)

    ax.xaxis.set_major_formatter(mpl.ticker.NullFormatter())
    fltfmt = mpl.ticker.FormatStrFormatter('%g')
    ax.yaxis.set_major_formatter(fltfmt)
    ax.set_ylabel(tseries['long_name'] + ' (' + tseries['units'] \
                       + ')')

    # Time ticks
    minval = int(round(24.*np.min(mdates)-0.5)+.001)
    maxval = int(round(24.*np.max(mdates)+0.5)+.001)
    nticks = 6
    dt = 3
    interv = max(int(round((maxval-minval)/nticks/dt)*dt), dt)

    majloc = mpl.dates.HourLocator(interval=interv)
    minloc = mpl.dates.HourLocator(interval=interv/3)
    ax.xaxis.set_major_locator(majloc)
    ax.xaxis.set_minor_locator(minloc)
    labels = ax.get_xticklabels()

    for label in labels:
        label.set_rotation(45)
        label.set_fontsize(6)
    timfmt = mpl.dates.DateFormatter('%Y-%m-%d\n%H:%M ' + dates[0].tzname(), tz=timezone)
    ax.xaxis.set_major_formatter(timfmt)

    canvas = FigureCanvas(fig)
    # the tight_layout method modifies the figure every time it is called
    fig.tight_layout()
    fig.canvas.print_png(plotbuf)
    plotbuf = plotbuf.getvalue()
    return str(plotbuf)

class FixedOffset(datetime.tzinfo):
    """Fixed offset in minutes east from UTC."""

    def __init__(self, offset, name):
        self.__offset = datetime.timedelta(minutes=offset)
        self.__name = name

    def utcoffset(self, dt):
        return self.__offset

    def tzname(self, dt):
        return self.__name

    def dst(self, dt):
        return datetime.timedelta(0)

def _(text, lang):
    """Poor mans translation method."""
    da_dict = {'Sea level above lowest astronomical tide': \
               'Vandstand over laveste astronomiske tidevand'}
    langdict = {'da': da_dict}
    if lang in langdict:
        if text in langdict[lang]:
            result = langdict[lang][text]
        else:
            result = text
    else:
        result = text
    return result

def validate_timezone(tz):
    """Validates timezone input and returns a fixed offset timezone."""
    tz = int(tz)
    if tz == 0:
        tzinfo = pytz.UTC
    else:
        tzh = int(float(tz)/60)
        if tz < 0:
            sign = '-'
        else:
            sign = '+'
        tzm = abs(tz - tzh*60)
        tzname = "UTC\n%s%02d:%02d" % (sign, abs(tzh), tzm)
        tzinfo = FixedOffset(tz, tzname)
    return tzinfo

def validate_time(start, end):
    """\
    Validates start and end time input. Returns UTC datetime objects
    corresponding to input.
    """
    timefmt = '%Y-%m-%dT%H:%M:%S'
    if start is not None:
        try:
            start = datetime.datetime.strptime(start[0:19], '%Y-%m-%dT%H:%M:%S')
            start = start.replace(tzinfo=pytz.UTC)
        except ValueError as err:
            msg = 'Invalid time specification: %s' % start
            raise cherrypy.HTTPError("400 Bad Request", msg)
    if end is not None:
        try:
            end = datetime.datetime.strptime(end[0:19], '%Y-%m-%dT%H:%M:%S')
            end = end.replace(tzinfo=pytz.UTC)
        except ValueError as err:
            msg = 'Invalid time specification: %s' % start
            raise cherrypy.HTTPError("400 Bad Request", msg)
    if start is not None and end is not None:
        if start >= end:
            msg = 'Invalid time range: %s, %s' % (start, end)
            raise cherrypy.HTTPError("404 Not Found", msg)
    return start, end

def mimetype(mtype):
    def decorate(func):
        def wrapper(*args, **kwargs):
            cherrypy.response.headers['Content-Type'] = mtype
            return func(*args, **kwargs)
        return wrapper
    return decorate

class TidalData(object):
    """\
    Class which holds the tidal data.
    """
    def __init__(self, datadir):
        """Initializes (loads) tidal data."""
        self.datadir = datadir
        self.stations = self.load_stations()

    def get_station(self, station, start, end, lang):
        """Returns station data in a dictionary."""
        t, v = self.stations[station]
        istart = np.searchsorted(t, start, side='left')
        iend = np.searchsorted(t, end, side='right')
        t = t[istart:iend]
        v = v[istart:iend]
        long_name = _('Sea level above lowest astronomical tide', lang)
        units = _('m', lang)
        tseries = {'datetimes': t,
                   'values': v,
                   'long_name': long_name,
                   'units': units}
        return tseries

    def load_stations(self, start=None, end=None):
        """Loads sea level data for all stations."""
        datafiles = glob.glob(self.datadir + '/*.txt')
        datafiles.sort()
        basenames = map(os.path.basename, datafiles)
        names = [bn.split('20')[0] for bn in basenames]
        stations = {}
        for name in names:
            if not name.isalpha():
                msg = 'Invalid station name: %s' % name
                raise cherrypy.HTTPError("403 Forbidden", msg)
            filetemplate = self.datadir + '/%(name)s2015_16_10min.txt'
            filename = filetemplate % {'name': name}
            npzfile = filename + '.npz'
            if os.path.isfile(npzfile):
                print('Loading %s' % npzfile)
                t, v = self.load_station_npz(npzfile)
            else:
                print('Loading %s' % filename)
                t, v = self.load_station_ascii(filename)
                # Store efficient version for subsequent server restarts
                print('Saving %s' % npzfile)
                np.savez(npzfile, t=t, v=v)
            stations[name] = t, v
        return stations

    def load_station_npz(self, filename):
        """Loads sea level data for a given station."""
        npzfile = np.load(filename)
        t = npzfile['t']
        v = npzfile['v']
        return t, v

    def load_station_ascii(self, filename, start=None, end=None):
        """Loads sea level data for a given station."""
        t = []
        v = []
        offset = datetime.timedelta(hours=3)
        with open(filename, 'r') as infile:
            # Skip header
            infile.next()
            for line in infile:
                #print line
                ldate, ltime, ldata = line.split()
                #ldatetime = datetime.datetime.strptime(ldate + ltime, '%d/%m/%Y%H:%M')
                dtargs = map(int, [ldate[6:10], ldate[3:5], ldate[:2],
                         ltime[:2], ltime[3:5]])
                ldatetime = datetime.datetime(*dtargs, tzinfo=pytz.UTC)
                #ldatetime = datetime.datetime.strptime(ldate + ltime, '%d/%m/%Y%H:%M')
                ldatetime += offset
                #ldatetime = ldatetime.replace(tzinfo=pytz.UTC)
                if start is None or start <= ldatetime:
                    if end is None or end >= ldatetime:
                        t.append(ldatetime)
                        v.append(float(ldata))
        t = np.array(t)
        v = np.array(v)
        return t, v

class TidePlotter(object):
    def __init__(self):
        self.tidal_data = TidalData(datadir='data')

    @cherrypy.expose
    @mimetype("image/png")
    def index(self, station, start=None, end=None, nx=600, ny=400, tz=0, lang='en'):
        """Returns sea level plot."""
        # Input validation and parsing
        if not station.isalpha():
            msg = 'Invalid station name: %s' % station
            raise cherrypy.HTTPError("403 Forbidden", msg)
        nx = min(1024, max(128, int(nx)))
        ny = min(1024, max(128, int(ny)))
        start, end = validate_time(start, end)
        tzinfo = validate_timezone(tz)
        # Load time series
        try:
            tseries = self.tidal_data.get_station(station, start=start, 
                                                  end=end, lang=lang)
        except:
            raise
            msg = 'No data present for %s' % station
            raise cherrypy.HTTPError("404 Not Found", msg)
        # Plot time series
        plotopts = {'xsize': nx, 'ysize': ny, 'tz': tzinfo}
        png = plot_timeseries(tseries, plotopts)
        return png

def main():
    cherrypy.config.update({'server.socket_port': 7070,
                            'server.socket_host': '0.0.0.0',
                            'environment': 'production'})
    cherrypy.quickstart(TidePlotter())

def test():
    tidal_data = TidalData(datadir='data')

if __name__ == '__main__':
    #import cProfile
    #cProfile.run('test()', 'tidal.prof')
    exit(main())
