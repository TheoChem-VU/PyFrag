try:
    from functools import lru_cache
except ImportError:
    # Python 2 does stdlib does not have lru_cache so let's just
    # create a dummy decorator to avoid crashing
    print ("WARNING: Cache for this example is available on Python 3 only.")
    def lru_cache():
        def dec(f):
            def _(*args, **kws):
                return f(*args, **kws)
            return _
        return dec
import os
from os.path import dirname, join
import pandas as pd
from bokeh.io import curdoc
from bokeh.layouts import row, column
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import PreText, Select
from bokeh.plotting import figure
import numpy as np
from bokeh.layouts import gridplot
from bokeh.plotting import figure, show, output_file




DATA_DIR = join(dirname(__file__), 'daily')


DEFAULT_TICKERS = ['p', 'r1', 'r2', 'rc', 'ts']


@lru_cache()
def load_ticker(ticker):
    fname = join(DATA_DIR, '%s.csv' % ticker.lower())
    data = pd.read_csv(fname, header=None, parse_dates=['step'],
                       names=['step', 'energy', 'energy_change', 'energy_change_', 'e_TF', 'constrained_gradient_max', 'constrained_gradient_max_', 'gm_TF', 'constrained_gradient_rms', 'constrained_gradient_rms_', 'gr_TF', 'cart_step_max', 'cart_step_max_', 'sm_TF', 'cart_step_rms', 'cart_step_rms_', 'sr_TF'])
    data = data.set_index('step')
    return pd.DataFrame({ticker: data.energy, 'E': data.energy_change, 'E_': data.energy_change_, 'e_TF': data.e_TF, 'c_max': data.constrained_gradient_max, 'c_max_': data.constrained_gradient_max_, 'gm_TF': data.gm_TF, 'c_rms': data.constrained_gradient_rms, 'c_rms_': data.constrained_gradient_rms_, 'gr_TF': data.gr_TF, 's_max': data.cart_step_max, 's_max_': data.cart_step_max_, 'sm_TF': data.sm_TF, 's_rms': data.cart_step_rms, 's_rms_': data.cart_step_rms_, 'sr_TF': data.sr_TF})

@lru_cache()
def get_data(t1):
    df1 = load_ticker(t1)
    data = pd.concat([df1], axis=1)
    data = data.dropna()
    data['t1']               = data[t1]
    data['t1_E']             = data['E']
    data['t1_E_']            = data['E_']
    data['t1_e_TF']          = data['e_TF']
    data['t1_c_max']         = data['c_max']
    data['t1_c_max_']        = data['c_max_']
    data['t1_gm_TF']         = data['gm_TF']
    data['t1_c_rms']         = data['c_rms']
    data['t1_c_rms_']        = data['c_rms_']
    data['t1_gr_TF']         = data['gr_TF']
    data['t1_s_max']         = data['s_max']
    data['t1_s_max_']        = data['s_max_']
    data['t1_sm_TF']         = data['sm_TF']
    data['t1_s_rms']         = data['s_rms']
    data['t1_s_rms_']        = data['s_rms_']
    data['t1_sr_TF']         = data['sr_TF']
    return data

# set up callbacks
def update(selected=None):
    t1 = ticker1.value
    data = get_data(t1)
    source.data = source.from_df(data[['t1',  't1_E', 't1_E_',  't1_e_TF', 't1_c_max',  't1_c_max_', 't1_gm_TF',  't1_c_rms', 't1_c_rms_',  't1_gr_TF', 't1_s_max',  't1_s_max_', 't1_sm_TF',  't1_s_rms', 't1_s_rms_', 't1_sr_TF']])
    source_static.data = source.data
    update_stats(data, t1)
    ts1.title.text = t1

def update_stats(data, t1):
    print (data, t1)
    stats.text = str(data[[t1, 'E', 'E_', 'e_TF', 'c_max', 'c_max_', 'gm_TF', 'c_rms', 'c_rms_', 'gr_TF', 's_max', 's_max_', 'sm_TF', 's_rms', 's_rms_', 'sr_TF']])


def ticker1_change(attrname, old, new):
    ticker1.options = DEFAULT_TICKERS
    update()

def selection_change(attrname, old, new):
    t1 = ticker1.value
    data = get_data(t1)
    selected = source.selected.indices
    if selected:
        data = data.iloc[selected, :]
    update_stats(data, t1)

# set up widgets
stats = PreText(text='', width=800)
ticker1 = Select(value='p', options=DEFAULT_TICKERS)

#set up plot

source        = ColumnDataSource(data=dict(step=[], t1=[], t1_E=[], t1_E_=[],  t1_e_TF=[], t1_c_max=[],  t1_c_max_=[], t1_gm_TF=[],  t1_c_rms=[], t1_c_rms_=[],  t1_gr_TF=[], t1_s_max=[],  t1_s_max_=[], t1_sm_TF=[],  t1_s_rms=[], t1_s_rms_=[], t1_sr_TF=[]))
source_static = ColumnDataSource(data=dict(step=[], t1=[], t1_E=[], t1_E_=[],  t1_e_TF=[], t1_c_max=[],  t1_c_max_=[], t1_gm_TF=[],  t1_c_rms=[], t1_c_rms_=[],  t1_gr_TF=[], t1_s_max=[],  t1_s_max_=[], t1_sm_TF=[],  t1_s_rms=[], t1_s_rms_=[], t1_sr_TF=[]))
tools = 'pan,wheel_zoom,xbox_select,reset'

ts1 = figure(plot_width=900, plot_height=200, tools=tools, active_drag="xbox_select")
ts1.line('step', 't1', source=source_static)
ts1.circle('step', 't1', size=1, source=source, color=None, selection_color="orange")


ticker1.on_change('value', ticker1_change)
source.on_change('selected', selection_change)


###########pyfrag figure part############

widgets = column(ticker1, stats)
main_row = row(widgets)
series = column(ts1)

current_path=os.path.abspath('.')
pyfrag_sign=os.path.exists(os.path.join(current_path,'stocks/PYFRAG.csv'))

if pyfrag_sign==False:
    layout = column(main_row, series)
else:
    from pyfrag import PYFRAG
    # Figure one for ASM which plot energy of total energy, interaction energy and total strain energy
    # againt bondlength change. You can add more data plot by adding a line like:
    # p1.line(PYFRAG['bondlength'], PYFRAG['something'], color='#A6CEE3', legend='something')
    p1 = figure(title="Activation strain analysis")
    p1.grid.grid_line_alpha=0.3
    p1.xaxis.axis_label = 'bondlength'
    p1.yaxis.axis_label = 'energy'

    p1.line(PYFRAG['bondlength'], PYFRAG['energytotal'], color='#A6CEE3', legend='energytotal')
    p1.line(PYFRAG['bondlength'], PYFRAG['inter'],       color='#B2DF8A', legend='inter')
    p1.line(PYFRAG['bondlength'], PYFRAG['straintotal'], color='#33A02C', legend='straintotal')
    p1.legend.location = "top_left"

    # Figure two for ASM which plot the decomposition of interaction energy into
    # orbital, pauli, electronic statistic energy terms.
    p2 = figure(title="Activation strain analysis")
    p2.grid.grid_line_alpha=0.3
    p2.xaxis.axis_label = 'bondlength'
    p2.yaxis.axis_label = 'energy'

    p2.line(PYFRAG['bondlength'], PYFRAG['inter'], color='#A6CEE3', legend='inter')
    p2.line(PYFRAG['bondlength'], PYFRAG['oi'], color='#B2DF8A', legend='oi')
    p2.line(PYFRAG['bondlength'], PYFRAG['pauli'], color='#33A02C', legend='pauli')
    p2.line(PYFRAG['bondlength'], PYFRAG['elstat'], color='#FB9A99', legend='elstat')
    p2.legend.location = "top_left"

    # Lay out the position of two figures and show figures.
    # You can add more figures simply by adding something like:
    # p3 = figure(title="Activation strain analysis")
    # p3.grid.grid_line_alpha=0.3
    # p3.xaxis.axis_label = 'energy'
    # p3.yaxis.axis_label = 'bondlength'
    # p3.line(PYFRAG['bondlength'], PYFRAG['inter'], color='#A6CEE3', legend='inter')
    # p3.line(PYFRAG['bondlength'], PYFRAG['inter'], color='#A6CEE3', legend='inter')
    # p3.legend.location = "top_left"
    # layout = column(main_row, series,gridplot([[p1,p2,p3]], plot_width=400, plot_height=400))

    layout = column(main_row, series,gridplot([[p1,p2]], plot_width=400, plot_height=400))

# initialize
update()
curdoc().add_root(layout)
curdoc().title = "Stocks"
