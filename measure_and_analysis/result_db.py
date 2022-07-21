#!/usr/bin/env python
#coding=utf8
#PYTHON_ARGCOMPLETE_OK 
from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
from functools import reduce
standard_library.install_aliases()
from builtins import input
from builtins import filter
from builtins import map
from builtins import str
from builtins import zip
from builtins import range
from builtins import *
from builtins import object
import argparse, argcomplete
import inspect 
import os, sys
import shutil 
import warnings 
import pprint
from tabulate import tabulate
import pickle as pickle
from collections import (OrderedDict, namedtuple)
from operator import itemgetter
import time
import scipy.integrate

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D

from scipy.signal import savgol_filter

import itertools 
from mpl_toolkits.mplot3d import Axes3D
try:
    import pandas as pd 
except:
    warnings.warn('pandas not imported')
import math 
from math import pi, cos, sin
import socket 
import platform 
import types

import importlib
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import sici
import unittest
from numpy import inf

from IPython.display import display

from merapy.utilities import dict_to_object , load , print_vars
from merapy.context_util import rpyc_load, rpyc_save, LOCAL_USERNAME, LOCAL_HOSTNAME 

#from mypy.brokest.task_center import submit_one 
#from mypy.brokest.brokest import queue, run_many, pack_runnable_msg, send_msg 


__all__ = ['MARKER_LIST', 'MARKER_CYCLE', 
    'ResultDB', 'ResultDB_idmrg', 'ResultDB_vmps', 'ResultDB_mera', 'ResultDB_tdvp', 
    'ResultDB_ed', 
    'BACKUP_STATE_DIR', 'RESULTDB_DIR', 'RESULTDB_ROOT', 'ResultDB_bethe_ansatz', 
    'display',  
    ]


HOME = os.path.expanduser('~')
HOME=HOME.replace('\\', '/') 
if 'cygwin' in HOME:  #properly deal with cygwin path
    assert 'cygwin64/home' in HOME
    HOME = HOME.replace('cygwin64/home', 'Users')

BACKUP_STATE_DIR = 'backup_tensor_dir'
RESULTDB_DIR = 'resultdb_dir'
if platform.system()=='Linux':
    #RESULTDB_ROOT = '/'.join([HOME, RESULTDB_DIR]) # attention,  I symbol linded RESULTDB_DIR from dropbox to home
    RESULTDB_ROOT = '/'.join([HOME,'dropbox', RESULTDB_DIR]) # attention,  I symbol linded RESULTDB_DIR from dropbox to home
else:
    RESULTDB_ROOT = '/'.join([HOME, 'Dropbox', RESULTDB_DIR], )
    #RESULTDB_ROOT = '/'.join([HOME, 'Onedrive', RESULTDB_DIR], )


def html_border(s, fontsize=20): 
    sss=  """
    <div id='page' style='width: 600px'>
      <h1 style='border:1.5px black solid; font-size:%(fontsize)dpx;'>
       %(s)s
      </h1>
    </div>
    """%vars() 
    return sss

def set_matplotlib_style():
    mpl.rcParams['figure.figsize'] = (4, 3)
    mpl.rcParams['axes.labelsize'] = 16
    #mpl.rcParams['legend.frameon'] = False  # whether or not to draw a frame around legend   
    v = mpl.__version__
    #if platform.system()=='Windows':
    if v<'2.0.0':
        raise  
        MATPLOTLIBRC = {
            u'grid.linestyle':'dotted', 
            #'legend.alpha': 0,    # this will make the box totally transanalysis_class 
            #'legend.edge_color': 'white',   # this will make the edges of the border white to match the background instead 
            'legend.frameon': 0,# whether or not to draw a frame around legend    
            'savefig.bbox': 'tight',  # default value is 'standard', it make figure craped when save 
            'figure.facecolor': 'white', 
            u'axes.titlesize':'x-large',
            }
    else:
        MATPLOTLIBRC = {
            'axes.grid':1, 
            u'legend.fontsize':'large', 
            u'grid.linestyle':'dotted', 
            #'Legend.alpha': 0,    # this will make the box totally transanalysis_class 
            u'legend.edgecolor': 'white',   # this will make the edges of the border white to match the background instead 
            u'legend.frameon': False,# whether or not to draw a frame around legend    
            u'savefig.bbox': 'tight',  # default value is 'standard', it make figure craped when save 
            u'figure.facecolor': 'white', 
            
            #u'axes.titlesize':'x-large',
            u'axes.titlesize':'20',
            
            }
        #MATPLOTLIBRC.pop('legend.alpha')
        #MATPLOTLIBRC.pop('legend.edege_color')
        #MATPLOTLIBRC.update({'legend.framealpha':0})
    mpl.rcParams.update(MATPLOTLIBRC)
    
    #print Line2D.marker
    marker_cycle = itertools.cycle(Line2D.markers)
    ## for reference 
    {
         None: 'nothing',
         0: 'tickleft',
         1: 'tickright',
         2: 'tickup',
         3: 'tickdown',
         4: 'caretleft',
         5: 'caretright',
         6: 'caretup',
         7: 'caretdown',
         '': 'nothing',
         ' ': 'nothing',
         '*': 'star',
         '+': 'plus',
         ',': 'pixel',
         '.': 'point',
         '1': 'tri_down',
         '2': 'tri_up',
         '3': 'tri_left',
         '4': 'tri_right',
         '8': 'octagon',
         '<': 'triangle_left',
         '>': 'triangle_right',
         'D': 'diamond',
         'H': 'hexagon2',
         'None': 'nothing',
         '^': 'triangle_up',
         '_': 'hline',
         'd': 'thin_diamond',
         'h': 'hexagon1',
         'o': 'circle',
         'p': 'pentagon',
         's': 'square',
         'v': 'triangle_down',
         'x': 'x',
         '|': 'vline'} 

    MARKER_LIST = [ 'o', 's', 'd',  'x', '*','+', 'D', '^',  '>', '.', '<', 'h']  #s=squre
    MARKER_CYCLE = itertools.cycle(MARKER_LIST)
    
    COLOR_LIST = [ 'black', 'red', 'green',  'yellow',  'olive',  'pink',
             'orange',  'brown',  'darkgrey',
            'silver',   'skyblue', 'gray',   'darkgreen',  'grey',  'violet',
            'blue', 'purple',  ]
    COLOR_CYCLE = itertools.cycle(COLOR_LIST)
    lines = ['', ' ', 'None', '--', '-.', '-', ':']
    LINE_CYCLE= itertools.cycle(lines)
    return MARKER_LIST, MARKER_CYCLE, MATPLOTLIBRC
    
MARKER_LIST, MARKER_CYCLE, MATPLOTLIBRC = set_matplotlib_style()


class AnalysisTools(object): 
    """
        an abstract base calss 
    """
    
    def filter_array(self, x, x_min=None, x_max=None, period=None, x1=None) :  
        """
            returns:
                arg
        """
        r = x1 
        if r is None: 
            if x_min is None: 
                #x_min = -np.inf  #this is not compatible with period below
                x_min = np.min(x)
            if x_max is None: 
                x_max = np.inf 
            temp = np.logical_and(x>=x_min, x<=x_max)
            #if period != 1:
            if period is not None: 
                #print_vars(vars(),  ['x-x_min', '(x-x_min)%period'])
                temp = np.logical_and(temp, (x-x_min)%period==0)
            arg = np.argwhere(temp).ravel()
        else:
            #arg = np.argwhere()
            arg = []
            a = 0
            for _r in r: 
                ii = 0
                for _x in x[a:]: 
                    if _r == _x:
                        arg.append(a + ii)
                        #print _r, a + ii   
                        a = ii 
                        break 
                    ii += 1  
        return arg 
    
    @staticmethod 
    def fit_curve(x, y, func, p0=None, method='lm', fit_range=None): 
        """
            just a wrapper 
            about curve_fit:
                Curve fitting is not always that straightforward.  
                p0 may be very useful in sometimes
            
        """
        if not isinstance(x, np.ndarray): 
            x = np.asarray(x)
        if not isinstance(y, np.ndarray): 
            y = np.asarray(y)
        #param, cov = curve_fit(func, x, y)    
        #res = func, param, cov 
        res= None
        try:
            param, cov = curve_fit(func, x, y, p0=p0, method=method)    
            res = func, param, cov 
        except RuntimeError as err:
            warnings.warn(str(err))
        except TypeError as err:
            warnings.warn(str(err))
        except Exception as err:
            warnings.warn(str(err))
            raise  
        return res
    
    def fit_line(self, line, func=None, data=None, x_min=None, x_max=None, 
            x_extra=None, plot_range=None,  period=None, x1=None, p0=None, method='lm'):
        if isinstance(func, str): 
            if func == 'EE_vs_L' : 
                func = lambda x, c, a: c/6.*np.log(x) + a
            elif func  == 'EE_vs_D': 
                #c=12*(1./k - 1)**(-2); k=1/(math.sqrt(12/c) + 1)
                #func = lambda x, c, a:  1/(math.sqrt(12/c) + 1)*np.log(x) + a 
                func = lambda x, k, a: k*np.log(x) + a 
            elif func == 'power' : 
                #func = lambda x, k, eta: k*x**(eta)
                func = lambda x, eta, k: k*x**(-eta)
            else: 
                raise 
        if data is not None:
            x, y = list(zip(*data))
        else:
            x, y = line.get_data()
        x = np.copy(x)
        y = np.copy(y)
        x = x.astype(float); y=y.astype(float)
        arg = self.filter_array(x, x_min, x_max, period, x1)
        x = x[arg]
        y = y[arg]
        if len(x) ==0 or len(y)==0:
            return None
        temp = self.__class__.fit_curve(x, y, func, p0=p0, method=method)
        if temp is None:
            return None
        else:
            func_, param, cov = temp #self.__class__.fit_curve(x, y, func)
        if x_extra is not None:  # append one point for expolatating 
            inc = True if np.all(np.diff(x)>0) else False
            if x_extra>max(x): 
                if inc:
                    x = np.append(x, x_extra)
                else:
                    x = np.insert(x, 0, x_extra)
            elif x_extra<min(x): 
                if inc:
                    x = np.insert(x, 0, x_extra)
                else:
                    x = np.append(x, x_extra)
            else: 
                raise ValueError((1./x, min(x), max(x)))
        if plot_range is not None:
            a, b = plot_range
            if a is not None and a < min(x):
                x = np.insert(x, 0, a)
            if b is not None and  b > max(x):
                x = np.append(x, b)
        
        y_fit = func_(x, *param)
        res = {'param': param, 'cov': cov, 'func': func, 'x': x, 'y': y_fit}
        return res 
    
    def fit_lines_many(self, ax, func=None, which_lines=None, plot_fit=True,  
            add_text=False, return_all_params=False, x_extra=None, plot_range=None,  rounding=4, fault_tol=True,  **kwargs): 
        """
            return:
                line_lable, param[0]
                
        """
        ll = list(ax.lines)
        which_lines = list(range(len(ll))) if which_lines is None else which_lines 
        #for l in ax.lines[:len(aa)]:
        temp = []
        fit_line_args= {a: kwargs.get(a) for a in 
                ['x_min', 'x_max', 'x_extra', 'period', 'x1'] }
        fit_line_args.update(x_extra=x_extra, plot_range=plot_range)
        for i, l in enumerate(ll): 
            if not i in which_lines: 
                continue 
            try: 
                a=eval(l.get_label())#[0]
            except: 
                a = None 
            if 0:
                try: 
                    res=self.fit_line(l, func, **fit_line_args) 
                    k = res['param'][0]
                except Exception as err: 
                    k = np.nan 
                    if fault_tol: 
                        warnings.warn(str(err))
                    else: 
                        raise
                    continue 
            else:
                res=self.fit_line(l, func, **fit_line_args) 
                
                if res is not None:
                    k = res['param'][0]
                else:
                    k = np.nan 
                    continue 
            if add_text: 
                if func == 'EE_vs_D' : 
                    k = 12*(1./k - 1)**(-2)  #central charge 
                #label = ' '.join([l.get_label(), '%s'%str(k)[:4]]) 
                label = ' '.join([l.get_label(), '%1.2f'%k])
                l.set_label(label)
            #temp.append((a, round(k, 2)))
            if return_all_params: 
                temp.append((a, res['param']))
            else: 
                temp.append((a, round(k, rounding)))
          
            if plot_fit: 
                args= kwargs.copy()
                color = kwargs['color'] if 'color' in kwargs else l.get_color()
                marker = kwargs.get('marker', None)
                args.update(color=color, marker=marker)
                #print_vars(vars(),  ['args'])
                #_ = self._plot.im_func(None, res['x'], res['y'], 
                #        ax=ax, color=l.get_color(), label='', marker=None)        
                _ = self._plot.__func__(None, res['x'], res['y'], 
                        ax=ax, **args)
                #if x_extra is not None:
                #    tic= ax.get_xticks().tolist()
                #    print_vars(vars(),  ['tic'])
                #    pass
                
        return temp 
    
    def fig_layout(self, ncol=1, nrow=1, size=(5, 4), dim=2): 
        if size is None and ncol == 1 and nrow == 1: 
            size = (6, 5)
        #size= size if size is not None else (3.5, 2.5)
        size= size if size is not None else (5, 4)
        size = ncol*size[0], nrow*size[1]
        fig = plt.figure(figsize=size)
        iii=0
        for i in range(nrow): 
            for j in range(ncol): 
                iii+=1
                fig.add_subplot(nrow, ncol, iii, projection=None if dim==2 else '3d')
        fig.tight_layout()
        axes = fig.axes 
        if len(axes)==1: 
            axes= axes[0]
        return fig, axes 

    def add_ax(self, fig=None, size=None, direct='h'): 
        if fig is None: 
            fig, ax=self.fig_layout(size=size)
            return ax 
        
        n = len(fig.axes)
        #the following excludes bars of imshows 
        #temp = [a for a in fig.axes if a.get_aspect()=='auto' or a.get_aspect()<5]  
        temp = []
        bars= []
        for i, a in enumerate(fig.axes):
            if a.get_aspect()=='auto' or a.get_aspect()<5: 
                temp.append((i, a))
            else: 
                bars.append((i, a))
        n1 = len(temp)
        
        if direct == 'h' : 
            if 1: 
                fig.set_figwidth(fig.get_figwidth()*(n + 1.)/n)
                for i, ax in enumerate(fig.axes): 
                    ax.change_geometry(1, n+1, i+1)
                ax = fig.add_subplot(1, n + 1, n + 1)
            if 0: 
                n = n1
                #fig.set_figwidth(fig.get_figwidth()*(n + 1.)/n)
                print_vars(vars(),  ['fig.get_figwidth()'])
                fw = fig.get_figwidth()*(n1 + 1.)/n1
                print_vars(vars(),  ['fw'])
                fig.set_figwidth(fw)
                #for i, ax in enumerate(fig.axes): 
                for i, ax in temp: 
                    ax.change_geometry(1, n+1, i+1)
                print_vars(vars(),  ['fig.get_figwidth()'])
                ax = fig.add_subplot(1, n + 1, n1 + 1)
                print_vars(vars(),  ['fig.get_figwidth()'])
            if 0: 
                #fig.set_figwidth(fig.get_figwidth()*(n + 1.)/n)
                fw = fig.get_figwidth()*(n1 + 1.)/n1
                fig.set_figwidth(fw)
                temp = [t[1] for t in temp]
                bars = [t[1] for t in bars]
                for i, ax in enumerate(temp): 
                    ax.change_geometry(1, n+1, i+1)
                ax = fig.add_subplot(1, n+1, n1 + 1)
                for i, ax in enumerate(bars): 
                    ax.change_geometry(1, n+1, n1 + 1 + i)
                
                print_vars(vars(),  ['fig.get_figwidth()'])
                
        elif direct == 'v' : 
            fig.set_figheight(fig.get_figheight()*(n + 1.)/n)
            for i, ax in enumerate(fig.axes): 
                ax.change_geometry(n+1, 1, i+1)
            ax = fig.add_subplot(n+1, 1, n + 1)
        
        return ax 
        
    def plot_alpha_2d(self, aa=None, xparam=None, yparam=None, rotate=False, **kwargs):
        aa = aa if aa is not None  else self.alpha_list 
        n = len(self.param_list)
        if len(aa)>0: 
            xparam = xparam if xparam is not None else self.param_list[0]
            yparam = yparam if yparam is not None else self.param_list[1]
            xparam_id = self.param_list.index(xparam)
            yparam_id = self.param_list.index(yparam)
            x = [a[xparam_id] for a in aa]
            y = [a[yparam_id] for a in aa]
            #aa = [a[: n] for a in aa]
            #x, y = zip(*aa)
        else: 
            x, y = np.nan, np.nan 
        #xlabel = self.param_list[0]
        #ylabel = self.param_list[1]
        xlabel = xparam
        ylabel = yparam
        #kwargs= kwargs.copy()
        if rotate: 
            temp = x
            x = y; y = temp 
            temp = xlabel
            xlabel = ylabel; ylabel = temp 
        #kwargs.update(lw=0, xlabel=xlabel, ylabel=ylabel, marker=kwargs_orig.get('marker'), 'x')
        args=dict(lw=0, xlabel=xlabel, ylabel=ylabel, marker='x')
        args.update(kwargs)
        self._plot(x, y, **args)

    def plot_derivative(self, line_list, **kwargs):
        pass
        ll = line_list
        for l in ll:
            x, y = l.get_data()
            y_diff = np.diff(y)
            x_diff = np.diff(x)
            y = y_diff/x_diff 
            label = l.get_label()
            self._plot(x[:-1], y, label=label,  **kwargs)
        
    plot_diff = plot_derivative #def plot_diff

    def find_lines_extreme(self, ax, which='max', add_text=True, 
            find_range=None, 
            zoom_scale=None, font_dict=None, rounding=3): 
        temp=[]        
        fd = {'color': 'r', 'size': 14} 
        if font_dict is not None: 
            fd.update(font_dict)
        for l in ax.lines:
            x, y= l.get_data()
            if find_range is not None: 
                x_ = np.logical_and(x>=find_range[0], x<=find_range[1])
                arg = np.argwhere(x_).ravel()
                x = x[arg]
                y = y[arg]
            if np.all(np.isnan(y)):
                warnings.warn("omit empty line %s in find_extreme "%(l.get_label()))
                continue 
            if which == 'max' : 
                #min_ind =y.argmax()
                min_ind =np.nanargmax(y)
            else: 
                #min_ind =y.argmin()
                min_ind =np.nanargmin(y)
            min_location = x[min_ind]
            min_val = y[min_ind]
            if add_text: 
                ax.text(min_location, min_val, round(min_location, rounding), fontdict=fd,  )
            temp.append((l.get_label(), round(min_location, rounding), min_val))
        if temp: 
            labels, loc, val= list(zip(*temp))
            labels= list(map(str, labels))  # by default, labels are unicode, u'...', I remove u signiture
        else: 
            labels, loc, val = None, None, None 
            
        if zoom_scale and len(ax.lines)==1: 
            #temp=temp[1][0]
            temp = loc[0]
            #delta=0.5
            delta = zoom_scale
            x1=temp-delta; x2=temp+delta
            ax.set_xlim(x1, x2)
            l=ax.lines[0]
            x,y=l.get_data()
            ind=np.bitwise_and(x>x1, x<x2)
            y=y[ind]
            y1=np.nanmin(y); y2=np.nanmax(y)
            #print y1, y2
            ax.set_ylim(y1, y2)            
        return labels, loc, val 
    
    def find_points_in_line_bac(self, line, ymin, ymax):
        x, y = line.get_data()
        arg = np.bitwise_and(y>=ymin, y<=ymax)
        return x[arg], y[arg]
    
    def find_points_in_line(self, line, range=None, close_to=None):
        x, y = line.get_data()
        if range is not None: 
            ymin, ymax = range
            arg = np.bitwise_and(y>=ymin, y<=ymax)
            return x[arg], y[arg]
        if close_to is not None:
            diff = np.abs(y - close_to)
            arg = np.nanargmin(diff)  #ignore nan 
            return x[arg], y[arg]
       

    def plot_sign(self, ax, line_id=0, center=0.0, magnitude=None, **kwargs):
        l = ax.lines[line_id]
        x, y = l.get_data()
        magnitude = magnitude if magnitude is not None else np.max(np.abs(y)) *0.8 
        yy = np.sign(y)*magnitude  + center 
        self._plot.__func__(None, x, yy, ax=ax, **kwargs)
    
    def compare_two_lines(self, line1, line2, plot_sign=False, **kwargs): 
        """
           plot y1-y2 
        """
        x1, y1 = line1.get_data()
        x2, y2 = line2.get_data()
        if 1:   #here asumes no duplicate elements in x1 or x2  
            #ind12 = np.intersect1d(x1, x2, )
            ind12 = np.in1d(x1.round(14), x2.round(14))   #arg of ele of x1 that is in x2
            ind21 = np.in1d(x2.round(14), x1.round(14))
            x_common = x1[ind12]
            y12 = y1[ind12]
            y21 = y2[ind21]
            y_diff = y12 -y21
        #axes= line1.get_axes()   
        axes= line1.axes
        if 'ylabel' not in kwargs: 
            ylabel = axes.get_ylabel()
            ylabel = '_'.join([ylabel, 'diff'])
            kwargs['ylabel'] = ylabel
        if 'xlabel' not in kwargs: 
            kwargs['xlabel'] = axes.get_xlabel()
        #self._plot.im_func(None, x1, y1-y2, **kwargs)
        if 'label' not in kwargs: 
            label = '-'.join([line1.get_label(), line2.get_label()])
            kwargs['label'] = label 
        if plot_sign: 
            kwargs.update(yfunc=np.abs, yscale='log')
        
        kwargs['yfunc'] = kwargs.get('yfunc', np.abs)
        kwargs['yscale'] = kwargs.get('yscale', 'log')
        
        fig=self._plot.__func__(None, x_common, y_diff, **kwargs)
        if plot_sign: 
            self.plot_sign(fig.axes[-1])
    
    
    def remark(self, s, style=None, color=None,  border=0): 
        import IPython.display as display
        if color is not None: 
            s= '<font color=%s> %s </font>'%(color, s)
        if border: 
            s= html_border(s)
        return display.HTML(s)
    
    @staticmethod 
    def ax_fewer_legend(ax, line_id_list):
        n = len(ax.lines)
        line_id_list = [l for l in line_id_list if l <  n]
        ll=[ax.lines[i] for i in line_id_list]
        temp=[(l, l.get_label()) for l in ll if '_line' not in l.get_label()]
        x, y= list(zip(*temp))
        ax.legend(x, y, loc=0)
    
    @staticmethod 
    def ax_hide_lines(ax, line_id_list, reverse=False):
        n = len(ax.lines)
        line_id_list = [i for i in line_id_list if i<n]
        if reverse: 
            ll = list(range(len(ax.lines)))
            temp = list(line_id_list)
            line_id_list = [i for i in ll if i not in temp]
        lines  = [ax.lines[i] for i in line_id_list]
        for l in lines:
            l.set_marker(None)
            l.set_alpha(0)
            l.set_label('')
        ax.legend(loc=0)        

    def search_str(self, input_str, str_list=None, 
            only_first=False, only_one=False, assert_found=False): 
        """
        """
        ss = input_str.split(' ')
        ss_len = len(ss)
        found = False
        res= []
        for fn in str_list: 
            match_count = 0
            for n in ss: 
                if n not in fn: 
                    break 
                else: 
                    match_count += 1
            if match_count == ss_len: 
                res.append(fn)
                found = True 
                if only_first: 
                    break 
        if only_one: 
            if len(res)>0: 
                res= [r for r in res if r == input_str]
        if assert_found: 
            assert len(res)>0, '"{}" matches no field'.format(input_str)
        return res 

    def set_xtics_for_inverse_N(self, ax):
        tic= ax.get_xticks().tolist()
        if 1:
            tic0 = tic[:1]  # it may be 0, 1/0 raises 
            tic = tic[1:]
            tic = ['%s\n%d'%(i, 1./(i)) for i in tic]
            _ = ax.set_xticklabels(tic0 + tic)    
        else:  # not good 
            if tic[0] == 0.0:
                tic[0] = 1e-10
            tic = ['%s\n%d'%(i, 1./(i)) for i in tic]
            _ = ax.set_xticklabels(tic)    

    def get_contour_line_data(self, cb, i):
        res= []
        assert i < len(cb.collections)
        #for i in range(n):   
        p = cb.collections[i].get_paths()[0]
        v = p.vertices
        #print v.shape
        x = v[:,0]
        y = v[:,1]
        
        return x, y

    def smooth(self, data, window_size, order):
        res = savgol_filter(data, window_size, order)
        return res 

class AnalyticFormular(object):
    """
       some formular for comparing with numerics or extracting values
    """
    @staticmethod
    def corr_luttinger( which='zz', nu=0.5):
        """  
            correlation for XXZ (or t-V) model
            ref. 
                Giamarchi's book p.168
                Haldane 1981 effect:  eq.7-9
        """
        cos = np.cos   # use math.cos may raise due to data type
        _nu={0.5:0.5, 0.33:1./3, 0.25:0.25}[nu]
        m = 0.5*(2*_nu - 1.0)   #magnetization
        # kf/pi = nu
        
        if which == 'zz':
            if nu==0.5:
                def func(x,  K, C1, C2): #eq.6.38 first line
                    res = 0
                    #res += C1/x**2   # C1 can be specified as below
                    res += 4*(-1.)*K/(2*pi**2)*(1./x**2)  
                    res += C2*(-1.)**x * (1/x)**(2*K)   #(-1.)**x i.e. cos(2kf)
                    #res += C2*(-1.)**x * (1./x)**(2*K)* np.log(x)**0.5   #with log correction
                    return res
            else:
                def func(x, K, C1, C2):
                    res=0
                    res += 4*m**2  # note diff with Giamarchi's book which uses m**2, because Sz takes value  +1/-1 here
                    #res += C1/x**2  
                    res += 4*(-1.)*K/(2*pi**2)*(1./x**2)  
                    #res += C2*cos(pi*(1+2*m)*x)*(1./x)**(2*K)  #since one just has (1+2*m)=2*nu, use the following
                    res += C2*cos(2*pi*_nu*x)*(1./x)**(2*K)
                    return res
        
        elif which =='pm':
            #if nu==0.5:
            #    def func(x,  K, C1, C2):   #eq.6.38 sec line  
            #        """  
            #            this has been verified thoroughly
            #        """
            #        res = 0
            #        res += C1*(1/x)**(2*K+1./(2*K)) 
            #        res += C2*(-1.)**x * (1./x)**(1/(2*K))
            #        return res
            #else:
            def func(x,  K, C1, C2 ):   
                """
                    
                """
                res = 0
                res += C1*cos(2*pi*m*x)*(1./x)**(2*K+1./(2*K))
                res += C2*(-1.)**x * (1./x)**(1./(2*K))  
                return res
                

        return func
    
    @staticmethod
    def exact_corr_cdag_c(r):   #exact_corr_xx_model
        """
            the one particle reduced tensity matrix
            for spinless free fermion ( a.k.a xx model)
            <c^+_{i} c_j> which can be mapped into spin operator sp sm with JS
            string in between 
            
            ref. 
                Latorre 2004 eq. C8 C9 C10
                Evenbly 2011 eq.46
                
        """
       
        if r != 0: 
            res=-sin(pi*r/2.)/(pi*r)
        else:
            #issue: I dont know here should be 0.5 or -0.5
            #res= -0.5
            res= 0.5
        return res
    
    @staticmethod
    def luttinger_param_theory(Jzz):
        """
            theoritical value of K for XXZ model
            Jzz  =  Delta
        """
        #luttinger_param_theory=lambda x: np.pi/(2* np.arccos(-x))
        
        res = np.pi/(2* np.arccos(-Jzz))
        return res

    def kac_rescale(N, alpha):
        """
            Kac rescaling for LR models with power law interaction 
        
        """
        if 0:
            temp = [abs(i-j) for i in range(N)  for j in range(N) if i<j]
            temp = np.asarray(temp)
            res = np.sum(temp**(-alpha)) / (N-1)
        else:
            rr = np.arange(1, N, 1)
            temp = rr**(-alpha)
            res= np.sum(temp)
            
        return res 

class ResultDB(OrderedDict, AnalyticFormular,  AnalysisTools): 
    """
        a DB to store measurement results 
        
        the structure of the db,  that gives the instruction of the key chains to fetch recs
        version log 
            1.0
                mera STRUCT = {
                            'entanglement': [shape, iter, layer, num_sites], 
                            }
                MPS STRUCT = {
                            'entanglement': ['shape', 'iter', num_sites], 
                            }
            1.2:  'iter' is removed 
                MPS STRUCT = {
                            'entanglement': ['shape', 'num_sites], 
                            }
            
        design: 
            fields are of two kinds,  
            first,  like 'version', 'status' is not dependent on sh 
            second, physical quantities, like 'energy', 'correlation'
        
    """
    DBNAME = 'RESULT.pickle.db'   #fix the name of db
    ALGORITHM = 'all'
    AUTO_UPDATE = True 
    VERSION = 1.2 
    ALL_FIELD_NAMES = ['energy', 'magnetization', 'concurrence', 'concurrence_2']
    def __init__(self, parpath, dbname=None, version=None, algorithm=None, model_param=None,  
            use_local_storage=False, create_empty_db=False,  upgrade=0, 
            analysis_class=None):
        """
            params: 
                algorithm: algorithm is necessary, or else it is difficult to determine which alg
                param: model_param 
        """
        self.use_local_storage = use_local_storage
        self.parpath = parpath
        self.analysis_class = analysis_class
        dbname = dbname if dbname is not None else self.DBNAME
        self.dbname = dbname
        self.model_param = model_param
        self.param = model_param # a shorter name 
        self.path = '/'.join([parpath , dbname])
        
        #self.inited = False  
       
        if 0: 
            if not os.path.exists(parpath):  # these lines are useful, when merge remote db
                self.create_parpath()
       
        try: 
            db = rpyc_load(self.path, use_local_storage=use_local_storage)
        except IOError as err: 
            db = {}
            if create_empty_db: 
                ResultDB.create_empty_db(self.path, use_local_storage) 
        except Exception as err: 
            raise 
        OrderedDict.__init__(self, db)   
        
        if algorithm is not None : 
            self['algorithm'] = algorithm
        
        #if 0:
        #    self.version = version
        #    self.upgrade = upgrade
        #   
        #    version = self.version if self.version is not None else self.get('version')
        #    self.version = version
        #    self['version'] = version
         
        if 'version' not in self: 
            self['version'] = self.VERSION 
        if 'status' not in self:   # means is the data correct。these make init slow, better make it lazy eval 
            self['status'] = 'good'   # default is good 

        #if self['version'] >= 1.0:  
        #    self.fetch_easy = self.fetch_easy_new
        if self['version'] < 1.0:  
            self.fetch_easy = self.fetch_easy_old
        
        self._last_modify_time = time.time()   # used for auto refresh db 

        if 0: 
            if upgrade:   #this will be deprecated 
                if self.get('version') <= 1.0 and dbname != self.dbname + '.new': 
                    self.upgrade_db()
                    
        elif self.AUTO_UPDATE:   
            if self['version']<self.VERSION: 
                self.update_db_structure()
    
    
    #@staticmethod  del this line
    def algorithm_name_to_rdb(alg_name):
        mapping = {'mera':ResultDB_mera, 
                'mps':ResultDB_vmps, 
                'vmps':ResultDB_vmps, 
                'idmrg':ResultDB_idmrg, 
                'tdvp':ResultDB_tdvp, 
                }
        return mapping[alg_name]

    @property
    def state_parpath(self): 
        dir = self.parpath.replace(RESULTDB_DIR, BACKUP_STATE_DIR)
        dir = dir.replace('dropbox', '')
        return dir 
    
    @property
    def dir_name(self):
        res = os.path.basename(self.parpath)
        return res 
    
    parpath_name = dir_name  # def parpath_name 
    
    def parpath_map(self, dir):
        if 'C:' in dir:
            xxx = [ ('C:', ''), 
                    ('Users', 'home'), 
                    ('zhihua', LOCAL_USERNAME),  
                    ('Dropbox', 'dropbox')]
            for a, b in xxx:
                dir=dir.replace(a, b)
            dir = os.path.normpath(dir)
            dir = dir.replace('\\', '/')
        return dir 
    
    @property
    def parpath_center(self):
        dir = self.parpath_map(self.parpath)
        return dir

    @property
    def last_modify_time(self):
        try:
            return os.path.getmtime(self.path)
        except OSError:
            return 0
        except Exception:
            raise  

    def __repr__(self): 
        return list(self.keys())
        temp = []
        temp.append(self.parpath)
        temp.append(str(list(self.keys())))
        res= '\n'.join(temp)
        return res 
    
    def __str__(self): 
        temp = []
        temp.append(self.parpath)
        temp.append(str(list(self.keys())))
        res= '\n'.join(temp)
        return res 
    
    def __reduce__(self):
        """
            without this, it will raise if picking or unpickle resultdb instantce 
            see also 
                https://www.google.com/webhp?sourceid=chrome-instant&ion=1&espv=2&ie=UTF-8#newwindow=1&safe=off&q=python+subclass+ordereddict+pickle+   
                http://stackoverflow.com/questions/22690729/pickling-an-ordereddict-derived-object
                http://stackoverflow.com/questions/19855156/whats-the-exact-usage-of-reduce-in-pickler          
                
            Attention: 我写这个的目的仅仅是为了在远程计算时，可以把db作为args之一序列化传过去。
            事实上，在pickle时，db中的多数内容都没有保存，只是保存了如何组装重建db, 这正好！！
            如果想保存db内容到pickle，需要修改本函数，参考上面的链接
                    
        """
        return (self.__class__, (self.parpath, ))
    
        
    if 0: 
        def __getitem__(self, k): 
            """
                
            """
            print('iiiii')
            if not self.inited: 
                self.init_db()
                print('iiiii')
            else: 
                #return self[k]
                return OrderedDict.__getitem__(self, k)

        def __getitem__(self, k): 
            return self._db[k]
        
        def __setitem__(self, k, v): 
            self._db[k] = v
            pass
    
    def load_S(self, sh, path=None, to_obj=False, use_local_storage=False):
        """
            todo: load 时，应该支持 update state 
        """
        if sh[0] == 'max' :
            N = sh[0]
            temp = self.get_shape_list( only_return_max=1, from_energy_rec=0)
            temp = [t[0] for t in temp]
            if not temp:
                N = -1
            else:
                N = max(temp)
            sh = (N, sh[1])
        
        if sh[1] ==  'max':
            if 0: # this is not precise 
                sh = sh[0], self['dim_max'].get(sh[0], 0)
            else: 
                N = sh[0]
                temp = self.get_shape_list(N=N, only_return_max=1, from_energy_rec=0)
                if not temp: # sh  == [] :
                    sh = (N, -1)
                else:
                    sh = temp[0]
        
        path = self.shape_to_backup_path(sh)
        try:
            res= rpyc_load(path, use_local_storage=use_local_storage)
        except IOError as err: 
            warnings.warn(str(err))
            res = None    
        if res is not None and to_obj: 
            res= dict_to_object(res)
        return res
    
    def reload_class(self): 
        module_name =self.__class__.__module__
        #print module_name
        module=importlib.import_module(module_name)
        reload(module)
        name = self.__class__.__name__
        self.__class__=getattr(module, name)
        print('class %s is reloaded'%(name, ))

    def create_parpath(self): 
        cmd = 'mkdir %s'%self.parpath
        os.system(cmd)
        msg = 'parpath not exist, so %s'%cmd

    def has_shape(self, sh, system_class=None, from_energy_rec=False): 
        if from_energy_rec: 
            return self.has_key_list(['energy', sh])
            
        if system_class is None:
            system_class = self.__class__ 
        fn = self.shape_to_backup_fn(sh)
        path = '/'.join([self.state_parpath, fn])
        return os.path.exists(path)

    def has_key_list(self, key_list, info=1): 
        temp = self
        for k in key_list: 
            #if not hasattr(temp, 'has_key'):   # val of rec
            #    if info>0: 
            #        print 'not has this key: %s'%(k, )
            #    return False
            #if not hasattr(temp, 'has_key') or not temp.has_key(k): 
            if k not in temp: 
                #if info>0: 
                #    print 'not has this key: %s'%(k, )
                return False
            temp = temp[k]
        
        return True
    
    def fetch_from_key_list(self, key_list, default=None, rec=None, fault_tol=True):
        if rec is None:
            rec = self
        try: 
            for key in key_list: 
                rec = rec[key]
        except KeyError as err: 
            if fault_tol: 
                warnings.warn(str(err))
                #print 'keys are', rec.keys()
                return default
            else: 
                raise
        #finally:
        #    return default
        return rec
    
    @staticmethod 
    def make_nested_dict(key_list, old_dic=None,  val=None, info=0): 
        if old_dic is None: 
            res = OrderedDict()
        else: 
            res= old_dic 
        temp = res 
        for k in key_list[: -1]: 
            if k not in temp: 
                temp[k] = OrderedDict()
            #temp[k] = OrderedDict()
            temp = temp[k]
        
        last_key = key_list[-1]
        if val is None: 
            temp[last_key] = OrderedDict()
        else: 
            temp[last_key] = val
        return res 
    
    make_key_chain  =  make_nested_dict   #def make_key_chain
    
    def add_key_list(self, key_list, val=None, verbose=0): 
        if 1: 
            temp = self
            for k in key_list[: -1]: 
                if k not in temp: 
                    temp[k] = OrderedDict()
                temp = temp[k]
            last_key = key_list[-1]
            if val is None: 
                temp[last_key] = OrderedDict()
            else: 
                temp[last_key] = val
            #print self.has_key_list(key_list)
            #exit()
            if verbose: 
                print('add_key_list %s'%(key_list, ))
        else:  # replace above with make_nested_dict in future 
            temp=self.__class__.make_nested_dict(key_list[1: ], val=val)
            
            self[key_list[0]] = temp 
    
    add_key_chain = add_key_list 
    insert_key_list = add_key_list   #def insert_key_list
    
    @classmethod
    def create_empty_db(cls, path, use_local_storage=False): 
        if 0: 
            out = open(path,  'wb')  #raise #useless , just use 'ab+' mode
            db = OrderedDict()
            pickle.dump(db, out)
            out.close()
        db = OrderedDict()
        rpyc_save(path, db, use_local_storage=use_local_storage)
        print('ResultDB %s has been created'%path)
    
    def get_fn_list(self, dir=None):
        dir = dir if dir is not None else self.state_parpath
        #issue: this could be slow if do this each time , find better way later 
        #dir = dir.replace('resultdb_dir', 'backup_tensor_dir')
        if not os.path.exists(dir):
            return []
        nlist = os.listdir(dir) 
        pickle_files = [i for i in nlist if "pickle" == i.split('.')[-1] and '=' in i]
        #if not self.get('algorithm')=='mps': 
        #    pickle_files= [i for i in pickle_files if 'trans' not in i]
        
        return pickle_files
    
    def get_time_serials(self, sh, attr_name=None, info=0, force=0): 
        #if rec is not None: 
        if 'time_serials' in self: 
            ts = self['time_serials'].get(sh)
        else: 
            self['time_serials'] = {}
            ts= None 
        if ts is None or force : 
            s= self.load_S(sh)
            if s is None:
                return None 
            ts= s['time_serials']
            if 0:
                ts = pd.DataFrame(ts).T
                self['time_serials'][sh] = ts 
                self.commit()
            else:
                self['time_serials'][sh] = ts 
        res = ts 
        if attr_name is not None: 
            res= [ts[i][attr_name] for i in ts]
            #res = ts[attr_name]
        return res 
    
    def parse_fn(self, fn):
        ND=fn.split('.')[0]
        N, D = ND.split('-')
        N = N.split('=')[1]; N=eval(N)
        D = D.split('=')[1]; D=eval(D)
        return N, D 

    def get_shape_list(self, N=None, D=None, sh_max=None, sh_min=None, 
            from_energy_rec=True, field='energy',  only_return_max=False, only_return_N=False, only_return_D=False, ): 
        if socket.gethostname()==LOCAL_HOSTNAME and not from_energy_rec: 
            fn = self.get_fn_list()
            sh = [self.parse_fn(f) for f in fn]
        else:
            #warnings.warn('get_shape_list using energy rec')
            sh = list(self.get(field, {}).keys())
       
        #if N is not None:
        #    sh_min = (N, 0)
        #    sh_max = (N, np.inf)
        #if D is not None: 
        #    sh_min = (0, D)
        #    sh_max = (np.inf, D)
        if N is not None:
            sh=[x for x in sh if x[0]==N]
        if D is not None:
            sh=[x for x in sh if x[1]==D]
        
        if sh_max is not None : 
            sh=[x for x in sh if x[0]<=sh_max[0] and x[1] <= sh_max[1]]
        if sh_min is not None : 
            sh=[x for x in sh if x[0]>= sh_min[0] and x[1]>= sh_min[1]]
        
        if only_return_max:
            NN = [i[0] for i in sh]
            NN = list(set(NN))
            sh = [(N, max([x for x in sh if x[0]==N])[1]) for N in NN]
            #sh = [(N, self.get_dim_max_for_N(N)) for N in NN]
        if only_return_N:
            sh = [i[0] for i in sh]
            sh = list(set(sh))
        if only_return_D:
            sh = [i[1] for i in sh]
            sh = list(set(sh))
        sh.sort()
        return sh
    
    get_mera_shape_list = get_shape_list 
    
    def get_dim_max_for_N(self, N, return_N=False, update_db=False, 
            from_energy_rec=True, 
            force=False, info=0): 
        """
            params:
                N: can be integer N or tuple (N, D)
        
        """
        
        if isinstance(N, tuple):  
            N = N[0] 
            return_N = True
        if N == 'max' :
            temp=self.get_shape_list(from_energy_rec=from_energy_rec, only_return_N=1)        
            N = max(temp) if temp else -1
            
        try: 
            if force: 
                raise KeyError 
            if not return_N:
                return self['dim_max'][N]
            else:
                return N, self['dim_max'][N]
        except KeyError: 
            temp=self.get_shape_list(N=N, from_energy_rec=from_energy_rec)
            if len(temp)>0:  
                D = max(temp)[1] 
            else: 
                D = 0 
            if update_db: 
                if 'dim_max' not in self: 
                    self['dim_max'] = {}
                self['dim_max'][N] = D
                self.commit(info=info)
            if not return_N:
                return D 
            else:
                return N, D 
    
    get_dim_max = get_dim_max_for_N  #def get_dim_max a shorter name 
    
    def get_energy_fluc_count(self, sh, bins=None): 
        if bins is None: 
            bins = [-1, 0,  1e-8, 1e-7, 1e-6, 1e-5, 1] 
     
        diff=self.get_energy( sh=sh, which='diff', iter_max=None)
        count, bin = np.histogram(diff, bins=bins)
        return count, bins
    
    
    def commit(self, use_local_storage=None, info=0):
        #issue: should change to 'ab+?'
        if 0: 
            out = open(self.path, 'wb')
            #out = open(self.path, 'wb')
            #pickle.dump(self, out)  #this is rather critical, or else got error
            pickle.dump(OrderedDict(self), out)
            out.close()
        
        use_local_storage = use_local_storage if use_local_storage is not None else self.use_local_storage
        path = self.path
        if use_local_storage:
            path = self.parpath_map(path)
        rpyc_save(path, OrderedDict(self), use_local_storage=use_local_storage)
        if info>0: 
            temp = str(path)
            print('\tmodification in ...%s has been commited'%(temp, ))
    
    def _fetch(self, dim, num_of_layer=None): 
        keys= iter(self.keys())
        #keys= filter(lambda x: x[0]==dim and x[1]==num_of_layer, keys)
        keys= [x for x in keys if x[0]==dim]
        if num_of_layer is not None: 
            keys= [x for x in keys if x[1]==num_of_layer]
        if len(keys)<1: 
            return None
        else: 
            #print 'kkk', keys
            #print self.keys()
            try: 
                res = max(keys, key=itemgetter(2))  #records with max S.iter
            except: 
                print('parpath %s'%(self.parpath, ))
                print('keys are, %s'%keys)
                raise
            return res
    
    def fetch(self, field_name, dim, num_of_layer=4, step=None, info=1): 
        key = self._fetch(dim, num_of_layer)
        if key is None: 
            if info>0: 
                msg = 'key (%s) not found in ResultDB  at %s\n keys are %s'%(
                    (dim, num_of_layer), self.parpath,  list(self.keys()))
                print(msg)
            return None
        else: 
            if field_name in self[key]: 
                return  self[key][field_name]
            else: 
                if info>0: 
                    msg = "rec not found for field_name '%s' at %s"%(field_name, key)
                    print(msg)
                return None

    def fetch_easy_old(self, field_name, mera_shape,  sub_key_list=None, mera_shape_optional=None, 
             auto_meassure=False, default=None,  info=0, fault_tolerant=1): 
        """
            interface to make fetch easier
        """
        
        dim, layer = mera_shape[:2]
        if 1: 
            mapper = {
                    'EE':  ['central_charge', 'EE'  ], 
                    'EE1': ['central_charge', 'EE', 1], 
                    'EE2': ['central_charge', 'EE', 2], 
                    }
            if field_name in mapper: 
                temp = mapper[field_name]
                field = temp[0]
                key_list = temp[1: ]
            else: 
                field = field_name
                key_list = None
            if sub_key_list is not None: 
                if key_list is None: 
                    key_list = sub_key_list
                else: 
                    key_list += sub_key_list
        
        
        rec = self.fetch(field, dim, num_of_layer=layer, info=info)
        if 0: 
            if key_list is not None and rec is not None: 
                try: 
                    for key in key_list: 
                        rec = rec[key]
                except: 
                    if fault_tolerant: 
                        print('error in fetch_easy_old')
                        print('keys are', list(rec.keys()))
                    else: 
                        raise
        else:
            rec=self.fetch_from_key_list(key_list, rec=rec, 
                    default=default, 
                    fault_tol=fault_tolerant)
        return rec

    def fetch_easy_v1p0(self, field_name, mera_shape,  sub_key_list=None, mera_shape_optional=None, 
             auto_meassure=False, info=0, fault_tolerant=1): 
        """
            fetch_easy for version 1.0
            for backward compatibility. remove this meth after some time
        """
        
        try: 
            if mera_shape[1] == 'max' : 
                mera_shape = mera_shape[0], self['dim_max'].get(mera_shape[0], 0)
               
            rec = self[field_name][mera_shape]
            #key = next(reversed(rec))  #last iter step
            #rec = rec[key]
        except KeyError as err: 
            if fault_tolerant: 
                if info>0: 
                    print('key not found:', err)
                return None
            else:
                raise
        except: 
            raise
        
        if 1: 
            mapper = {
                    'EE':  ['central_charge', 'EE'  ], 
                    'EE1': ['central_charge', 'EE', 1], 
                    'EE2': ['central_charge', 'EE', 2], 
                    }
            if field_name in mapper: 
                temp = mapper[field_name]
                field = temp[0]
                key_list = temp[1: ]
            else: 
                field = field_name
                key_list = None
            if sub_key_list is not None: 
                if key_list is None: 
                    key_list = sub_key_list
                else: 
                    key_list += sub_key_list
        
        if key_list is not None and rec is not None: 
            try: 
                for key in key_list: 
                    rec = rec[key]
            except: 
                if fault_tolerant: 
                    if info>0: 
                        print('error in fetch_easy')
                        print('existent keys are %s; required key_list is %s'%(list(rec.keys()), key_list))
                    return None 
                else: 
                    raise
        return rec

    def fetch_easy(self, field_name, mera_shape,  sub_key_list=None, mera_shape_optional=None, 
             auto_meassure=False, default=None,  info=0, fault_tolerant=1): 
        """
            fetch_easy for version >= 1.2
        """
        try: 
            if mera_shape[1] == 'max' : 
                mera_shape = mera_shape[0], self['dim_max'].get(mera_shape[0], 0)
            elif mera_shape[1] == 'field_max':  #the max dim for each field may differ, then may use this
                N = mera_shape[0]
                ss = list(self.get(field_name, {}).keys())
                dd = [s[1] for s in ss if s[0]==N ]
                D = max(dd) if dd else 0     
                mera_shape = N, D
            rec = self[field_name][mera_shape]
        except KeyError as err: 
            if fault_tolerant: 
                if info>0: 
                    print('key not found:', err)
                return default
            else:
                raise err
        except: 
            raise
        
        if sub_key_list is not None and rec is not None: 
            rec=self.fetch_from_key_list(sub_key_list, rec=rec, 
                    default=default, 
                    fault_tol=fault_tolerant)
            
        return rec
    
    def search_field_name_del(self, search_str, all_field_names=None, only_first=False, only_one=False, assert_found=False): 
        """
            #issue:  use AnalysisTools.search_str to implementate this
            avoid duplicated code 
        """
        #if isinstance(self, type): 
        #    all_field_names = self.ALL_FIELD_NAMES 
        #else: 
        #    all_field_names = self.iterkeys()
        all_field_names = all_field_names if all_field_names is not None else iter(self.keys())
        ss = search_str.split(' ')
        ss_len = len(ss)
        found = False
        res= []
        for fn in all_field_names: 
            match_count = 0
            for n in ss: 
                if n not in fn: 
                    break 
                else: 
                    match_count += 1
            if match_count == ss_len: 
                res.append(fn)
                found = True 
                if only_first: 
                    break 
        if only_one: 
            if len(res)>0: 
                res= [r for r in res if r == search_str]
        if assert_found: 
            assert len(res)>0, '"{}" matches no field'.format(search_str)
        return res 

    def search_field_name(self, name, has_sh=None, only_first=False, only_one=False, assert_found=False): 
        """
           
        """
        all_field_names =  iter(self.keys())
        res = self.search_str(name, all_field_names, only_first=only_first, 
                only_one = only_one, assert_found = assert_found,) 
        if has_sh is not None: 
            res= [r for r in res if isinstance(self[r], dict) and has_sh in self[r]]
        return res 

    def fetch_easier(self, name_str, sh, sub_key_list=None): 
        """
            status: 
                not completed 
                
            compare with fetch_easy, it allows fuzzy field_name, default sub_key_list, etc 
            the interface is easier than fetch_easy, but it is slower 
        """
        
        #raise NotImplemented  #at present search_field_name is only implemented in plot_field_vs_alpha 
        field_name = self.search_field_name(name_str, only_first=True, assert_found=True)[0]
        if sub_key_list is None: 
            if field_name in ['energy' ]: 
                pass 
            elif field_name in ['entanglement_entropy', 'entanglement']: 
                if self.ALGORITHM == 'mera' : 
                    sub_key_list = [0, 1]
                elif self.ALGORITHM == 'vmps' : 
                    sub_key_list = ['EE', sh[0]//2, 1] 
                else: 
                    raise NotImplemented
            else: 
                raise NotImplemented 
        res = self.fetch_easy(field_name, sh, sub_key_list)
        return res 
    
    def insert(self, field_name, sh, val, sub_key_list=None): 
        """
            remove param iter
        """
        
        if field_name not in self: 
            self[field_name] = OrderedDict()
        if sh not in self[field_name]: 
            self[field_name][sh] = OrderedDict()
            
        if sub_key_list is None: 
            self[field_name][sh] = val
        else:
            dic=self.__class__.make_nested_dict(
                    sub_key_list, old_dic=self[field_name][sh],val=val)
            self[field_name][sh] = dic 

    def rename_field(self, field_name_old, field_name_new): 
        print('rename key %s -> %s'%(field_name_old, field_name_new))
        self[field_name_new] = self[field_name_old]
        self.pop(field_name_old)
        self.commit()
    
    def rename_folder_surfix(self, new_surfix=None, dry_run=1): 
        """
            opeation on the parpath !! 
        """
        warnings.warn('still bug in some cases, so better dry_run first')
        
        db_root = os.path.dirname(self.parpath)
        state_root = os.path.dirname(self.state_parpath)
        
        old_name = os.path.basename(self.parpath)
        temp = old_name.split('-')[-1]
        
        if '=' in temp: 
            old_sur = ''
            new_name = '-'.join([old_name, new_surfix])
        else: 
            old_sur = temp
            new_name = old_name.replace(old_sur, new_surfix)
        
        db_parpath_new = '/'.join([db_root, new_name])
        state_parpath_new = '/'.join([state_root, new_name])

        if db_parpath_new[-1] == '-':
            db_parpath_new = db_parpath_new[:-1]

        if state_parpath_new[-1] == '-':
            state_parpath_new = state_parpath_new[:-1]
        
        print('rename ...%s -> ...%s'%(
                self.parpath[-100:], db_parpath_new[-100: ]))
        if not dry_run: 
            try:
                os.rename(self.parpath, db_parpath_new)
            except Exception as err:
                print(err)
        
        print('rename ...%s -> ...%s'%(
                self.state_parpath[-100:], state_parpath_new[-100: ]))
        if not dry_run: 
            try:
                os.rename(self.state_parpath, state_parpath_new)
            except Exception as err:
                print(err)
    
    def move_folder(self, local_root): 
        dir = self.parpath 
        cmd = 'mv %s %s/'%(dir, local_root)
        status=os.system(cmd)
        msg = 'mv ...%s -> ...%s/'%(dir[-20: ], local_root[-20: ])
        if status== 0:  
            msg += '   ... done' 
        else: 
            raise 
        print(msg) 
    
    def get_surfix(self): 
        pass 
        root = os.path.dirname(self.parpath)
        name = os.path.basename(self.parpath)
        temp = name.split('-')[-1]
        if '=' in temp: 
            sur = ''
        else: 
            sur = temp
        return sur 
    
    def mark_bad(self): 
        if 0: 
            sur = self.get_surfix()
            if sur  == '' : 
                new_sur = '-bad'
            else: 
                new_sur = '_'.join([sur, 'bad'])
            print_vars(vars(),  ['new_sur'])
            db.rename_folder(new_surfix=new_sur)

    def get_correaltion_exact(self, r, model_name=None, which='shastry'): 
        """
        """
        #r = r if r is not None else np.array([3**i for i in range(1, 6)])
        #if which == 'shastry': 
        if model_name in ['haldane_shastry', 'HS']: 
            #the exact value is: [(3, -0.044424436468189332), (9, -0.014200837573249929), (27, -0.0046643666206658257), (81, -0.0015470704943921049), (243, -0.00051483226149443202)]
            pi = np.pi
            res = 1./(1.0*pi*r)*sici(r*pi)[0]*(-1)**r
        elif model_name  == 'ising' : 
            #ref:  Bikas et.al. 1996,  Quantum Ising ... 
            #only support C^{xx}_{0, r} now 
            # res for 0<r<20:  [(1, 0.63661977236758138), (2, 0.54037964609246814), (3, 0.48926772236438937), (4, 0.4556470945476353), (5, 0.43107225309972047), (6, 0.41194225211020358), (7, 0.39641407232806986), (8, 0.38342749061536185), (9, 0.37232072902492203), (10, 0.36265500306854692), (11, 0.35412552042367873), (12, 0.34651258246402439), (13, 0.33965298163178659), (14, 0.33342240261115208), (15, 0.32772413252704091), (16, 0.32248156011136442), (17, 0.31763304115495111), (18, 0.3131282920367997), (19, 0.30892579903943168), (20, 0.3049909202249208)] 
            if 0:  #symblic result,  may be useful later 
                import sympy
                r=sympy.symbols('r')
                pi=sympy.pi 
                G=lambda r: 2/pi * (-1)**r/(2*r+1)
                sympy.pretty_print(G(r))
                G_mat=lambda n: sympy.Matrix(n,n, lambda i,j: G(j-i))
                G_mat(5)                
                
            def corr_ising_xx(n):
                G_func = lambda j: 2/np.pi * (-1)**j/(2*j+1)
                x, y = np.meshgrid(np.arange(n), np.arange(n))
                G_mat = G_func(x-y)
                res = np.linalg.det(G_mat)
                return res            
            if not hasattr(r, '__iter__'):  
                res = corr_ising_xx(r)
            else: 
                f_vec = np.vectorize(corr_ising_xx)
                res= f_vec(np.asarray(r))
        else: 
            raise 
        
        return res 
    
    def _get_Css(self, sh, **kwargs):
        db = self
        field = kwargs.get('field', 'correlation')
        field = kwargs.get('field', field)
        info = kwargs.get('info', 0)
        algorithm=db['algorithm']
        #field= 'correlation'
        if algorithm != 'mps':
            
            cx = db.fetch_easy(field, mera_shape=sh, sub_key_list=['xx'])
            cz = db.fetch_easy(field, mera_shape=sh, sub_key_list=['zz'])
            
            if cx is None or cz is None:
                return None
            rx, cx = list(zip(*list(cx.items()))); cx=np.asarray(cx)
            rz, cz = list(zip(*list(cz.items()))); cz=np.asarray(cz)
            assert rx==rz
            r=rx
            cxz= cx+cz
            res = OrderedDict(list(zip(r, cxz.tolist())))
        else:
            cx = db.fetch_easy(field, mera_shape=sh, sub_key_list=['xx'])
            cz = db.fetch_easy(field, mera_shape=sh, sub_key_list=['zz'])
            if cx is None or cz is None:
                if info>0: 
                    print('eiher cx or cz is None, so return None ')
                return None
            r0 = [i[0] for i in cx]
            r1 = [i[1] for i in cx]
            cx=[i[2] for i in cx]; cx=np.asarray(cx)
            cz=[i[2] for i in cz]; cz=np.asarray(cz)
            cxz= cx+cz
            res = list(zip(r0, r1, cxz.tolist()))
            
            pass
            
        return res  
 
    def _get_mag(self, sh, **kwargs):
        site = kwargs.get('site', sh[0]//2)
        mx = self.fetch_easy('magnetization', sh, ['x', site])
        mz = self.fetch_easy('magnetization', sh, ['z', site])
        try:
            return math.sqrt(abs(mx**2 + mz**2) )
        except TypeError:
            return None 
        except Exception:
            raise 

    def _get_EE_ln_from_log2(self, sh, i=None, info=0, **kwargs):
        """
            temporary work around 
        """
        assert i is not None 
        if i != 'all': 
            try:
                qn_spect_dic=self.fetch_easy('entanglement_entropy', sh, ['spectrum', i, 1])
                res = 0.0
                for v in list(qn_spect_dic.values()): 
                    res += -np.sum(v*np.log(v))
            except Exception as err:
                if info>0: 
                    warnings.warn(str(err))
                res=None
        else: 
            res=[(0, 0.0)]
            try:
                rec=self.fetch_easy('entanglement_entropy', sh, ['spectrum'])

                for i, dic in rec[1:]:
                    EE=0.0
                    for v in list(dic.values()): 
                        EE += -np.sum(v*np.log(v))
                    res.append((i, EE))
            except Exception as err:
                warnings.warn(str(err))
                res=None
        return res 

    def _get_concurence(self, sh, Jzz): 
        """
            this is only used for xxz model !
            refs: 
                Gu sj 2003 ent.. pra 
        """
        eng = self.fetch_easy('energy', sh)
        if self.ALGORITHM == 'mera' : 
            Czz = self.fetch_easy('correlation_extra', sh, sub_key_list=['zz'])
        else: 
            Czz = self.fetch_easy('correlation', sh, sub_key_list=['zz', 10, 0, 1])
        #print_vars(vars(),  ['Czz'])
        #raise
        res= 0.5*max(0.0, abs(eng-Jzz*Czz)-Czz-1.0)
        return res 
    
    def _get_corr_conn(self, sh, clip=0, spinfull=False, connected=True, fault_tol=True):
        """
            corr_conn(i, j) = <n_i n_j> - <n_i><n_j>
        """
        
        db = self
        N=sh[0]
        if not spinfull:
            corr = db.fetch_easy('correlation_all', sh, ['zz'])
            ni = db.fetch_easy('magnetization', sh, ['z'])
        else:
            corr = db.fetch_easy('correlation_all', sh, ['nn'])
            ni = db.fetch_easy('fermion_number', sh)
            
        if corr is None or ni is None:
            if not fault_tol:
                raise ValueError(
                "type(corr)=%s, type(ni)=%s"%(type(corr), type(ni))) 
            if corr is None:
                warnings.warn('corr is None')
            if ni is None:
                warnings.warn('ni is None')
            return None
        ni = list(ni.values())
        
        if 0<clip<1:
            clip = int(clip*N)
        L=N-2*clip
        if L<= 0: #clip too much
            return None
        corr_mat=corr
        if corr_mat.shape[0]<N:
            L = corr_mat.shape[0]
            start, end = (N - L)//2, (N + L)//2  
            ni = ni[start:end]

        ni=np.asarray(ni)
        if clip!=0:
            L=N-2*clip
            corr_mat=corr_mat[clip:-clip, clip:-clip]
            ni=ni[clip:-clip]
        
        ni=ni.reshape(L, 1)
        ni_mat=ni*ni.T  # it equals with ni.T*ni
        
        if not spinfull:
            if connected:
                corr_conn = 0.25*(corr_mat-ni_mat)
            else:
                corr_conn = 0.25*corr_mat
        else:
            corr_conn = corr - ni_mat
        return (corr_conn, L)

    def _get_corr_k(self, sh, k=None, clip=0, spinfull=False,  connected=True, per_site=False, fault_tol=1):
        """
            for spin models
            this formular returns the static structure factor
            S(k) = 1/N*sum(
                    exp(ik(i-j))*(<n_i n_j> - <n_i><n_j>))
        """
        temp = self._get_corr_conn(sh, clip=clip, spinfull=spinfull, connected=connected)
            
        if temp is None:
            if k is None:
                return None, None
            else:
                return None
        corr_conn_ij, L = temp 
        if L is None:
            raise 
        expij = lambda k, L: np.fromfunction(
                lambda i,j:np.exp(-1j*k*(i-j)), (L, L))
        
        func = lambda k: 1./L*np.sum(corr_conn_ij * expij(k, L))
        if k is None:
            return func, L 
        else:
            res = func(k)
            if per_site:
                res = res/L 
                #res = res/sh[0]
            return res 

    def _get_Luttinger_param(self, sh, clip=0, q_expr=None, 
            which=None, 
            spinfull=0, fault_tol=1):
        """
            params:
                q_expr: 
                    use this expresion to adjust the smallest momentum. 
                    I found q = 4*pi/N may be more accurate than 2*pi/N
                
            note  the fomular differs for spinfull and spinless models by a
            factor of 2.  
        """
        func, L = self._get_corr_k(sh, clip=clip, spinfull=spinfull, fault_tol=fault_tol) 
        pi = np.pi 
        if L is None:
            if not fault_tol:
                raise  
            return None 
        
        N = sh[0]
        assert L == N, (L, N) 
        q_expr = q_expr if q_expr is not None else '4*pi/N'
        q = eval(q_expr)  
        which = which if which is not None else 2
        if which == 1: # this choice is rather not accurate, I dont know why 
            res = pi*func(q)/q
            #print func(0)  #I checked,  this equal 0 quite well
        elif which == 2: 
            res = pi*(func(2*q)-func(q))/q
        if not spinfull :
            res *= 2
        return res
    
    #deprecated
    def _get_corr_k_fermion(self, sh, 
            x0=0, x1=1.1, clip=0,  num_points=100, 
            return_Sk = False, fault_tol=1, ):
        """
            for e.g. fermion hubbard model
        """
        db = self
        N=sh[0]
        L=N-2*clip
        corr=db.fetch_easy('correlation_all', sh, ['nn'])
        ni=db.fetch_easy('fermion_number', sh)
        if corr is None or ni is None:
            if not fault_tol:
                raise  
            return None,None
        ni=list(ni.values())
       
        if corr.shape[0]<N:
            L = corr.shape[0]
            start, end = (N - L)//2, (N + L)//2  
            ni = ni[start:end]
        if N==0:
            N=L=corr.shape[0]
            period=len(ni)
            ni=ni*(L/period)
        ni=np.asarray(ni)
        if clip!=0:
            L=N-2*clip
            corr=corr[clip:-clip, clip:-clip]
            ni=ni[clip:-clip]
        ni=ni.reshape(L,1)
        ninj=ni*ni.T # it equals with ni.T*ni
        corr_conn_ij = corr - ninj
        expij=  lambda k,  L: np.fromfunction(
               lambda i,j:np.exp(-1j*k*(i-j)), (L, L))
        
        if return_Sk:
            res= lambda k: np.sum(corr_conn_ij * expij(k, L))/L , L
            return res 
        x=np.linspace(x0, x1, num=num_points)
        y = [np.sum(corr_conn_ij * expij(np.pi*i, L))/L for i in x]    
        return x, y
    
    #deprecated
    def _get_Luttinger_param_fermion(self, sh, clip=0, fault_tol=1):
        #raise # deprecated use _get_Luttinger_param instead  
        pi = np.pi
        func, L = self._get_corr_k_fermion(sh, clip=clip, return_Sk=1, fault_tol=fault_tol) 
        if L is None:
            if not fault_tol:
                raise  
            return None 
        q = 2*pi/L
        #q = PI/L
        return pi*func(q)/q

    def _get_Luttinger_param_inf_len(self, sh=None, 
            x_min=40, x_max=10000, q_expr=None, which=None,  spinfull=0, fault_tol=1):
        db = self
        NN=db.get_shape_list(only_return_N=1, sh_min=(0, 1), from_energy_rec=1)
        data=[]
        for N in NN:
            rec=db._get_Luttinger_param((N, 'max'), 
                    q_expr=q_expr, which=which,  spinfull=spinfull)
            if rec is not None:
                data.append((1./N, rec))
        if not data:
            return None
        data.reverse()
        K=None
        try:
            func = lambda x, a0, a1, a2, a3: a0 + a1*x  + a2*x**2 +a3*x**3
            res=db.fit_line(None, func=func, data=data, 
                x_min=1./x_max, x_max=1./x_min, )
            K=res['param'][0] #except RuntimeError as err:
        except Exception as err: 
            warnings.warn(str(err))
        return K
   
    
    def _get_SF_k(self, sh, k=np.pi, per_site=True, fault_tol=1):  #delete this
        """
            
            calc magnetization by S(pi)
            attention:
                One question is do one need to calc correlation_all for full
                length?  
                The answer is: For a state that has good translation
                invariance, one dont need to, e.g. when Jzz=0 for the XXZ
                model; otherwise it is needed. 
                
        """
        raise #deprecated 
        func, L = self._get_corr_k(sh)
        if L is None:
            if not fault_tol:
                raise  
            return None 
        res = func(k)
        if per_site:
            res= res/sh[0]
        return res

    
    def _get_charge_gap(self, sh, N=1, z=0.0, fault_tol=True):
        p = self.__class__(self.parpath + '-Np%d'%N)
        h = self.__class__(self.parpath + '-Nm%d'%N)
        temp =  [x.fetch_easy('energy', sh, fault_tolerant=fault_tol) for x in [p, h, self] ]
        if all(temp):
            ep, eh, eg = temp
            #note eng is energy per bond
            N = sh[0]
            gap = (N-1)*(ep + eh - 2*eg)
            gap *= N**z 
        else:
            gap = None
        
        return gap 
    
    def _get_charge_gap_2(self, sh, surfix='Np1', z=0.0):
        p = self.__class__(self.parpath + '-' + surfix)
        temp = [x.fetch_easy('energy', sh) for x in [p, self] ]
        if all(temp):
            ep, eg = temp
            N = sh[0]
            #issue: 下面正好反了，
            if surfix in ['Np1', 'Np2']:  #vmps算法Np1其实是Nm1，这里将错就错
                gap = (ep - eg)
            else:
                gap = (eg - ep)
            gap = 2*(N-1)*gap #note eng is energy per bond
            
            gap *= N**z 
        else:
            gap = None
            
        return gap 
    
    def _get_chemical_potential(self, sh, surfix='Np1'):
        p = self.__class__(self.parpath + '-' + surfix)
        temp = [x.fetch_easy('energy', sh) for x in [p, self] ]
        if all(temp):
            e1, eg = temp
            #note eng is energy per bond
            N = sh[0]
            if surfix in ['Np1', 'Np2']:  
                gap = e1 - eg
            elif surfix in ['Nm1', 'Nm2']:
                gap = eg - e1
            else:
                raise  
            gap=(N-1)*gap
        else:
            gap = None
            
        return gap 
    
    _get_mu = _get_chemical_potential  #def _get_mu
    
    def _get_charge_gap_inf_len(db, sh=None, x_min=None, x_max=None):
        """
            sh is only a place holder
        """
        NN=db.get_shape_list(only_return_N=1, sh_min=(0, 1), from_energy_rec=1)
        data=[]
        for N in NN:
            rec=db._get_charge_gap_2((N, 'max'))
            if rec is not None:
                data.append((1./N, rec))
        if not data:
            return None #np.nan 
        data.reverse()
        K = None
        try:
            func = lambda x, a0, a1, a2, a3: a0 + a1*x  + a2*x**2 +a3*x**3
            res=db.fit_line(None, func=func, data=data, 
                        x_min=x_min, x_max=x_max, )
            K=res['param'][0] #except RuntimeError as err:
        except Exception as err: 
            warnings.warn(str(err))
        return K

    def _get_field_inf_length(self, sh, field_name, D=None, 
             x_min=20, x_max=10000, 
            fault_tol=True, args=None):
        """
            sh is only a place holder
            fitting by 3 order polynomial
        """
        mapper = {
                'Ck':self._get_corr_k, 
                'charge_gap':self._get_charge_gap,  
                'charge_gap_2':self._get_charge_gap_2,  
                'charge_stiff':self._get_charge_stiff, 
            }
        db = self 
        D = D if D is not None else  sh[1]
        NN=db.get_shape_list(only_return_N=1, sh_min=(0, 0), from_energy_rec=1)
        data=[]
        for N in NN:
            _sh = (N, D)
            if field_name not in mapper:
                rec = db.fetch_easy(field_name, _sh)
            else:
                args= args if args is not None else {}
                rec_getter = mapper[field_name].__func__
                rec = rec_getter(db, _sh, **args)
            if rec is not None:
                data.append((1./N, rec))
        if not data:
            #return np.nan 
            return None 
        data.reverse()
        try:
            func = lambda x, a0, a1, a2, a3: a0 + a1*x  + a2*x**2 +a3*x**3
            res=db.fit_line(None, func=func, data=data, 
                        x_min=1./x_max, x_max=1./x_min)
            K=res['param'][0] #except RuntimeError as err:
        except Exception as err: 
            if fault_tol:
                warnings.warn(str(err))
                #K=np.nan 
                K = None 
            else:
                raise  
        return K
    
    def _get_momentum_distribution(self, sh, rdm_one_fermion=None):  # def _get_nk 
        """
            Calc n(k) for 1D spinless fermions. Sum n(k) for k in the brllion
            zone shold equal tot particle number.
            
            note1: one can also use e.g. 
                kk = [-n*pace for n in range(L//2)]  + [n*pace for n in range(1, L//2  + 1)]
                It is simply a reordering.  Because they are identicle modulo 2*pi
        """
        if rdm_one_fermion is None:
            rdm_one_fermion = self.fetch_easy('rdm_one_fermion', sh)
        corr = rdm_one_fermion
        L = sh[0]
        exp = np.exp
        Ck = lambda k: 1/L*np.sum([exp(-1j*k*(x-y))*corr[x, y]  
            for x in range(L) for y in range(L)])
        pace = 2*np.pi/L
        kk = [n*pace for n in range(L)]  #note1
        kk.sort()
        kk = np.asarray(kk)
        
        nk = [Ck(k)  for k in kk]
        nk = np.asarray(nk)
        return kk, nk 
    
    
    def delete_rec(self, field_name_list, sh_list, dry_run=0, need_comfirm=0): 
        """
            params:
                field_name_list: 
                    it can take parameter 'all'
                sh_list: 
                    e.g. 'all', [(20, 'max')]
            pitfall:
                
        """
        msg = 'delete fields "%s" for sh in %s'%(field_name_list,  sh_list)
        if need_comfirm:
            i=input(msg + '?')
            if i.lower()!='y':
                print('canceled')
                return 
        else:
            print(msg)
        
        if field_name_list == 'all' : 
            field_name_list = list(self.keys())
        if sh_list == 'all' : 
            for i in field_name_list: 
                self.pop(i)
        else:
            for sh in sh_list: 
                if sh[1] == 'max' :
                    sh = sh[0], self.get_dim_max_for_N(sh[0]) 
                for f in field_name_list: 
                    temp = self[f]
                    if isinstance(temp, dict) and sh in temp: 
                        temp.pop(sh)
                    if not temp:  # all sh deleted 
                        self.pop(f)
            
            #the above wont affect dim_max, so change it here
            NN = [sh[0] for sh in sh_list]        
            for N in NN:
                D=self.get_dim_max_for_N(sh[0], force=1) 
                if D>0:
                    self['dim_max'][N] = D
                else:
                    if N in self['dim_max']:
                        self['dim_max'].pop(N)
        if not dry_run:       
            self.commit()
        print('done')

    def delete_db(self, info=1): 
        os.remove(self.path)
        if info>0: 
            msg = 'db is removed'
            print(msg)
    
    def delete_file(self, sh_list, delete_rec=True, need_comfirm=1, info=1):
        """
            params:
                sh_list: can be 'all'
        
        """
        if sh_list == 'all' : 
            sh_list = self.get_shape_list(from_energy_rec=False) 
            
        msg = 'delete files for sh in %s'%(sh_list)
        if need_comfirm:
            i=input(msg + '?')
            if i.lower()!='y':
                print('canceled')
                return 
        else:
            print(msg)
        
        #dir = self.parpath 
        dir = self.state_parpath
       
        for sh in sh_list: 
            fn=self.__class__.shape_to_backup_fn(sh)
            path = '/'.join([dir, fn])
            msg='rm file %s'%(path[-30: ])
            try: 
                os.remove(path)
                msg += '  ... done' 
            except Exception as err: 
                msg += '  ... failed. '  +  str(err)
            if info>0: 
                print(msg) 
        
        if delete_rec:
            NN = [sh[0] for sh in sh_list]
            sh_min = (min(NN), 0)
            sh_max = (max(NN), 10000)
            sh_list = self.get_shape_list(sh_min=sh_min, sh_max=sh_max, 
                    from_energy_rec=1)
            self.delete_rec('all', sh_list, need_comfirm=need_comfirm)
        
    def dump(self, file=None): 
        print('---'*30)
        print('all records in %s'%self.dbname)
        #print self.items()
        pprint.pprint(self, stream=sys.stdout)

    def merge(self, other, info=0): 
        """
            merge other into self
        """
        if self.get('verion') is not None: 
            raise # need modify 
        self.update(other)
        self.commit()
        if info>0: 
            msg = 'other rdb is merged into self'
            print(msg)
        
    def merge_remote(self, parpath, remote_host=None): 
        raise  # deprecated 
        remote_host = 'zhihuali@211.86.151.102:'
        #remote_root = 'current/run-long-better/'
        remote_path = remote_host + parpath   + '/RESULT.pickle.db'
        locale_parpath = self.parpath
        locale_path = locale_parpath + 'temp.pickle'
        cmd_str = 'scp %s %s'%(remote_path, locale_path)
        print(cmd_str)
        status=os.system(cmd_str)
        if status != 0: 
            msg = 'scp %s ---> %s'%(remote_path[-20:], locale_path[-20: ])                     
            print('scp failed')
        else: 
            other = ResultDB(self.parpath, 'temp.pickle')
            self.merge(other)
            os.system('rm %s'%locale_path)
            #os.system('ls %s'%self.parpath)
            print('merge succeed')
    
    def _plot(self, x, y, show_nan=False, **kwargs): 
        figsize = kwargs.get('figsize')
        figsize = (4, 3) if figsize is None else figsize
        fig = kwargs.get('fig', None)
        ax = kwargs.get('ax', None)
        fig_type = kwargs.get('fig_type')
        if x is None or y is None:
            warnings.warn('input array x or y is None')
            if kwargs.get('return_ax', False): 
                return fig, ax
            else: 
                return ax.figure
        x = np.asarray(x)
        y = np.asarray(y)
        if not show_nan:
            arg = ~np.isnan(y.astype(float))
            if not all(arg):
                x = x[arg]
                y = y[arg]
            if x.size<1:
                if kwargs.get('return_ax', False): 
                    return fig, ax
                else: 
                    return ax.figure
        
        if fig is None and ax is None: 
            fig=plt.figure(figsize=figsize)
        if ax is None: 
            if len(fig.axes)==0: 
                ax=fig.add_subplot(111)
            else: 
                ax_found = False
                #ax = fig.axes[-1]   # improved by below
                for aa in fig.axes: 
                    if len(aa.lines)==0: 
                        ax = aa
                        ax_found = True
                        break 
                if not ax_found: 
                    raise Exception('ax is not defined; examine the layout of fig')
        if 1:  
            style_dic = {
                    'marker':next(MARKER_CYCLE), 
                    'ms': 5,    #'markersize': 7, 
                    'mfc': 'None',  #'w',
                    'linestyle': None, 
                    #'dashes':(None, None), 
                    #'mec': 'r', 
                    'mew': 1,  
                    'label': '', 
                    #'color': None,   #color默认值不知到应该是什么， 才使其随机变化
                    'lw': 1,  
                    }    
            dic = style_dic.copy()
        
        temp = list(style_dic.keys())  + ['color']
        if kwargs.get('label') and kwargs.get('label_surfix'): 
            kwargs['label'] = '-'.join([str(kwargs['label']), str(kwargs['label_surfix'])])
        if kwargs.get('ls'):
            kwargs['linestyle'] = kwargs['ls']
            
        dic_temp = {i:kwargs.get(i) for i in temp if i in kwargs  }
        dic.update(dic_temp)
        linestyle = dic.get('linestyle')
        if linestyle and not isinstance(linestyle, str): # linestyle can be dashes tuple 
            dic.update(linestyle=None, dashes=linestyle)
        
        if kwargs.get('yfunc'): 
            #y = kwargs['yfunc'](y)
            yfunc = kwargs['yfunc']
            y = list(map(yfunc, y))
        if kwargs.get('xfunc'): 
            x = kwargs['xfunc'](x)
        if fig_type is None: 
            lines = ax.plot(x, y, **dic)
        elif fig_type == 'scatter' : 
            lines = ax.scatter(x, y, **dic)
        else: 
            raise 
        
        temp = ['xlabel', 'ylabel', 'xlim', 'ylim', 'xscale', 'yscale']
        for t in temp: 
            tt = kwargs.get(t)
            if tt is not None :
                if isinstance(tt, tuple): 
                    ax.__getattribute__('set_' + t)(*tt)
                else: 
                    ax.__getattribute__('set_' + t)(tt)
        #if kwargs.has_key('title') and ax.is_first_row(): 
        if 'title' in kwargs:   
            ax.set_title(kwargs.get('title')) 
            
        for l in lines:  #line stype   #matploblib 可能有个bug，mfc = 'w', 则 mec总是黑色，和line color 不同，故这里要重新set mec
            mec = kwargs.get('mec', l.get_color())
            l.set_mec(mec)
        #note the position to invoke .legend maters,  it must be placed after line.set_mec,  or else line.set_mec wont change mec in the legend
        ax.legend(title=kwargs.get('legend_title'), loc=kwargs.get('legend_loc', 0))
        #ax.grid(1)
        
        #ax.figure.set_facecolor('white')   # matplotrc not working, so set at here
        
        if kwargs.get('show_fig'): 
           plt.show() 
        if kwargs.get('return_ax', False): 
            return fig, ax
        else: 
            return ax.figure
    
    def plot_field_vs_dim(self, field_name, sh_list=None, sub_key_list=None, 
            inverse_D=False, sh_min=None, sh_max=None, xfunc=None, 
            yfunc=None, rec_getter=None, rec_getter_args=None,    **kwargs): 
        sub_key_list = sub_key_list if sub_key_list is not None else []
        info = kwargs.get('info', 0)
        sh_list = sh_list if sh_list is not None else self.get_shape_list(sh_max=sh_max, sh_min=sh_min)
        #aa=list(self.alpha_list)
        not_found=[]
        data=[]
        for sh in sh_list:  
            if rec_getter is None: 
                rec=self.fetch_easy(field_name, sh, sub_key_list=sub_key_list)
            else: 
                rec_getter_args=rec_getter_args if rec_getter_args is not None else {}
                if type(rec_getter)==types.MethodType:
                    rec_getter = rec_getter.__func__
                
                rec = rec_getter(self, sh, **rec_getter_args)
            dim = sh[1]
            if rec is not None:
                data.append((dim, rec))
            else:
                not_found.append(dim)
            
        if not_found and info: 
            print('not_found is ', not_found)
        #x, y=zip(*data)       
        data_is_empty = False 
        try: 
            x,y=list(zip(*data))        
            x = np.asarray(x)
            y = np.asarray(y)
            if xfunc is None: 
                if inverse_D:  #inverse_N takes effects provided xfunc is None
                    x = 1./x
            else: 
                x = xfunc(x)
        
        except ValueError as err: 
            print(err)
            #x, y = np.nan, np.nan 
            return 
        
        if 1:  # some per field specific treatment
            if field_name == 'energy' and kwargs.get('yscale') and yfunc is None: 
                y = y -np.min(y)
        x = np.asarray(x); y=np.asarray(y)
        if yfunc is not None : 
            y = yfunc(y)
        xlabel = '1/D' if inverse_D else 'D'
        kwargs.update(xlabel=xlabel, ylabel=field_name)
        fig = self._plot(x, y, **kwargs) 
        if kwargs.get('return_fig'): 
            return fig 
    
    def plot_field_vs_size(self, field_name, sh_list, sub_key_list=None, 
            data=None, inverse_N=True, xfunc=None, yfunc=None, 
            exponent=1, 
            rec_getter=None, rec_getter_args=None, **kwargs): 
        
        sub_key_list = sub_key_list if sub_key_list is not None else []
        info = kwargs.get('info', 0)
        #sh_list = sh_list if sh_list is not None else self.get_shape_list(sh_max=sh_max, sh_min=sh_min)
        #aa=list(self.alpha_list)
        not_found=[]
        if data is None:  
            data=[]
            for sh in sh_list:  
                if rec_getter is None: 
                    rec=self.fetch_easy(field_name, sh, sub_key_list=sub_key_list)
                else: 
                    rec_getter_args=rec_getter_args if rec_getter_args is not None else {}
                    if type(rec_getter)==types.MethodType:
                        rec_getter = rec_getter.__func__
                    rec = rec_getter(self, sh, **rec_getter_args)
                size = sh[0]
                if rec is not None:
                    data.append((size, rec))
                else:
                    not_found.append(size)
        if not data:
            warnings.warn('empty data')
        else:
            try: 
                x,y=list(zip(*data))        
                x = np.asarray(x)
                y = np.asarray(y)
                #y = y*x**exponent 
                if xfunc is None: 
                    if inverse_N:  #inverse_N takes effects provided xfunc is None
                        x = 1./x
                        x = x**exponent
                else: 
                    x = xfunc(x)
            except ValueError as err: 
                print(err)
                #x, y = np.nan, np.nan 
                return 
            
            if yfunc is not None : 
                y = yfunc(y)
            xlabel = '$1/N$' if inverse_N else '$N$'
            ylabel = field_name if rec_getter is None else (
                    rec_getter.__name__.replace('_get', ''))
            kwargs.update(xlabel=xlabel, ylabel=field_name, return_ax=1)
            fig, ax = self._plot(x, y, **kwargs) 
            
            #below may cause problem
            #tic= ax.get_xticks().tolist()
            #tic = ['%s\n%d'%(i, 1./(i)) for i in tic]
            #_ = ax.set_xticklabels(tic)    
        
        if kwargs.get('return_fig'): 
            return fig 
    
    def plot_field_vs_position(self, field_name, sh, key_list=None,
            r=None, r_min=None,  r_max=None,  period=None, rec=None,  fit=None, plot_fit=False, 
            fault_tolerant=1, rec_getter=None, rec_getter_args=None, 
            info=0, algorithm=None,  **kwargs): 
        rec_getter_args = rec_getter_args if rec_getter_args is not None else {}
        if rec is None: 
            if rec_getter is None: 
                rec = self.fetch_easy(field_name, mera_shape=sh, sub_key_list=key_list, 
                        fault_tolerant=fault_tolerant, info=info-1) 
            else: 
                if type(rec_getter)==types.MethodType:
                    rec_getter = rec_getter.__func__
                rec = rec_getter(self, sh, **rec_getter_args)
                field_name=rec_getter.__name__.replace('_get_', '').replace('_', ' ')  #for ylabel
                
        if rec is None: 
            return None
        algorithm = self['algorithm']
        #if algorithm in ['idmrg', 'mera']: 
        if algorithm  ==  'idmrg':  
            rec = list(rec.items())
            x, y = list(zip(*rec))
            if isinstance(x[0], tuple):
                x = [a[1] for a in x]        
                
        elif algorithm == 'mera' : 
            x, y = list(zip(*rec))
            
        elif 'mps' in algorithm:   #elif self['algorithm'] == 'vmps':
            if field_name in ['correlation', 'correlation_bond'] : 
                x, y = list(zip(*rec))
                r0 = key_list[-1]
                x = np.asarray(x)
                x -= r0 
            else: 
                if isinstance(rec, OrderedDict):
                    rec = list(rec.items())
                x, y = list(zip(*rec))
                
        elif algorithm == 'proj_qmc': 
            x, y = list(zip(*rec))
        else:
            x, y = list(zip(*rec))
        
        #if self.get('algorithm')=='idmrg' and 'corr' in field_name : 
        #    x = [a[1] for a in x]
        x = np.array(x); y=np.array(y)
        
        arg=self.filter_array(x, r_min, r_max, period)
        x = x[arg]
        y = y[arg]
        kwargs['return_ax'] = 1
        if 'xlabel' not in kwargs:
            kwargs['xlabel'] = '$i$'
        if 'ylabel' not in kwargs:
            kwargs['ylabel'] = field_name
        fig, ax=self._plot(x, y, **kwargs)
        if fit: 
            fit_succeed = 1
            if fit == 'power' : 
                func = lambda x, C, eta: -C*x**eta*(-1)**x
            else: 
                raise 
            try:
                param,cov = curve_fit(func, x, y)
                if isinstance(cov, float): 
                    cov = np.array([[cov]])
            except Exception as err: 
                print('fitting faild', err) 
                fit_succeed = 0
            if fit_succeed:     
                l = ax.lines[-1]
                label = l.get_label()
                label =  ''.join([label, ' %1.2f'%param[1]])
                l.set_label(label)
                
                
            if plot_fit: 
                color = l.get_color()
                if fit_succeed: 
                    #label = '' if not show_label else param[1]
                    dic = dict(color=color)
                    
                    fig=self._plot(x, func(x, *param), #xscale='log',yscale='log', 
                            return_fig=1, marker=None, ax=ax, **dic) #,  **kwargs)
        
        if kwargs.get('return_fig'): 
            return fig 
    
    plot_field_vs_length = plot_field_vs_position

    def plot_field_vs_k(self, field_name, sh,  key_list=None,
            rec=None, fault_tol=1, rec_getter=None, rec_getter_args=None, 
            info=0,  **kwargs): 
        raise NotImplemented
        rec_getter_args = rec_getter_args if rec_getter_args is not None else {}
        if rec is None: 
            if rec_getter is None: 
                rec = self.fetch_easy(field_name, mera_shape=sh, sub_key_list=key_list, 
                        fault_tolerant=fault_tol, info=info-1) 
            else: 
                if type(rec_getter)==types.MethodType:
                    rec_getter = rec_getter.__func__
                rec = rec_getter(self, sh, **rec_getter_args)
                field_name=rec_getter.__name__.replace('_get_', '').replace('_', ' ')  #for ylabel
                
        if rec is None: 
            return None
        algorithm = self['algorithm']
        #if algorithm in ['idmrg', 'mera']: 
        x, y = list(zip(*rec))
        
        #if self.get('algorithm')=='idmrg' and 'corr' in field_name : 
        #    x = [a[1] for a in x]
        x = np.array(x); y=np.array(y)
        
        kwargs['return_ax'] = 1
        if 'xlabel' not in kwargs:
            kwargs['xlabel'] = '$i$'
        if 'ylabel' not in kwargs:
            kwargs['ylabel'] = field_name
        fig, ax=self._plot(x, y, **kwargs)
        
        if kwargs.get('return_fig'): 
            return fig 
   
    def plot_field_vs_iter(self, attr_name, sh, data=None, x_min=None, x_max=None, period=None, sh_min=None, sh_max=None, 
            xfunc=None, yfunc=None,  sub_key_list=None, rec_getter=None, fault_tol=True,  rec_getter_args=None,    **kwargs): 
        db = self
        sub_key_list = sub_key_list if sub_key_list is not None else []
        
        info = kwargs.get('info', 0)
        #sh_list = sh_list if sh_list is not None else self.get_shape_list(sh_max=sh_max, sh_min=sh_min)
        #aa=list(self.alpha_list)
        
        try:
            if data is None:
                ts=db.get_time_serials(sh)#[:100]
                if hasattr(ts, 'N'):  # idmrg 
                    x=ts.N.values.astype(int)
                else: 
                    x =  ts.index.astype(int) 
                #y=ts.EE.values.astype(float)
                y=ts[attr_name].values.astype(float)
            else:
                x,y=list(zip(*data))        
                x = np.asarray(x)
                y = np.asarray(y)
            
            if x_min or x_max or period:
                arg = db.filter_array(x, x_min, x_max, period)
                x=x[arg]
                y=y[arg]
        except Exception as err:
            if fault_tol: 
                if info>0: 
                    print(err) 
                return 
            else: 
                print(sh, self.parpath) 
                raise 

        #x, y=zip(*data)       
        
        #if 1:  # some per field specific treatment
        #    if attr_name == 'energy' and kwargs.get('yscale') and yfunc is None: 
        #        y = y -np.min(y)
        #x = np.asarray(x); y=np.asarray(y)
        if yfunc is not None : 
            y = yfunc(y)
        if xfunc is not None: 
            x = xfunc(x)
        kwargs.update(xlabel='N', ylabel=attr_name)
        fig = self._plot(x, y, **kwargs) 
        if kwargs.get('return_fig'): 
            return fig 
    
    def plot_time_seirials(self, sh=None, attr=None, sh_list=None, **kwargs): 
        """
            attr: 
                'dt_real', 'dt_cpu', 
        """
        if sh is None: 
            sh_list = sh_list if sh_list is not None else self.get_shape_list()
        else: 
            sh_list = [sh]
        count = 0
        for sh in sh_list: 
            ts= self.load_S(sh)['time_serials']
            df = pd.DataFrame(ts).T 
            #print df.columns 
            temp = df[attr]
            x = list(temp.keys())
            y = temp.values 
            kwargs.update(return_ax=1, marker=None)
            fig, ax=self._plot(list(range(count, count+y.size)), y, **kwargs)
            kwargs['ax'] = ax 
            count += y.size  
    
    def plot_spectrum(db, field_name, sh, key_list=None, 
            num_level=20, rec_getter=None, rec_getter_args=None, fault_tol=True,  **kwargs):
        
        if rec is None: 
            if rec_getter is None: 
                rec = self.fetch_easy(field_name, mera_shape=sh, sub_key_list=key_list, 
                        fault_tolerant=fault_tolerant, info=info-1) 
            else: 
                if type(rec_getter)==types.MethodType:
                    rec_getter = rec_getter.__func__
                rec = rec_getter(self, sh, **rec_getter_args)
                field_name=rec_getter.__name__.replace('_get_', '').replace('_', ' ')  #for ylabel
        
        rec =db.fetch_easy(field_name, sh, key_list)
        if rec is None:
            warnings.warn('no record')
        else:
            n=num_level
            y=-np.log10(rec[:n])
            x=list(range(min(n, y.size)))
            kwargs.update(return_ax=1)
            kwargs['lw'] = kwargs.get('lw', 0)
            fig, ax=db._plot(x, y, 
                    **kwargs)    
    
            
    def plot_entanglement_scaling_old(self, alpha=None, sh_list=None, sh_min=None, sh_max=None, 
            log2=False, num_site=None, aver=False, linear_fit=True,  fit_points= None, 
            **kwargs): 
        #fig=None, show_res=0, return_fig=0, info=0): 
        fig = kwargs.get('fig')
        ax = kwargs.get('ax')
        show_res= kwargs.get('show_res')
        return_fig = kwargs.get('return_fig')
        info = kwargs.get('info', 0)
        if sh_list is None: 
            sh_list = self.get_mera_shape_list(sh_max=sh_max, sh_min=sh_min)
            if len(sh_list)==0: 
                raise
        if fig is None : 
            fig=plt.figure(figsize=(5,3))
            ax=fig.add_subplot(111)
         
        if len(fig.axes)==0: 
            ax=fig.add_subplot(111)
        else:
            for i in fig.axes: 
                if len(i.lines)==0:
                    ax = i
                    break
            #if ax is None: 
            #    ax = fig.axes[-1]
                
            
        not_found= []
        algorithm =  self.get('algorithm') 
        
        for sh in sh_list: 
            if self.get('algorithm')=='mps': 
                rec = self.fetch_easy('entanglement_entropy', sh, ['EE'], info=info-1)
                if rec is None: 
                    not_found.append(sh)
                    continue
                rec = rec[: len(rec)/2-1]
                x,y=(list(zip(*rec)))
                
            else: 
                name = 'entanglement_brute_force'
                if num_site is not None :   # in [6, 9]
                    #if num_site >= 9: 
                    #name += '_9' 
                    name += '_%d'%num_site 
                if aver: name += '_aver' 
                rec=self.fetch_easy(name, sh, ['EE'],)
                if rec is None: 
                    not_found.append(sh)
                    continue
                x,y=(list(zip(*list(rec.items()))))
            
            if show_res: 
                print(alpha, sh, rec)
            x=np.array(x); 
            y=np.array(y); 
            label = ''
            if alpha is not None : label += '%s'%(alpha, )
            if len(sh_list)>1: label += ' %s'%(sh, )
            
            if fit_points is None: 
                xx = x;  yy = y
            else: 
                xx = x[fit_points]
                yy = y[fit_points]
            
            if linear_fit: 
                c, b = np.polyfit(np.log2(xx.real[1: -1]), yy.real[1: -1], deg=1)
                c = c*3
                
            else: 
                func = lambda x, c, const: c/3.*np.log2(x) + const
                param,cov = curve_fit(func, xx.real[1: -1], yy.real[1: -1])
                #print alpha, sh, param.round(3),  cov.diagonal()
                c = param[0]
            
            label += ' %1.2f'%c 
            ax.plot(x,y, '.-', label=label)
            
        l=ax.legend(loc='lower right')    
        if rec is not None : 
            l.get_frame().set_alpha(0) # this will make the box totally transanalysis_class
            l.get_frame().set_edgecolor('white') # this will make the edges of the border white to match the background instead
        
        if log2 and rec is not None: 
            #ax.set_xticklabels(2**x)
            #tic = x.tolist()
            ax.set_xscale('log')
            if algorithm != 'mps':  
                tic = x[1: ]   # x[0] = 0, wont show in log scale
                ax.set_xticks(tic)
                ax.set_xticklabels(tic)
        ax.grid(1)
        
        if return_fig: 
            return fig
           
    def plot_entanglement_scaling(self, alpha=None, sh_list=None, sh_min=None, sh_max=None, 
            log2=False, num_site=None, site_start=0, site_period=1, site_end=None, aver=False, linear_fit=True,  fit_points= None, 
            **kwargs): 
       
        
        show_res= kwargs.get('show_res')
        return_fig = kwargs.get('return_fig')
        info = kwargs.get('info', 0)
        if sh_list is None: 
            sh_list = self.get_mera_shape_list(sh_max=sh_max, sh_min=sh_min)
            if len(sh_list)==0: 
                raise
        
        not_found = []
        algorithm =  self.get('algorithm') 
        
        for sh in sh_list: 
            if self.get('algorithm')=='mps': 
                rec = self.fetch_easy('entanglement_entropy', sh, ['EE'], info=info-1)
                if rec is None: 
                    not_found.append(sh)
                    continue
                #rec = rec[: len(rec)/2-1]
                site_end = len(rec)/2-1  if site_end is None else site_end
                rec = rec[site_start: site_end: site_period]
                x,y=(list(zip(*rec)))
                
            else: 
                name = 'entanglement_brute_force'
                if num_site is not None :   # in [6, 9]
                    #if num_site >= 9: 
                    #name += '_9' 
                    name += '_%d'%num_site 
                if aver: name += '_aver' 
                rec=self.fetch_easy(name, sh, ['EE'],)
                if rec is None: 
                    not_found.append(sh)
                    continue
                x,y=(list(zip(*list(rec.items()))))
            
            if show_res: 
                print(alpha, sh, rec)
            x=np.array(x); 
            y=np.array(y); 
            label = ''
            if alpha is not None : label += '%s'%(alpha, )
            if len(sh_list)>1: label += ' %s'%(sh, )
            
            if fit_points is None: 
                xx = x;  yy = y
            else: 
                xx = x[fit_points]
                yy = y[fit_points]
            
            if linear_fit: 
                c, b = np.polyfit(np.log2(xx.real[1: -1]), yy.real[1: -1], deg=1)
                c = c*3
                
            else: 
                func = lambda x, c, const: c/3.*np.log2(x) + const
                param,cov = curve_fit(func, xx.real[1: -1], yy.real[1: -1])
                #print alpha, sh, param.round(3),  cov.diagonal()
                c = param[0]
            
            label += ' %1.2f'%c 
            #ax.plot(x,y, '.-', label=label)
            kwargs['label'] = label
            kwargs.update(label=label, marker='.', return_ax=1)
            fig, ax=self._plot(x,y, **kwargs)
            kwargs['fig'] = fig
            kwargs['ax'] = ax
        
        if 'ax' in locals():          
            l=ax.legend(loc='lower right')    
            if rec is not None : 
                l.get_frame().set_alpha(0) # this will make the box totally transanalysis_class
                l.get_frame().set_edgecolor('white') # this will make the edges of the border white to match the background instead
            
            if log2 and rec is not None: 
                #ax.set_xticklabels(2**x)
                #tic = x.tolist()
                ax.set_xscale('log')
                if algorithm != 'mps':  
                    tic = x[1: ]   # x[0] = 0, wont show in log scale
                    ax.set_xticks(tic)
                    ax.set_xticklabels(tic)
            ax.grid(1)
        lll=locals()
        if return_fig and 'fig' in lll: 
            return fig
        if kwargs.get('return_ax'): 
            return lll.get('fig'), lll.get('ax')
       
    def plot_entanglement_scaling_compare(self, num_site_max=6): 
        """
            only for mera
        """
        for n, a in [(None, False),(6, False), (6, True), (9, False), (9,True)]:
            if n>num_site_max: 
                continue
            fig=self.plot_entanglement_scaling(num_site=n, sh_max=(12,6), aver=a, log2=1, return_fig=1)
            fig.axes[0].set_title('%s-%s'%(n,a))
    
    def plot_entanglement_vs_chi(self,  sh_list=None, sh_min=None, sh_max=None, filter_func=None,  **kwargs): 
        if sh_list is None: 
            sh_list = self.get_shape_list(sh_min=sh_min, sh_max=sh_max)
        #data = [(np.nan, np.nan)]
        data = []
        not_found = []
        for sh in sh_list: 
            ee = self.fetch_easy('entanglement', sh)
            if ee is None: 
                not_found.append(sh)
                continue
            data.append((sh[1], ee))
        if filter_func is not None: 
            data = list(filter(filter_func, data))
        if not_found: 
            print('rec not found for %s'% not_found)
        x, y = list(zip(*data))
        
        label=kwargs.get('label', '')
        x = np.log(x)
        try: 
            k, b= np.polyfit(x,y, deg=1)
            #print k, b 
            c=12*(1./k - 1)**(-2)
            label  += ' k=%1.2f c=%1.2f'% (k, c)
            kwargs['label'] = label 
        except Exception as err: 
            print(err)
       
        kwargs['return_ax'] = 1
        marker = kwargs.get('marker', '.')
        kwargs.update(marker=marker, xlabel='$\\ln(\\chi)$', ylabel='$EE$')
        fig, ax=self._plot(x, y,  **kwargs)
        ax.plot(x, k*x+b,)
        
        #ax.set_xticklabels()
        #tic = ax.get_xticklabels()
        #tic = [t.get_text() for t in tic]
        #print tic 
        tic= ax.get_xticks().tolist()
        tic = ['%s\n%d'%(i, math.exp(i)) for i in tic]
        ax.set_xticklabels(tic)


        #ax.set_xscale( 'log')
        #ax.set_xlim(x[0], x[-1])
        if kwargs.get('return_ax'): 
            return fig, ax

    def plot_structure_factor(self, sh, direct=None, dist_max=100, show_fig=0,  **kwargs): 
        fft = np.fft.fft
        algorithm = self.get('algorithm')
        if algorithm == 'mps': 
            name = 'correlation' 
            direct = 'xx' if direct is None else direct
            ind = 1
            key_list = [direct, 10]
        elif algorithm == 'idmrg' : 
            direct = 'xx' if direct is None else direct
            ind = 1
            key_list = [direct]
        elif algorithm == 'mera': 
            name = 'correlation_extra'
            direct = 'zz' if direct is None else direct
            ind = 1
            key_list = [direct]
        else: 
            raise 
        rec = self.fetch_easy(name, sh, sub_key_list=key_list)
        if rec is None: 
            msg =  'rec is None for sh %s, exit'%(sh, )
            warnings.warn(msg)
            return
        
        temp = [i[ind] for i in rec]
        temp = np.array(temp)[: dist_max]
        sf = fft(temp)
        
        fig = self._plot(list(range(len(sf))), sf, **kwargs)
        if kwargs.get('return_fig'): 
            return fig
        if show_fig: 
            plt.show()

    def plot_energy_diff_vs_time_old(self, sh, system_class=None, **kwargs):
        if self.get('version') is None: 
            if system_class is None: 
                from merapy.hamiltonian import System as system_class
            sh = sh[0], sh[1]-1
            path = self.parpath  +  system_class.shape_to_backup_fn(*sh)
            #print path; exit()
            res = system_class.load(path)
            recs = res.energy_record
        else: 
            sh = sh[0], sh[1]
            recs= self['energy'][sh]
        
        fig = kwargs.get('fig')
        if fig is None:
            fig = plt.figure()
        if len(fig.axes)==0: 
            ax = fig.add_subplot(111)
        else: 
            ax = fig.axes[-1]
        if 1:
            t_eng = [(t, v[0])  for t, v in recs.items()]
            t, eng = list(zip(*t_eng))
            t_diff = [(t_eng[i+1][0], (t_eng[i+1][1]-t_eng[i][1])/(t_eng[i+1][0]-t_eng[i][0])) 
                    for i in range(len(t_eng)-1)]
            t, diff = list(zip(*t_diff))
            diff = np.array(diff)
            diff_log = np.log10(abs(diff))  # 必须在log下画，否则什么都看不清楚
        ax.plot(t, np.sign(diff), label='$sign(E_{t+1}-E_t)$')           
        ax.plot(t, diff_log, label='$log(E_{t+1}-E_t)$')
        ax.set_xlabel('t',fontsize=16)
        l=ax.legend(prop={'size':14})  #fontsize
        if l is not None : 
            l.get_frame().set_alpha(0) # this will make the box totally transanalysis_class
            l.get_frame().set_edgecolor('white') # this will make the edges of the border white to match the background instead
       
        
        if kwargs.get('show_fig'): 
            plt.show()
    
    def plot_energy_vs_time(self, sh, switch=None, system_class=None, **kwargs):
        """
            kwargs: xlim, 
        """
        #print 'kkk', kwargs
        self._plot_energy_vs_time(sh, switch=switch, **kwargs)

    def update_energy_temp(self, sh): 
        """
            a func,  used temporaryly
        """
        aa = ['energy', tuple(sh)]
        if not self.has_key_list(aa): 
            if system_class is None: 
                from merapy.hamiltonian import System as system_class
            if len(sh)==2: 
                sh = sh[0], sh[1]-1
            else: 
                sh = sh[0], sh[1]-1, sh[2]
            path = self.parpath  +  system_class.shape_to_backup_fn(*sh)
            #print path; exit()
            res = system_class.load(path)
            recs = res.energy_record
            self.add_key_list(aa, val=recs, verbose=1)
            self.commit(info=1)

    def calc_fidelity(self, sh, other, force=0):
        algorithm = self['algorithm']
        if algorithm in ['vmps', 'mps'] : 
            from vmps.measure_and_analysis.measurement_vmps import fidelity
        elif algorithm == 'idmrg':
            from vmps.measure_and_analysis.measurement_idmrg_mcc import fidelity_direct  as fidelity
        else: 
            raise 
        which = 'fidelity'
        rec = self.fetch_easy(which, sh, sub_key_list=[other], fault_tolerant=1) 
        
        if rec is not None and not force:  
            return rec
        else:
            s1 = self.load_S(sh)
            #print_vars(vars(),  ['self.analysis_class'])
            db2 = self.analysis_class[other]
            #print_vars(vars(),  ['db2.parpath'])
            s2 = db2.load_S(sh)
            if s1 is None or s2 is None:
                res = None
            else:
                res = fidelity(s1, s2)        
                if 1:
                    ff = res 
                    a2 = other 
                    db1 = self 
                    N, D = sh
                    D1 = db1.get_dim_max_for_N(N=N) if D == 'max' else D 
                    D2 = db1.get_dim_max_for_N(N=N) if  D == 'max' else D 
                    #D = min(D1, D2)
                    db1.insert(which, (N, D1),  sub_key_list=[a2], val=ff)
                    db1.commit(info=0)
                    #db2.insert(which, (N, D),  sub_key_list=[a1], val=ff)
                
            
        return res 
    
    def count_eng_fluctuation(self, sh, iter_min=0, iter_max=1e10, info=0, **kwargs): 
        try: 
            recs = self.get_energy(sh)
        except IOError as err: 
            if info>0: 
                print(err)
            return None
        temp = np.asarray( [(k, v[0]) for k, v in recs.items() if iter_min< k < iter_max])
        iter = temp[: , 0].astype(int)
        eng = temp[: , 1]
        
        
        diff = eng[1: ] -eng[0: -1]
        
        ind_inc = np.where(diff>0)[0]
        
        num_pos = np.sign(diff) + 1 
        num_pos =np.count_nonzero(num_pos)
        if 0: 
            #print num_pos, len(ind_inc), type(ind_inc), ind_inc, len(ind_dec)
            print(ind_inc)
            print(iter[ind_inc])
            print(diff[ind_inc])
        mean_fluc = np.mean(diff[ind_inc[5:]])
        
        return num_pos, mean_fluc
    
    def _plot_energy_vs_time(self, sh, switch='diff', diff_lim=None, only_sign=False, sign_center=0, system_class=None, **kwargs):
        try:         
            recs = self.get_energy(sh)
        except Exception as err: 
            warnings.warn(str(err))
            return 
            
        if 1:
            t_eng = [(t, v[0])  for t, v in recs.items()]
            t, eng = list(zip(*t_eng))
        if 1:   
            if switch == 'diff':  
                t_diff = [(t_eng[i+1][0], (t_eng[i+1][1]-t_eng[i][1])/(t_eng[i+1][0]-t_eng[i][0])) 
                        for i in range(len(t_eng)-1)]
                t, diff = list(zip(*t_diff))
                t = np.array(t, dtype=int)
                diff = np.array(diff)
                if diff_lim is not None: 
                    arg = np.where(np.abs(diff)>=diff_lim)
                    t = t[arg]
                    diff = diff[arg]
                diff_log = np.log10(abs(diff))  # 必须在log下画，否则什么都看不清楚
                show_fig = kwargs.get('show_fig')
                kwargs['show_fig'] = 0
                if 'label' not in kwargs: 
                   kwargs['label']  = '$sign(E_{t+1}-E_t)$' 
                fig=self._plot(t, 0.8*np.sign(diff) + sign_center , **kwargs)           
                kwargs['show_fig'] = show_fig
                kwargs['fig'] = fig
                if not only_sign: 
                    kk = kwargs.copy()
                    if 'label' in kk: 
                        kk.pop('label')
                    self._plot(t, diff_log, label='$log(E_{t+1}-E_t)$', xlabel='t', **kk)
            elif switch  == 'err': 
                eng = np.array(eng)
                eng_min = np.min(eng)
                err = eng - eng_min
                #ax.plot(t, np.log10(err) , label='$E-E_{min}$')           
                self._plot(t, np.log10(err), xlabel='t', **kwargs)
  
    def show_scaling_dim(self, sh, alpha_remap={}, alpha_shape_dic={}, **kwargs): 
        data=[]
        if 1: 
            rec= self.fetch_easy('scaling_dim', sh)
            #print 'rrrr', rec
            if rec is not None:
                data.append(('a',) + rec[0][:5] + rec[1][:3] +rec[2][:3] )
                #data.append( rec[0][:5] + rec[1][:3] +rec[2][:3] )
            else: 
                print('rec not found for sh = %s'%(sh, ))
        print(tabulate(data, headers=['']+[str(i) for i in [0,0,0,0,0,1,1,1,2]], tablefmt='simple'))

    def upgrade_db(self): 
        if 0: 
            print(self.__version__); exit()
            if not hasattr(self, '__version__'): 
                self.__version__ = 0.9
        if self.get('version')==1.0: 
            print('already newest version')
            return 
        self.backup()
        new = ResultDB(parpath=self.parpath, dbname=self.dbname + '.new')
        
        if 1: 
            for sh_iter, V in self.items(): 
                sh = sh_iter[: -1]
                it = sh_iter[-1]
                for field_name, val in V.items(): 
                    if field_name not in new: 
                        new[field_name] = OrderedDict()
                    if sh not in new[field_name]: 
                        new[field_name][sh] = OrderedDict()
                    new[field_name][sh][it] = val
            
            #pprint.pprint(new['entanglement_entropy'][12, 5])

            new['version'] = 1.0
            new.commit(info=1)
            cmd = 'mv  %s %s'%(new.path, self.path)
            os.system(cmd)
            print('db is upgraded')
    
    def time_cost(self): 
        pass
        
    def backup(self): 
        #bac = ResultDB(self.parpath, dbname=self.dbname +'.bac')
        #bac.update(self)
        #bac.commit(info=1)
        
        #cmd = 'cp %s %s'%(self.path, self.path + '.bac')
        shutil.copyfile(self.path, self.path + '.bac')
        #os.system(cmd)
        print('db is backuped to RESULTDB.pickle.bac')

    def show_fig(self): 
        plt.show()
    
    def measure(self, sh_list, which=None, state=None, measure_func=None, 
            param=None, server=None,  field_surfix='', 
            force=0, fault_tolerant=1, priority=1,  num_of_threads=2, 
            use_local_storage=1, 
            job_group_name=None, 
            submit=1, info=1,  **kwargs): 
        """
            params:
                sh_list:
                    it can be 'all',  a tuple or a list of sh tuples
                which: can be str or a list of str
                    e.g. which='energy', which=['energy', 'variance']
        """
        from merapy.measure_and_analysis.measurement import measure_S 
        from mypy.brokest.task_center import submit_one, LOCAL_IP
        if isinstance(which, list) and len(which)==1:
            print_vars(vars(),  ['which'])
            which = which[0]
            print_vars(vars(),  ['wwwwww', 'which', 'len(which)'])
        if param is not None and isinstance(which, str):
            param = {which:param}  # to make it compatible with the interface of measure_S
        
        job_group_name = job_group_name if job_group_name is not None else str(which)
        job_group_name = '-'.join(['mea-', job_group_name])
        
        if sh_list == 'all' :
            sh_list = self.get_shape_list()
        if isinstance(sh_list, tuple):
            sh_list= [sh_list]
        for sh in sh_list:
            #if sh[1] == 'max' : 
            #    if '2site' in self.parpath:  #this is a temporary workaroud. alg info should be also stored in db in future 
            #        D = 0
            #    else:
            #        D = self.get_dim_max_for_N(sh[0])
            #        if D == 0: 
            #            print('max D for %s is zero, skip it'%(sh, ))
            #            continue
            #    sh = sh[0], D 
            if sh[1] == 'max':
                Dmax = self.get_dim_max_for_N(sh[0])
                sh = _sh = sh[0], Dmax
                if Dmax == 0 and not force:                
                    print('\tmax D for %s is zero, skip it'%(sh, ))
                    continue
                print('\tconverted %s->%s'%((sh[0], 'max'), _sh))
            elif sh[1] == 0:
                Dmax = self.get_dim_max_for_N(sh[0], force=force)
                if Dmax == 0 and not force:                
                    print('\tmax D for %s is zero, skip it'%(sh, ))
                    continue
                _sh = sh[0], self.get_dim_max_for_N(sh[0])
            else:
                _sh = sh
            if isinstance(which, str) and self.has_key_list([which, _sh]) and not force:
                print('   found, skip it')  
                continue
            
            ##################### make args  ############################
            parpath = self.parpath 
            path = self.shape_to_backup_path(sh)
            if 'C:' in parpath:
                xxx = [('C:', '/home/%s/dropbox/'%(LOCAL_USERNAME, )), ('Users', ''), ('zhihua', ''),  ('Dropbox', '')]
                for a, b in xxx:
                    parpath=parpath.replace(a, b)
                    path=path.replace(a, b)
                path=path.replace('dropbox', '')
               
            rdb = self 
            temp = ['parpath', 'path', 'rdb', 'which', 'measure_func', 'field_surfix', 'param', 
                    'force', 'fault_tolerant', 'num_of_threads']
            dic = vars()
            temp = {t: dic[t] for t in temp}
            temp.update(kwargs)
            temp['NUM_OF_THREADS'] = num_of_threads
            if not submit: 
                measure_S(**temp)
            else: 
                server = server if server is not None else (LOCAL_IP, 90999)
                temp['rdb'] = None # or else it raises 
                temp['fault_tolerant'] = 0
                temp.update(use_local_storage=use_local_storage)
                status = submit_one(measure_S, kwargs=temp,
                        server=server, 
                        job_info = {'priority': priority, 
                            'job_group_name': job_group_name, 
                            'job_description': (os.path.dirname(path), sh), 
                            'delay_send': 10, 
                            })
                
                print('\tresponse: ', status)

    def run(self):
        from vmps.minimize import VMPS 
        parpath = self.state_parpath
        #v = VMPS()
    
    @staticmethod
    def shape_to_backup_fn(sh): 
        N, D = sh
        res = 'N=%(N)d-D=%(D)d.pickle'%vars()
        return res
    
    def shape_to_backup_path(self, sh): 
        fn = self.__class__.shape_to_backup_fn(sh)
        dir = self.parpath.replace('Dropbox', '').replace('dropbox', '')
        dir = dir.replace(RESULTDB_DIR, BACKUP_STATE_DIR)
        return '/'.join([dir, fn])

    def update_db_structure(self, update_what='all'):
        """
            params:
                update_what:
                    update is not only controled by version number, 
                    also by a specific name 'update_what'
        """
        if self['version'] < 1.2 and update_what in ['all', 'remove_iter']:
            #this is big change, removed  -1 from the key list
            kk = list(self.keys())
            kk = [k for k in kk if k not in ['algorithm', 'status', 'version', 'dim_max']]
            ss = self.get_shape_list(from_energy_rec=1)  #here I use from_energy_rec because state file may be lost, while only db file retained
            for k in kk:
                for sh in ss:
                    if self.has_key_list([k, sh, -1]):
                        self[k][sh] = self[k][sh][-1]
            self['version']=1.2
            self.commit(info=1)    
            
class ResultDB_mera(ResultDB): 
    #AUTO_UPDATE = False 
    ALGORITHM = 'mera'
    ALL_FIELD_NAMES= list(ResultDB.ALL_FIELD_NAMES)
    ALL_FIELD_NAMES.extend([ 'entanglement_entropy', 'entanglement_special', 'correlation_extra'])
    def __init__(self, parpath,  **kwargs): 
        kwargs.update(algorithm='mera')
        ResultDB.__init__(self, parpath,  **kwargs)
    
    def parse_fn_bac(self, fn):
        return System.backup_fn_parse(fn)
    
    def parse_fn(self, fn):
        nqn = None
        if "transition" in fn or 'transtion' in fn: #'12-transition.pickle'
            dim = "0"
        elif "_" in fn:  #like "5_13.pickle or 5_13-4lay.pickle"
            dim = fn.split("-")[0].split("_")[1].split(".")[0]
            nqn = fn.split("-")[0].split("_")[0]
        else:
            if "lay" in fn: #"like "8-4lay.pickle""
                dim = fn.split("-")[0]
            else:   
                 #like "8.pickle"
                dim = fn.split(".")[0]
        
        i = fn.find("lay")
        layer = fn[i-1] if i>0 else "3"
        
        #return eval(dim), eval(layer) + 1  #  layer= mera.num_of_layer large than actual layer by 1 
        if nqn is None : 
            return eval(dim), eval(layer) + 1  #  layer= mera.num_of_layer large than actual layer by 1 
        else: 
            return eval(dim), eval(layer) + 1, eval(nqn)#  layer= mera.num_of_layer large than actual layer by 1 

    @staticmethod
    #def shape_to_backup_fn(sh,  nqn=None):
    def shape_to_backup_fn(sh ):  
        #layer, dim = sh
        dim, layer = sh[: 2]
        if len(sh)==3: 
            nqn = sh[2]
            nqn = "5_" if nqn == 5 else "" 
        else: 
            nqn = ""
        #layer = "-%dlay"%layer if layer != 3 else ""
        #layer = "-%dlay"%layer if layer != 4 else ""
        layer = "-%dlay"%(layer-1) if layer != 4 else ""
        res= "%(nqn)s%(dim)d%(layer)s.pickle"%vars()
        return res
    backup_fn_gen = shape_to_backup_fn 
    
    
    def get_correlation(self, sh, which='correlation_extra', direct='pm', 
            r_min=None,  r_max=None, period=None,  field_surfix=None, fault_tolerant=1, info=0): 
        algorithm = self.get('algorithm')
        if which == 'correlation_mixed' and self.get('algorithm') != 'mps': 
            rec1=self.fetch_easy('correlation_extra', mera_shape=sh, sub_key_list=[direct])
            rec2=self.fetch_easy('correlation', mera_shape=sh, sub_key_list=[direct, 'val'])
            if rec1 is None  or rec2 is None: 
                #not_found.append(sh)
                return None
            period  = period if period is not None else 2
            x1, y1=list(zip(*rec1[0:-1:period]))
            x2, y2 = list(zip(*rec2))
            x1 = np.array(x1);    x2 = np.array(x2);  
            y1 = np.array(y1);    y2 = np.array(y2);  
            x1_max = np.max(x1);  temp = np.where(x2>x1_max) 
            x = np.append(x1, x2[temp])
            y = np.append(y1, y2[temp])
        else: 
            if which == 'correlation': 
                key_list = [direct, 'val'] 
                period  = period if period is not None else 1                
            elif which == 'correlation_extra': 
                key_list = [direct] 
                period  = period if period is not None else 2            
            if self.get('algorithm') == 'mps': 
                which = 'correlation'  #not extra
                key_list = [direct]
                period = period if period is not None else 2            
            rec = self.fetch_easy(which, mera_shape=sh, sub_key_list=key_list, 
                    fault_tolerant=fault_tolerant, info=info-1) 
            if rec is None: 
                return None
            
            if algorithm == 'mps': 
                #rec = map(lambda x: (x[1]-x[0], x[2]), rec) 
                x, y = list(zip(*rec))
                r0 = key_list[-1]
                x = np.asarray(x)
                x -= r0 
                
            elif algorithm == 'idmrg' : 
                #rec = map(lambda x: (x[1]-x[0], rec[x]), rec.iteritems()) 
                rec = list(map(lambda k, v: (k[1]-k[0], v), iter(rec.items()))) 
                
            
            x,y=list(zip(*rec[0:-1:period]))
        
        if direct == 'pm':   
            y=np.array(y);  y=y*2;  
        
        x = np.array(x)
        a = 0
        b = -1
        if r_min is not None : 
            temp = np.where(x>=r_min);  a=temp[0][0]
        if r_max is not None : 
            temp = np.where(x<=r_max);  b=temp[0][-1] + 1
            x = x[a: b]
            y = y[a: b]
           
                   
        return x, y    

    def get_energy(self, sh, which='eng', iter_min=0, iter_max=1e6, system_class=None, **kwargs):
        aa = ['time_serials', tuple(sh)]
        iter_min = 0 if iter_min is None else iter_min
        iter_max = 1e6 if iter_max is None else iter_max
        if not self.has_key_list(aa, info=1): 
            
            #if len(sh)==2: 
            #    sh = sh[0], sh[1]-1
            #else: 
            #    sh = sh[0], sh[1]-1, sh[2]
            res = self.load_S(sh)
            if isinstance(res, dict): 
                recs= res['energy_record']
            else: 
                recs = res.energy_record
            #self.add_key_list(aa, val=recs, verbose=1)
            #self.commit(info=1)
        else:
            sh = sh[0], sh[1]
            recs = self['time_serials'][sh]
        if which in ['energy', 'eng']: 
            return recs
            
        elif which == 'diff': 
            pass
            temp = np.asarray( [(k, v[0]) for k, v in recs.items() if iter_min< k < iter_max])
            iter = temp[: , 0].astype(int)
            eng = temp[: , 1]
            diff = eng[1: ] -eng[0: -1]
            return diff
        else: 
            
            raise
    
    def fit_correlation(self, sh, func_type='power', which='correlation_extra', 
            direct='pm', r_min=None, r_max=None, period=None,  plot=0, plot_orig=1, show_label=1, **kwargs): 
        if 1: 
            power_func= lambda x, C, eta: -C*x**eta*(-1)**x
            #power_si_func= lambda x,  C,eta, m: -C*x**eta*sici(pi*m*x)[0]*(-1)**x
            power_si_func= lambda x,  C,eta, m: -C*x**eta*sici(m*x)[0]*(-1)**x
            power_cos_func= lambda x,  C,eta, m, b, k: -C*x**eta*(-1)**x * (b+ np.sin(1.*x/k))   #not usefull, only test
            def exp_func_1(t, a, b, c):                                                         
                return a * np.exp(b * t) + c   

            def exp_func_2(t, A, B):                                                         
                #return A * np.exp(-B * t) 
                return A * np.exp(- t/B) 

            if func_type == 'power': 
                func = power_func; 
            elif func_type == 'exp' : 
                func = exp_func_2
        if 1: 
            xy = self.get_correlation(sh, which, direct=direct, 
                    r_min=r_min, r_max=r_max,  period=period)
            if xy is None: 
                return None
            else: 
                x, y = xy
        
        if kwargs.get('info'): 
            print('xxx', x[0], x[-1], len(x))
        x=np.array(x); y=np.array(y);  
        y=np.abs(y)
        fit_succeed = 1
        try:
            param,cov = curve_fit(func, x, y)
            if isinstance(cov, float): 
                cov = np.array([[cov]])
            
        except Exception as err: 
            print('fitting faild', err) 
            fit_succeed = 0
      
        if fit_succeed: 
            if plot: 
                label = '' if not show_label else param[1]
                fig=self._plot(x, func(x, *param), xscale='log',yscale='log', 
                        label=label, return_fig=1, marker=None,  **kwargs)
                if plot_orig: 
                    if 'fig' not in kwargs: kwargs['fig'] = fig
                    fig=self._plot(x, y, xscale='log', yscale='log',   **kwargs)
            return param
       
        else: 
            return None
            
    def plot_correlation(self, sh,  direct='pm', r_min=None, r_max=None, period=None, which='correlation_mixed', 
            log_scale=1, fit=None, plot_fit=False,   return_fig=0, draw=0, **kwargs): 
        
        if 1: 
            assert which in ['correlation_extra', 'correlation', 'correlation_mixed']
            sh_list = sh
            fig = kwargs.get('fig')
            ax = kwargs.get('ax')
            if ax is None: 
                if fig is None :
                    fig=plt.figure(figsize=(6, 4))
                if len(fig.axes)==0: 
                    ax=fig.add_subplot(111)
                else: 
                    ax = fig.axes[-1]
        not_found=[]
        is_many_sh = 0
        if len(sh_list)==2 : # and isinstance(sh_list[0], int): 
            sh_list = [sh_list]
            #is_many_sh = 1
        #direct_orig = direct
        #if direct == 'pm': direct = 'xx' 
        sub_key_list = [direct] if which == 'correlation_extra' else [direct, 'val'] 
        if period is None: 
            period = 2 if log_scale else 1
        for sh in sh_list: 
                
            rec = self.get_correlation(sh, which, direct, r_min=r_min, r_max=r_max, period=period )
            if rec is not None: 
                x, y = rec
            else: 
                continue
            
            if direct == 'pm':   
                y=np.array(y);  y=y*2;  

            label = kwargs.get('label')
            
            label = label if label is not None else ''
            #label = ''  #if len(sh_list)<2: label += ' %s  '%(a, )   #if len(alpha_list)<2: label += ' %s  '%(sh, )
            if is_many_sh: label += ' %s  '%(sh, )
            marker = kwargs.get('marker')
            marker = marker if marker is not None else '--'
            if fit is None and len(marker)==1: marker  += '-'
            
            if log_scale: 
                y=np.abs(y)
                #ax.loglog(x, y, marker, label=label)
                ax.loglog(x, y, marker)
            else: 
                ax.plot(x, y, marker)
            
            line = ax.lines[-1]
            color = line.get_color()
            line.set_markersize(5)
            line.set_mfc('w')
            line.set_markeredgewidth(1)
            line.set_mec(color)
            if fit is not None : 
                    #line.set_mfc('w')
                    #line.set_markersize(5)
                res=self.fit_correlation(sh, func_type=fit,  which=which, direct=direct, 
                        r_min=r_min, r_max=r_max, period=period,  plot=plot_fit, plot_orig=0, 
                        show_label=0, ax=ax, fig=fig,  color=color)
                if res is not None : 
                    #if fit == 'power': 
                    A, eta  = res
                    if fit == 'power':  
                        name = '\\eta' 
                    else: 
                        name = '\\xi'
                    label += ' $%s=%1.2f$'%(name, eta)
                       
            line.set_label(label)
            if fit: ax.set_xlim(x[0], x[-1]*10)
        l=ax.legend()   #l=ax.legend(title='$\\alpha$')
        if l is not None : 
            l.get_frame().set_alpha(0) # this will make the box totally transanalysis_class
            l.get_frame().set_edgecolor('white') # this will make the edges of the border white to match the background instead
           
        ax.grid(1)
        if kwargs.get('show_fig'): plt.show()
        if draw: plt.show()
        if len(not_found)>0: print('not_found is ', not_found)
        if return_fig: return fig
    
    def plot_central_charge_vs_layer(self, sh, alpha_remap={}, alpha_sh_map={}, layer='all',  **kwargs):
        if layer=='top':
            ll= [sh[1]-2 ] 
        elif layer=='all':
            ll=  list(range(sh[1]-1))
        else: 
            ll = layer
        
        #ll= [sh[1]-2 ] if layer is None else range(sh[1]-1)
        not_found = []
        cc=[]
        for layer in ll:
            if 1: 
                #rec = db.fetch_easy('entanglement_entropy', shape, [layer])
                rec = self.fetch_easy('entanglement_entropy', sh)
                if rec is not None : 
                    rec = rec.get(layer)
                if rec is None: 
                    #not_found.append((a, layer))
                    not_found.append( layer)
                    continue
                #print 'rrrr', rec
                e1, e2 = rec[1], rec[2]
                #print a, e2, e1
                c= (e2-e1)*3
                cc.append((layer,c))
        #if len(cc)==0: 
        #    fig, ax = None, None
        #else: 
        if 1: 
            x,y=list(zip(*cc))    
            kwargs.update(marker='.', return_ax=1)
            fig, ax = self._plot(x, y, **kwargs)
            #ax.plot(x,y, '.-')    
            ax.legend()
            ax.set_xticks(ll)
            #ax.invert_xaxis()
            ax.grid(1)
        if len(not_found)>0: 
            print('rec for central_charge not_found is :%s'%not_found)
        if kwargs.get('return_ax'): 
            return fig, ax
        if kwargs.get('return_fig'): 
            return fig
        if kwargs.get('show_fig'): 
            plt.show()

class ResultDB_vmps(ResultDB): 
    VERSION = 1.2
    ALGORITHM = 'vmps'
    def __init__(self, parpath,  **kwargs): 
        #kwargs['version'] = 1.0    #version should det
        kwargs.update(algorithm='mps')
        
        ResultDB.__init__(self, parpath,  **kwargs)
        #super(ResultDB_vmps, self).__init__(parpath,  **kwargs)
    
    @staticmethod
    def shape_to_backup_fn(sh): 
        N, D = sh
        res = 'N=%(N)d-D=%(D)d.pickle'%vars()
        return res
    backup_fn_gen = shape_to_backup_fn 
    
    def update_db_structure(self): 
        if self['version'] < 1.01:   
            print('update version from %s to %s'%(self['version'], 1.01))
            def update_correlation(db):
                if db['version']==1.0:
                    ss= db.get_shape_list()
                    for sh in ss: 
                        rec=db.fetch_easy_v1p0('correlation', sh)
                        if rec is None: 
                            continue 
                        #db.insert('correlation_bac_v1.0', sh, -1, rec)
                        #db.insert('correlation_bac_v1.0', sh, rec)
                        rec_new=OrderedDict()
                        for dir in list(rec.keys()):
                            temp=rec[dir]
                            rec_new[dir]=OrderedDict()
                            x0, x1,y=list(zip(*temp))
                            r=list(zip(x0, x1))
                            r=[x[1]-x[0] for x in r]
                            
                            r0=x0[0]
                            rec_new[dir][r0]=list(zip(x1, y))

                        #db.insert('correlation', sh, -1, rec_new)
                        db.insert('correlation', sh, rec_new)
                    db['version']=1.01
                    db.commit(info=1)    
            update_correlation(self)
        
        if self['version'] <  1.02: 
            print('update version from %s to %s'%(self['version'], 1.02))
            sh_list = self.get_shape_list()
            for sh in sh_list: 
                N = sh[0]
                self.get_dim_max_for_N(N)  
            self['version'] = 1.02
            self.commit(info=1)
        
        if self['version']<1.2:
            ResultDB.update_db_structure(self, 'remove_iter')
        
    def calc_correlation(self, sh, direct, i0, j_list=None, force=False, info=0):
        """
            this is a temporary workaround 
        """
        
        from vmps.measure_and_analysis.measurement_vmps import correlation  as func 
       
      
        if j_list is None:
            j_list = list(range(i0+1, sh[0]))
            
        try: 
            rec = self['correlation'][sh][-1][direct][i0]
        except KeyError: 
            self.add_key_list(['correlation', sh, -1, direct, i0])
            #rec = self['correlation'][sh][-1][direct][i0]
            #print_vars(vars(),  ['rec'])
            #raise
            rec = []
        
        j_list_old = [a[0] for a in rec]
        if not force: 
            j_list = set(j_list)-set(j_list_old)
        ij_list = [(i0, j) for j in j_list]
        print_vars(vars(),  ['j_list_old', 'ij_list'])
        if len(j_list)>0: 
            state = self.load_S(sh)
            temp=func(state, [direct], ij_list=ij_list)
            print_vars(vars(),  ['temp'])
            rec.extend(temp[direct][i0])
            rec.sort(key=lambda x:x[0])
            self['correlation'][sh][-1][direct][i0] = rec
            self.commit(info=info)
    
    def calc_base(self, sh, func, key_list):
        pass
    
    def calc_current(self, sh, force=False):
        from vmps.mpo import MPO_Fatcory 
        key_list = ['current', sh]
        N = sh[0]
        if sh[1] == 'max' :
            D = self.get_dim_max_for_N(N)
            _sh = (N, D)
            key_list[1] = _sh
        else:
            _sh = sh
            
        res = self.fetch_from_key_list(key_list, default=None)
        if res is None or force:
            op = MPO_Fatcory.current_operator(N)
            state = self.load_S(sh)
            if state is not None: 
                mps= state['mps']
                res = mps.expectationvalue(op)
                self.insert_key_list(key_list, val=res)
                self.commit()
    
    def _get_energy_better(self, sh, N_diff=2):
        """
            better energy per bond. this solves the problem caused by OBC of VMPS
        """
        N = sh[0]
        e  = self.fetch_easy('energy', sh)
        e2 = self.fetch_easy('energy', (sh[0]+N_diff, sh[1]))
        #print_vars(vars(),  ['e, e2'])
        if e is not None  and e2 is not None:
            res=(e2*(N-1+N_diff)-e*(N-1))/N_diff
        else:
            res= None
        return res
     
    def _get_corr_len_bac(self, sh, k0=np.pi, z=1.0, connected=1, per_site=1):
        """
            extract correlation length from structure factor
            it can be used to locate QCP by scaling of \\xi 
            ref: Sandvik 2010 p.38 eq.70
        """
        L = sh[0]
        
        s0 = self._get_corr_k(sh, k=k0, connected=connected, per_site=0) 
        if s0 is None:
            return None
        q1 = k0 + 2*pi/L
        s1 = self._get_corr_k(sh, k=q1, connected=connected,  per_site=0) 
        res = 1/q1 * np.sqrt(s0/s1 - 1)  #note this res is already per site
        if not per_site:
            res *= L
        return res
    
    def _get_corr_len(self, sh, k0=np.pi, z=1.0, connected=1, per_site=1):
        """
            intro:
                extract correlation length from structure factor
                it can be used to locate QCP by scaling of \\xi 
                ref: 
                    Sandvik 2010 p.38 eq.70
            params:
                z: an critical exponent
        """
        L = sh[0]
        q1 = k0 + 2*pi/L
        s0 = self._get_corr_k(sh, k=k0, connected=connected, per_site=0) 
        if s0 is None:
            return None
        s1 = self._get_corr_k(sh, k=q1, connected=connected,  per_site=0) 
        res = 1/q1 * np.sqrt(s0/s1 - 1)  #note this res is already per site
        #if not per_site:
        #    res *= L**z
        #return res
        
        res = res*L**(z-1)
    
        return res
        
    
    def _get_central_charge(self, sh, fix_N=True, x_min=40, x_max=60, period=2):
        """
            
        """
        db=self
        if fix_N:
            data = self.fetch_easy('entanglement_entropy', sh, ['EE'])
        else :
            NN=db.get_shape_list(only_return_N=1)
            data=[(N, db.fetch_easy('entanglement_entropy', (N,'max'), 
                        sub_key_list=['EE', N//2, 1])) for N in NN]
            period=None
            data = [i for i in data if i[1] is not None ]
        if not data:
            return None 
        
        func = lambda x, c, b: c/6.*np.log2(x) + b
        try:
            res=self.fit_line(None, func=func, data=data, x_min=x_min, 
                    x_max=x_max, period=period)
            c=res['param'][0]
        except RuntimeError as err:
            warnings.warn(str(err))
            c=None
        except TypeError as err:
            warnings.warn(str(err))
            c=None
        return c 
    _get_c = _get_central_charge     #def _get_c
    
    def _get_charge_stiff(self, sh, dphi=0.05, fault_tol=True):
        """
            here assumes num of fermion to be odd,  
            so that E(phi) has minimu at phi=0.0  
            note1:
                there are different definations of D up to a factor of pi,
                here I dont multiply it by pi 
                
            ref:
                Giamarch's book p.405 and p.226 eq.7.80  
            examine:
                for free spinless fermion 
                D = 2/pi        at half filling
                D = sqrt(3)/pi  at 1/3 filling 
        
        """
        path = self.parpath.replace('phi=0.0', 'phi=%s'%(dphi))
        #print_vars(vars(),  ['path'])
        db_phi = self.__class__(path)
        temp =  [x.fetch_easy('energy', sh, fault_tolerant=fault_tol) 
                for x in [self, db_phi] ]
        if all(temp):
            N = sh[0]
            E_0,  E_phi = temp
            E_0 = (N-1)*E_0   #note eng is energy per bond
            E_phi = (N-1)*E_phi 
            #D = pi* N* 2* (E_phi - E_0)/ dphi**2  
            D =  N* 2* (E_phi - E_0)/ dphi**2   #note1
        else:
            D = None
        
        return D
    
    def _get_charge_stiff_2(self, sh, fault_tol=True):
        """
            phase sensitivity 
        """
        db_phi = self.__class__(self.parpath + '-anti_pbc')
        temp =  [x.fetch_easy('energy', sh, fault_tolerant=fault_tol) 
                for x in [self, db_phi] ]
        if all(temp):
            N = sh[0]
            E_0,  E_phi = temp
            #note eng is energy per bond
            res =(-1)**N *N*(N-1)*(E_0 - E_phi)
            res= res/pi*2
        else:
            res = None
        
        return res 
    
    def _get_current(self, sh, dphi=0.05, times_N=True,  fault_tol=True):
        """
        
        """
        temp = self.parpath.split('-')
        temp = [i for  i in temp if 'phi' in i][0]
        phi_str = temp.split('=')[1]
        phi = eval(phi_str.replace('m', '-'))
        
        phi2 = phi + dphi
        if 0:
            n = int(phi/dphi)
            phi2 = (n + 2)*dphi
        if 0:
            print_vars(vars(),  ['phi2/dphi'])
            phi2 = round(phi2/dphi, 14)*dphi 
        phi2 = round(phi2, 6)
        phi2_str = str(phi2).replace('-', 'm')
       
        path = self.parpath.replace(phi_str, phi2_str)
        db2 = self.__class__(path)
        temp =  [x.fetch_easy('energy', sh, fault_tolerant=fault_tol) 
                for x in [self, db2] ]
        if all(temp):
            N = sh[0]
            E_0,  E_phi = temp
            E_0 = (N-1)*E_0   #note eng is energy per bond
            E_phi = (N-1)*E_phi 
            I = -(E_phi-E_0)/dphi 
            if times_N:
                I *= N
        else:
            I = None 
        
        return I

class ResultDB_idmrg(ResultDB): 
    VERSION = 1.21
    ALGORITHM = 'idmrg'
    ALL_FIELD_NAMES= list(ResultDB.ALL_FIELD_NAMES)
    ALL_FIELD_NAMES.extend(['length', 'entanglement', 'correlation_length'])
    ALL_FIELD_NAMES= list(set(ALL_FIELD_NAMES))
    
    def __init__(self, parpath,  **kwargs): 
        #kwargs['version'] = 1.0
        kwargs.update(algorithm='idmrg')
        #ResultDB.__init__(self, parpath,  **kwargs)
        super(ResultDB_idmrg, self).__init__(parpath,  **kwargs)
    
    @staticmethod 
    def shape_to_backup_fn(sh): 
        """
        """
        N, D = sh
        res = 'N=%(N)d-D=%(D)d.pickle'%vars()
        return res
    
    backup_fn_gen =shape_to_backup_fn 
    
   
    def _measure(self, func, dim):
        res = OrderedDict()
        try:
            path='/'.join([self.parpath, 'N=0-D=%d.pickle'%dim])
            #S=iDMRG.load( path)
            S= self.load_S(path)
            res = func.__func__(None, S)
            
        except:
            #print a, path
            raise
            
        return res

    def plot_entanglement_normalized_vs_chi(self, sh_min=None, sh_max=None, filter_func=None,  **kwargs): 
        sh_list = self.get_shape_list(sh_min=sh_min, sh_max=sh_max)
        """
            only for test use 
        """
        #data = [(np.nan, np.nan)]
        data = []
        not_found = []
        for sh in sh_list: 
            ee = self.fetch_easy('entanglement_normalized', sh)
            if ee is None: 
                not_found.append(sh)
                continue
            data.append((sh[1], ee))
        if filter_func is not None: 
            data = list(filter(filter_func, data))
        if not_found: 
            print('rec not found for %s'% not_found)
        x, y = list(zip(*data))
        
        label=kwargs.get('label', '')
        x = np.log(x)
        try: 
            k, b= np.polyfit(x,y, deg=1)
            #print k, b 
            c=12*(1./k - 1)**(-2)
            label  += ' k=%1.2f c=%1.2f'% (k, c)
            kwargs['label'] = label 
        except Exception as err: 
            print(err)
       
        kwargs['return_ax'] = 1
        marker = kwargs.get('marker', '.')
        kwargs.update(marker=marker, xlabel='$\\ln(\\chi)$', ylabel='$EE$')
        fig, ax=self._plot(x, y,  **kwargs)
        ax.plot(x, k*x+b,)


        #ax.set_xscale( 'log')
        #ax.set_xlim(x[0], x[-1])
        if kwargs.get('return_ax'): 
            return fig, ax
    
    def calc_EE_from_finite_mps(self, sh, small_size=0, force=0, which_log='ln', fault_tol=1, info=0):
        field = 'EE_from_finite_mps'
        field = '_'.join([field, which_log])
        if small_size: 
            sh_finite = small_size, sh[1]
        else: 
            sh_finite = sh
        
        #rec=self.fetch_easy(field, sh_finite, fault_tolerant=fault_tol)
        rec=self.fetch_easy(field, sh_finite)
        
        if rec is None or force:
            try:
                state=self.load_S(sh, to_obj=1)
                from vmps.mps import MPS
                from vmps.idmrg_mcc import iDMRG_mcc
                state.N_max=1
                temp=iDMRG_mcc.generate_finite_mps_func.__func__(state)
                
                size = len(temp)
                if state.combine_2site: 
                    size = size*2 
                if small_size: 
                    if 1: #注意这个函数仅仅是临时使用, 参数都很特殊 
                        state.is_calc_Heff = 0 
                        state.site_dim = state.Wi.shape[2].copy()
                        state.which_eig = 'eigsh'
                        state.eigs_tol = 1e-12 
                        state.eigs_failure_lim = 0 
                        state.is_calc_v0_overlap = 0
                        state.trunc_err_tol = 1e-9 
                    temp = iDMRG_mcc.reduce_finite_mps_size.__func__(state, small_size, temp)
                    
                    
                    
                mps = temp.to_normal_form('right', inplace=0)
                #assert abs(abs(MPS.overlap(temp, mps)) -1.0 ) <= 1e-13
                if 1:
                    mmm = []
                    d=state.mpo[1].shape[2]
                    for A in mps: 
                        a1, a2 = MPS.split_2sites(A, 'right', d=d) 
                        mmm.extend([a1, a2])
                    mps = MPS(temp.N*2, None, d, mps=mmm, symmetry=state.symmetry) 
                import vmps.measure_and_analysis.measurement_vmps as mv
                if which_log == 'ln':  
                    rec=mv.entanglement_entropy_ln({'mps':mps})
                elif which_log == 'log2': 
                    rec=mv.entanglement_entropy({'mps':mps})
                else: 
                    raise ValueError(which_log)
                    
                rec = OrderedDict(rec['EE'])
                #self.insert(field, sh_finite, -1, rec)
                self.insert(field, sh_finite, rec)
                self.commit(info=info)
            except:
                if not fault_tol:
                    raise 
        #print_vars(vars(),  ['rec'])
        return rec

    def update_db_structure(self): 
        #if not self.has_key('version')
        if self['version'] <  1.01: 
            ver = 1.01
            print('update version from %s to %s'%(self['version'], ver))
            sh_list = self.get_shape_list()
            self.get_dim_max_for_N(0, update_db=1)  
            self['version'] = ver 
            self.commit(info=1)
        
        if self['version']<1.2:
            ver = 1.2
            print('update version from %s to %s'%(self['version'], ver))
            try:
                ResultDB.update_db_structure(self, 'remove_iter')
                succeed = True
            except Exception as err:
                succeed = False
                msg = 'failed to update, err is {}'.format(str(err))
                warnings.warn(msg)
                raise  
            if succeed:
                self['version'] = ver 
                self.commit(info=1)
            
        
        if self['version']<1.21:
            ver = 1.21
            print('update version from %s to %s'%(self['version'], ver))
            def update_correlation(db):
                #ss= db.get_shape_list()
                ss = list(db.get('correlation', {}).keys())
                for sh in ss: 
                    rec=db.fetch_easy('correlation', sh)
                    if rec is None: 
                        continue 
                    rec_new = OrderedDict()
                    for dir in list(rec.keys()):
                        if not isinstance(rec[dir], OrderedDict): 
                            continue
                        old = rec[dir]
                        rec_new[dir] = OrderedDict()
                        rr = list(old.keys())
                        rr = [r[1] for r in rr]
                        temp = [(r, old[(0, r)]) for r in rr]
                        rec_new[dir][0] = OrderedDict(temp)
                        #print_vars(vars(),  ['rec_new["zz"][0].keys()'])
                    db.insert('correlation', sh, rec_new)
            try:
                update_correlation(self) 
                succeed = True
            except Exception as err:
                succeed = False
                msg = 'failed to update, err is {}'.format(str(err))
                warnings.warn(msg)
                raise  
            if succeed:
                self['version'] = ver 
                self.commit(info=1)

    def calc_correlation(self, sh, direct, i0_list=None, r_max=None, r_list=None, force=False, 
            use_local_storage=False, parpath_center=None, test_run=False,  
            info=0):
        """
            this is a temporary workaround 
            params:
                parpath_center:  needed only in remote compute
                    
        """
        from vmps.measure_and_analysis.measurement_idmrg_mcc import correlation  as func 
        if r_list is None: 
            if r_max is None:
                D = sh[1]
                r_max = 400 if D<= 160 else (1000 if D<=320 else 4000) 
            r_max_log10 = np.log10(r_max)    
            r_list = np.arange(0, r_max_log10, 0.02)
            r_list = (10**r_list).astype(int)
            r_list = list(set(r_list))
            r_list = [r for r in r_list if 0<r<10000]
            r_list.sort()
        
        #if use_local_storage:  #this seems strange, but, if without it, when run it distributively, only an empty db is transfered to the remote machine
        if self is not None:
            db = self
        else:
            db = ResultDB_idmrg(parpath_center, use_local_storage=1)
            print_vars(vars(),  ['db'])
        state = db.load_S(sh, use_local_storage=use_local_storage)
        period = state['period']
        i0_list = i0_list if i0_list is not None else list(range(period))
        calc_aver = len(i0_list)==period
        changed = False
        for i0 in i0_list:
            kk = ['correlation', sh, direct, i0]
            if not db.has_key_list(kk):
                db.add_key_list(kk) 
            res_old = db.fetch_from_key_list(kk)
            r_old = list(res_old.keys())
            
            r_list_calc = set(r_list)-set(r_old) if not force  else r_list
            print_vars(vars(),  ['len(r_old)', 'len(r_list)', 'len(r_list_calc)'])
            
            if len(r_list_calc)>0: 
                changed = True
                temp=func(state, [direct], [i0], r_list=r_list_calc, r_max=None, info=info)
            
                res_old.update(temp[direct][i0])
                res_new = OrderedDict(sorted(res_old.items()))
                
                db['correlation'][sh][direct][i0] = res_new 
        if changed: 
            if calc_aver:
                corr = db['correlation'][sh]
                if 'aver' not in corr:
                    corr['aver'] = OrderedDict()
                aver = corr['aver']
                temp = [(r, np.mean([corr[direct][i0][r] for i0 in i0_list]))
                            for  r in r_list]
                aver.update(temp)
                aver = OrderedDict(sorted(aver.items()))
                corr[direct]['aver'] = aver
            if not test_run:   
                db.commit(info=1, use_local_storage=use_local_storage)
        else:
            print('nothing changed')

    def calc_rdm_one_fermion(self, sh, i0_list=None, r_max=None, r_list=None, force=False, 
            use_local_storage=False, parpath_center=None, test_run=False,  
            info=0):
        """
            this is a temporary workaround 
            params:
                parpath_center:  needed only in remote compute
                    
        """
        from vmps.measure_and_analysis.measurement_idmrg_mcc import rdm_one_fermion as func 
        if r_list is None: 
            if r_max is None:
                D = sh[1]
                r_max = 400 if D<= 160 else (1000 if D<=320 else 4000) 
            r_max_log10 = np.log10(r_max)    
            r_list = np.arange(0, r_max_log10, 0.02)
            r_list = (10**r_list).astype(int)
            r_list = list(set(r_list))
            r_list = [r for r in r_list if 0<r<10000]
            r_list.sort()
        
        #if use_local_storage:  #this seems strange, but, if without it, when run it distributively, only an empty db is transfered to the remote machine
        if self is not None:
            db = self
        else:
            db = ResultDB_idmrg(parpath_center, use_local_storage=1)
            #print_vars(vars(),  ['db'])
        state = db.load_S(sh, use_local_storage=use_local_storage)
        period = state['period']
        i0_list = i0_list if i0_list is not None else list(range(period))
        calc_aver = len(i0_list)==period
        changed = False
        for i0 in i0_list:
            #kk = ['correlation', sh, direct, i0]
            kk = ['rdm_one_fermion', sh, i0]
            if not db.has_key_list(kk):
                db.add_key_list(kk) 
            res_old = db.fetch_from_key_list(kk)
            r_old = list(res_old.keys())
            
            r_list_calc = set(r_list)-set(r_old) if not force  else r_list
            print_vars(vars(),  ['len(r_old)', 'len(r_list)', 'len(r_list_calc)'])
            
            if len(r_list_calc)>0: 
                changed = True
                temp=func(state, [i0], r_list=r_list_calc, r_max=None)
            
                res_old.update(temp[i0])
                res_new = OrderedDict(sorted(res_old.items()))
                
                db['rdm_one_fermion'][sh][i0] = res_new 
        if changed: 
            if calc_aver:
                corr = db['rdm_one_fermion'][sh]
                if 'aver' not in corr:
                    corr['aver'] = OrderedDict()
                aver = corr['aver']
                temp = [(r, np.mean([corr[i0][r] for i0 in i0_list]))
                            for  r in r_list]
                aver.update(temp)
                aver = OrderedDict(sorted(aver.items()))
                corr['aver'] = aver
            if not test_run:   
                db.commit(info=1, use_local_storage=use_local_storage)
        else:
            print('nothing changed')

    #def calc_correlation_submit(self, sh, direct, i0_list=None, r_list=None):
    def calc_correlation_submit(self, sh, which, server=None, param=None):
    
        """
            this is a convenient function
        """
        
        from mypy.brokest.task_center import submit_one, LOCAL_IP
        from mypy.brokest.brokest import run_many 
            
        #db=db.__class__(parpath, use_local_storage=1)
        param = param if param is not None else {}
        
        kwargs = dict(use_local_storage=False, 
                parpath_center = self.parpath_center, 
                )
        kwargs.update(param)
        if which == 'correlation':
            func = self.__class__.calc_correlation.__func__
            direct = kwargs['direct']
            kwargs.pop('direct')
            args= (None, sh, direct)
        elif which == 'rdm_one_fermion' :
            func = self.__class__.calc_rdm_one_fermion.__func__
            args = (None, sh)
        else:
            raise  
            
        if 1:
            server = server if server is not None else (LOCAL_IP, 90999)
            #server = server if server is not None else ('localhost', 90999)
            status=submit_one(func, 
                    server = server, 
                    args=args, kwargs=kwargs, 
                    job_info = {'priority': 1, 'job_group_name': 'calcl-correlation', 'delay_send': 20, })
            print('\tresponse: ', status)
        else:
            tasks= [(func, args, kwargs)]
            run_many(tasks, servers=('localhost', None))
            #run_many(tasks, servers=('qtg7501', None))
            #run_many(tasks, servers=('qtg7501', 90901))


    def calc_correlation_bond(self, sh,  r_list=None, info=0):
        """
            this is a temporary workaround 
        """
        from vmps.measure_and_analysis.measurement_idmrg_mcc import correlation_bond  as func 
        from merapy.common_util import set_num_of_threads 
        set_num_of_threads(4) 
        r_list = [] if r_list is None else r_list 
        res_old = self['correlation_bond'][sh][-1]
        r_old = list(res_old.keys())
        r_list = set(r_list)-set(r_old)
        
        if len(r_list)>0: 
            state = self.load_S(sh)
            temp=func(state, r_list=[r[1] for r in r_list])
        
            res_old.update(temp)
            res_new = OrderedDict(sorted(res_old.items()))
            self['correlation_bond'][sh][-1] = res_new 
            self.commit(info=info)
    
    def calc_corr_conn(self, sh, force=False):
        if 1:
            if sh[1] == 'max' : 
                sh = sh[0], self['dim_max'].get(sh[0], 0)
            elif sh[1] == 'field_max':  #the max dim for each field may differ, then may use this
                N = sh[0]
                ss = list(self.get('correlation', {}).keys())
                dd = [s[1] for s in ss if s[0]==N ]
                D = max(dd) if dd else 0     
                sh = N, D
        name = 'corr_conn'
        rec = self.fetch_easy(name, sh, ['aver'], default={})
        if 1:
            corr = self.fetch_easy('correlation', sh, ['zz', 'aver'], default={}) 
            #print_vars(vars(),  ['len(rec), len(corr)'])
            if len(rec) < len(corr):
                print('updateing corr_conn')
                force = 1
                
        #return np.split(data, np.where(np.diff(data) != 1)[0]+1)
        if rec is None or force:
            corr = self.fetch_easy('correlation', sh, ['zz', 0], default={})
            rr = list(corr.keys())
            func = self._get_corr_conn(sh, aver=True)
            rec = OrderedDict()
            if func is not None and rr:
                if 0 not in rr:
                    rr.insert(0, 0)
                for r in rr:
                    try:
                        rec[r] = func(r)
                    except KeyError:
                        pass
            if len(rec)>1:
                self.insert(name, sh, rec, ['aver'])
                self.commit(info=1)
            else:
                warnings.warn('failed to calc corr_conn')
                rec = None
        
        return rec
    
    
    def _get_corr_conn(self, sh, aver=False, connected=True):
        """
            params:
                connected: if false,  turn off extract magnetization
        """
        corr = self.fetch_easy('correlation', sh, ['zz'])
        mag = self.fetch_easy('magnetization', sh, ['z'])
        if corr is None:
            warnings.warn('corr is needed')   
            return None
        if mag is None:
            warnings.warn('mag is needed')   
            return None
       
        p = len(list(mag.keys()))
        if not connected:
            mag = np.zeros(len(mag))
        def func(i0, r):
            if r == 0: #i=j  
                res = 0.25  #this value affects C(k=0) !!
            elif r>0:   #i<j 
                res = 0.25*(
                    corr[i0][r] - mag[i0]*mag[(i0+r)%p] )
            elif r<0:  #j<i
                res = 0.25*(
                    corr[(i0+r)%p][-r] - mag[(i0+r)%p]*mag[i0] )
            return res
        
        i0_list = list(range(p))
        def func_aver(r):
            res = 1./p*np.sum([func(i0, r) for i0 in i0_list])                
            return res
        if not aver:
            return func
        else:
            return func_aver
    
    def _get_corr_k(self, sh, k=None, x_min=-100, x_max=100, 
            connected=True, per_site=False,  use_stored=True, i0='aver'):
        """
            ref. Karrasch 2012  eq.17
            
            get Ck i.e. the static structure factor
            applicable to XXZ or spinless fermion like models
            
            mapping
            ni = (zi - 1)/2    #here minus sign, |0> corresponding to spin up 
            < ni nj > - <ni><nj>  = 0.25 * (<zi zj> -<zi><zj>)
            where zi is the pauli mat sigma^z_i
            
            note1:
                if the summation is from 0 to L, the result will be half valued
        """
        db = self
        if 1:
            aver = True if i0 == 'aver' else False
            if i0 == 'aver': 
                if not use_stored:  # this is somewhat slow, so I store it in db
                    corr_conn_func = db._get_corr_conn(sh, aver=True, connected=connected)
                    if corr_conn_func is None:
                        warnings.warn('failed to get corr_conn')
                        return None
                else:
                    corr_conn = db.calc_corr_conn(sh)
                    #note that for aver of corr: corr[r] = corr[-r] exactly!
                    corr_conn_func = lambda r: corr_conn[abs(r)]
                    if corr_conn is None:
                        warnings.warn('corr_conn is failed to calc')
                        return None
            else:
                temp = db._get_corr_conn(sh, aver=False, connected=connected)
                corr_conn_func = lambda r: temp(i0, r)
        #note1
        #Ck = lambda k, x_min=0, x_max=100: np.sum( 
        #    [ np.exp(-1j*k*x)* corr_conn_func(x)  for x in range(x_min+1, x_max)]   )
        N = x_max - x_min
        if per_site:
            Ck = lambda k: 1./N*np.sum( [ np.exp(-1j*k*x)* corr_conn_func(x)  for x in range(x_min+1, x_max)]   )
        else:
            Ck = lambda k: np.sum( [ np.exp(-1j*k*x)* corr_conn_func(x)  for x in range(x_min+1, x_max)]   )
            
        if k is None:
            return  Ck
        else:
            try:
                return Ck(k) 
            except KeyError as err:
                warnings.warn(str(err))
                return None

    def _get_Luttinger_param_fit_corr(self, sh, which, nu, 
            i0='aver', f=None, x_min=20, x_max=40, period=1, 
            strict=True, fault_tol=True):
        """
            params:
                strict: only for check data
        """
        db = self
        #if which == 'zz' :
        #    rec = db.calc_corr_conn(sh, force=0)
        #else:
        rec=db.fetch_easy('correlation', mera_shape=sh, 
                sub_key_list=[which, i0])
        if rec is not None:
            rec=list(rec.items())
            x,y=list(zip(*rec))
            data=list(zip(x,y))
        else:
             return None
        if strict:
            temp = [i for i in x if x_min<=i<=x_max ]
            if len(temp)-1 != x_max-x_min: 
                msg= 'data not complete: %s %d'%(self.parpath, len(temp))
                if fault_tol:
                    warnings.warn(msg)
                else:
                    raise  
            
        func= f if f is not None else self.corr_luttinger(which, nu)
        res=db.fit_line(None, func=func, data=data, period=period,
                        x_min=x_min, x_max=x_max, )
        if res is not None:
            K=res['param'][0]
        else:
            K = np.nan
       
        return K    

    def _get_K_by_corr_k(self, sh, L=100, i0='aver', 
            k_min=0.02*pi, k_max=0.08*pi):
        """
            get luttinger parameter form fitting of C(k) near k~0
            ref. Karrasch 2012  eq.19
        
        """
        db = self
        Ck=db._get_corr_k(sh, x_min=-L, x_max=L, i0=i0)
        if Ck is None: 
            return None
        kk=np.linspace(k_min, k_max, 100)
        #kk=np.linspace(0.0, 0.08*pi, 201)
        try:
            y=[Ck(k) for k in kk]
        except KeyError as err:
            warnings.warn(str(err))
            return np.nan
            
        data = list(zip(kk, y))
        func = lambda x,K, b: K/(2*pi)*x+b
        #func = lambda x,K, b: K/2*x+b
        try:
            res=db.fit_line(None, func=func, data=data, 
                        x_min=k_min, x_max=k_max, )
            K=res['param'][0]
        except RuntimeError as err:
            warnings.warn(str(err))
            K=np.nan
        except TypeError as err:
            warnings.warn(str(err))
            K=np.nan
        return K
    
    def _get_K_by_nk(self, sh, nu, x_min=20, x_max=48, period=None):
        """
            ref. Karrasch 2012  eq.8 and 9
            attention: 
                acctually not by nk, but by G(x) whose Fourier trans is nk. 
                I use the name nk just because it is easier to remember
        
        """
        db = self
        rec = self.fetch_easy('rdm_one_fermion', sh, ['aver']) 
        if rec is None:
            return None
        rec=list(rec.items())
        x,y=list(zip(*rec))
        x = np.asarray(x)
        y = np.asarray(y)
        p = period if period is not None else {0.5:2, 0.33:3}[nu]
      
        sin, pi = np.sin, np.pi
        if 1:
            
            func = lambda x, K_, C: C*x**(-K_)        
            arg=self.filter_array(x, x_min=1, x_max=None, period=p)
            x = x[arg]
            y = y[arg]
            y = np.abs(y)  # note abs at here !
            
            res=db.fit_line(None, func=func, data=list(zip(x, y)), period=1,
                            x_min=x_min, x_max=x_max, )
            K = res['param'][0] if res is not None else np.nan
            yfunc=lambda x:x - np.sqrt(x**2-1)   # this operation sqares error around K=1.0  !
            K = yfunc(K)
        else:
            #func = lambda x, K, C:   -C*x**(-K/2.-1./(2.*K))
            #func = lambda x, K, C:  -C*sin(pi*x/p)/pi*x**(-K/2.-1./(2.*K))
            func = lambda x, K, C:  -C*(-1)**(x+1)*sin(pi*x/p)/pi*x**(-K/2.-1./(2.*K))
            #func = lambda x, K:  -sin(pi*x/2.)/pi*x**(-K/2.-1./(2.*K))
            #p0 is must, as there seems exists local minima!
            res=db.fit_line(None, func=func, data=list(zip(x, y)), period=1,
                    method = 'lm', 
                    p0 = (0.8, 1.0), 
                    x_min=x_min, x_max=x_max, )
            K = res['param'][0] if res is not None else np.nan
        #print_vars(vars(),  ['res["param"]', 'res["cov"][0, 0]'])
        
       
        return K    

    def reload_class(self): 
        module_name =self.__class__.__module__
        #print module_name
        module=importlib.import_module(module_name)
        reload(module)
        name = self.__class__.__name__
        self.__class__=getattr(module, name)
        print('class %s is reloaded'%(name, ))

    def _get_N(db, sh):
        try:
            s=db.load_S(sh)
        except Exception as err:
            warnings.warn(str(err))
            s={}
        N=s.get('N', None)
        return N


class ResultDB_tdvp(ResultDB): 
    def __init__(self, parpath,  **kwargs): 
        kwargs['version'] = 1.0
        kwargs.update(algorithm='tdvp')
        ResultDB.__init__(self, parpath,  **kwargs)
    
    def get_shape_list(self,  N=None, D=None, field='the_time', **kwargs):
        return super(ResultDB_tdvp, self).get_shape_list(N=N, D=D, field=field, **kwargs)
    
    def get_t_vs_field(self, field_name, sh, tlim=None, ts=None, integrate=False,  kac_rescale=None,  force=False):
        """
             
        """
        if ts is None:
            ts=self.get_time_serials(sh, force=force)
        if ts is None:
            return None, None
        temp = ts.values()
        tt = [i['the_time'] for i in temp]
        #val = [i[field_name] for i in temp]
        val = [i.get(field_name, None) for i in temp]
        tt = np.asarray(tt)
        val = np.asarray(val)
        if integrate:
            val = scipy.integrate.cumtrapz(val, tt, initial=0)
        
        if tlim is not None:
            arg = np.where(tt<tlim)
            tt = tt[arg]
            val = val[arg]
            
        if kac_rescale is not None:
            tt = tt/kac_rescale
            
        
        return tt, val 
    
    def get_current_from_mag(self, sh, mu=None, integrate=False, offset=True, force=False):
        """
            An indirect measure of current. 
            
            Note that the current is AVERAGED over interval of dt. Its value
            may be not completely concide with direct measure of current,  as
            the later is INSTANEOUS.  The descrypency is larger when the
            current changes with time faster. 
            
        """
        #tt, current_aver = self.fetch_easy('current_aver', sh)
        
        tt, mm = self.get_t_vs_field('magnetization',  sh,  force=force)
        if tt is None:
            return None
        
        N_half = sh[0]//2
        
        temp = [mm[i]['z'].values() for i in range(len(tt))]
        temp = [list(t) for t in temp]
        temp = np.asarray(temp, dtype=float)
        
        temp_l = temp[:, :N_half]
        mag_l = np.sum(temp_l, axis=1)
        mag_l *= 0.5  # magnetization was measure with sigma_z. Here multiply 0.5 to match with the def of the current 
        mag_diff = np.diff(mag_l)
        t_diff = np.diff(tt)
        current_aver = mag_diff/t_diff
        if offset: # use midpoint of t + 0.5dt as the time for the current 
            tt = tt[:-1]
            tt += 0.5*t_diff
            current_aver = np.insert(current_aver, 0, 0)  # I assume at t=0 the current_aver=0
            tt = np.insert(tt, 0, 0)
        else:
            current_aver = np.insert(current_aver, 0, 0)
        
        if integrate:
            current_aver = scipy.integrate.cumtrapz(current_aver, tt, initial=0)
        
        return tt, current_aver
    
    def convert_t_to_temperature(self, time_list):
        res = [1/(-i.imag*2) for i in time_list]
        return res 
    T  =  convert_t_to_temperature

    def get_imbalance(self, sh):
        db = self
        #s=db.load_S(sh)
        #if s is None:
        #    return None, None
        #ts=s['time_serials']
        ts=db.get_time_serials(sh)
        if ts is None:
            return None, None
        temp = ts.values()
        tt = [i['the_time'] for i in temp]
        mag=[list(i['magnetization']['z'].values()) for i in temp]
        mag=np.asarray(mag).real
        L=mag.shape[1]
        ii = [(-1)**(i+1) for i in range(L)]
        ii = np.asarray(ii)
        mag=mag*ii
        imb = np.sum(mag, axis=1)/L 
        return tt, imb

    def get_ee(self, sh):
        db = self
        ts=db.get_time_serials(sh)
        if ts is None:
            return None
        temp = ts.values()
        tt = [i['the_time'] for i in temp]
        ee=[i['EE_middle'] for i in temp]
        #mag=np.asarray(mag).real
        if 0:
            L=mag.shape[1]
            ii = [(-1)**(i+1) for i in range(L)]
            ii = np.asarray(ii)
            mag=mag*ii
            imb = np.sum(mag, axis=1)/L 
        return tt, ee

    def calc_diffuse_const(db, sh, tlim=None, force=False):
        """
            calc diffuse const through measure integrate of j(t)
            ref:
               Ljubotina, Prosen 2017  eq.4 and Fig. 2
        
        """
        tt, yy = db.get_t_vs_field('current', sh, tlim=tlim, force=force)
        if tt is None:
            return None, None
        yy = scipy.integrate.cumtrapz(yy, tt, initial=0)

        tt_log=np.log(tt)
        yy=np.log(np.abs(yy))

        tt_diff= np.diff(tt_log)
        yy_diff=np.diff(yy)
        yt_diff= yy_diff/tt_diff
        return tt[:-1], yt_diff
    
    def calc_diffuse_coeff(self, sh, tlist, force=False):
        """
            params:
                tlist: can either float or list of float
            ref:
                Ljubotina 2017 inset of Fig.2 
            
        """
        pass
        tt, jj = self.get_t_vs_field('current', sh, force=force)
        if tt is None:
            return None
        tt1, mm = self.get_t_vs_field('magnetization',  sh)
        if tt1 is None:
            return None
        N = sh[0]
        i = N//2-1 
        DD = []
        if isinstance(tlist,  float) or isinstance(tlist, int):
            tlist = [tlist]
        for t in tlist:
            ind = np.where(tt==t)[0]
            if not ind:
                D = None 
            else:
                mag=mm[ind][0]['z']
                    
                j = jj[ind][0]
                j = j.real 
                m_diff = mag[i+1] - mag[i]
                D = j/m_diff
            DD.append(D)
            
        DD = np.asarray(DD)
        if len(DD)==1:
            DD = DD[0]
            
        return DD

    def get_magnetization(self, sh, tlist, ts=None, sub_key_list='z',  return_t=False):
        """
            tince magnetization is frequently used. I write this function. 
        """
        tt, mm = self.get_t_vs_field('magnetization',  sh, ts=ts)
        if tt is None:
            return None, None
        tmax = abs(max(tt))
        #if isinstance(tlist,  float) :
        if isinstance(tlist,  float) or isinstance(tlist, int):
            tlist = [tlist]
        tlist = [t for t in tlist if t < tmax]
        ii = []
        for t in tlist:
            ind = np.where(tt==t)[0]
            assert len(ind) == 1, (t, tt)
            ii.append(ind[0])
        #tt = [tt[i][0].real for i in ii]
        tt = [tt[i].real for i in ii]
        tt = np.asarray(tt)
        
        if isinstance(mm, np.ndarray) and mm.ndim==2:
            temp = mm[ii, :]
        else:
            if sub_key_list == 'z' :
                #temp = [mm[i][0]['z'].values() for i in ii]
                temp = [mm[i]['z'].values() for i in ii]
            else:
                temp = [mm[i].values() for i in ii]   # e.g. for measure boson numbers
            
            temp = [list(t) for t in temp]
            temp = np.asarray(temp, dtype=float)
            
        if len(tlist)==1:
            temp = temp[0]
        if len(temp)==0:
            temp = None
        if return_t:
            return tt, temp
        else:
            return temp

    def get_ball_center(tt, data, tmin=None, tmax=None, xmin=None, xmax=None):
            arg = None
            if tmin is not None:
                arg = tt>= tmin
            if tmax is not None:
                arg = arg * (tt<=tmax)
            
            if arg is not None:
                arg = np.where(arg)[0]
                tt = tt[arg]
                data = data[arg, :]
            
           
            data = data.copy()
            data[:, :]  += 1  
            
            N = data.shape[1]
            
            xmin = xmin if xmin is not None else 0
            xmax = xmax if xmax is not None else N
            x0 = np.ndarray((len(tt), xmax-xmin), dtype=int)
            x0[:] = range(xmin, xmax)
            aver0 = np.average(x0, axis=1, weights=data[:, xmin:xmax])
            return tt, aver0


class ResultDB_ed(ResultDB_tdvp): 
    def __init__(self, parpath,  **kwargs): 
        kwargs['version'] = 1.0
        kwargs.update(algorithm='ED')
        ResultDB.__init__(self, parpath,  **kwargs)
    
    def get_energy_lowest(self, sh, i, qn_exclude=None, qn_only=None): 
        qn_exclude = qn_exclude if qn_exclude is not None else [] 
        rec = self.fetch_easy('energy', sh)
        temp = []
        if rec is not None:  
            for k, v in list(rec.items()): 
                if k not in qn_exclude: 
                    temp.extend(v)
        else: 
            temp = [None]*10
        temp.sort()
        return temp[i]



class ResultDB_proj_qmc(ResultDB): 
    def __init__(self, parpath,  **kwargs): 
        kwargs['version'] = 1.0
        kwargs.update(algorithm='proj_qmc')
        ResultDB.__init__(self, parpath,  **kwargs)
    
    @staticmethod
    def shape_to_backup_fn(sh): 
        N, m = sh
        res = 'N=%(N)d-m=%(m)d.pickle'%vars()
        return res
    backup_fn_gen = shape_to_backup_fn 

    def get_Svb(self, sh):
        N = sh[0]
        data = self.fetch_easy('svb_sample', sh)
        #data=(db['svb_sample'][(N, N*20)])
        if data is not None: 
            nn,yy=list(zip(*data))
            n=np.sum(nn)
            x=np.arange(N//2)
            y=reduce(lambda x,y:x+y,yy, np.zeros(N//2, dtype=int))*np.log(2.)/n 
            res = list(zip(x, y))
        else:
            res= None
        return res 

class ResultDB_bethe_ansatz(ResultDB):
    """
        for bethe ansatz
    """
    def calc_energy(self, L, tol, k0=0):
        from bethe_ansatz_xxz.luttinger_parameter_xxz import calc_energy as calc_energy_func
        sh = L, tol
        energy = self.fetch_easy('energy', sh)
        if energy is None:
            parpath = self.parpath
            fn = os.path.basename(parpath)
            Jzz = fn.split('-')[0].split('=')[1]
            Jzz = eval(Jzz)
            n = fn.split('-')[1].split('=')[1]
            n = eval(n)
            n = {0.5:0.5, 0.33:1./3}[n]
            N = int(L*n)
            assert L%N == 0, (L, N)
            
            try:
                energy=calc_energy_func(L, N, Jzz, k0=k0, tol=tol)
                self.insert('energy', sh, energy)
                self.commit(info=1)
            except Exception as err:
                print('calc energy failed')
                print(err)
        return energy 

    def calc_K(self, L, tol, k0=0, force=False):
        #from bethe_ansatz.xxz.luttinger_parameter_xxz import calc_energy as energy_func
        from bethe_ansatz_xxz.luttinger_parameter_xxz import calc_K  as calc_K_func
        
        sh = L, tol
        K = self.fetch_easy('K', sh)
        if K is None or force:
            parpath = self.parpath
            fn = os.path.basename(parpath)
            Jzz = fn.split('-')[0].split('=')[1]
            Jzz = eval(Jzz)
            n = fn.split('-')[1].split('=')[1]
            n = eval(n)
            n = {0.5:0.5, 0.33:1./3}[n]
            N = int(L*n)
            assert L%N == 0, (L, N)
            try:
                K=calc_K_func(L, N, Jzz, k0=k0, tol=tol)
                self.insert('K', sh, K)
                self.commit(info=1)
            except Exception as err:
                print('calc K failed')
                K = None
                print(err)
        return K 
        
class TestResultDB(unittest.TestCase): 
    def setUp(self):
        self.seq = list(range(10))

        args= {}
        #parpath = '/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=2.0'
        parpath = '/tmp/'
        sh = (12, 5);  
        args['dir'] = parpath
        args['mera_shape_list'] = sh
        args['draw'] = 1
        which='fit_correlation'
        args['which'] = which
        parpath = args['dir'] 
        parpath_idmrg = '/home/zhli/mps_backup_folder/run-heisbg-long/alpha=2.0'
        if 0: 
            self.db = ResultDB(parpath)
        
    def test_temp(self): 
        N = 128
        
        #from vmps.run_heisenberg.analysis import an_tdvp 
        #from mps_wigner_crystal.analysis import an_tdvp 
        from mps_wigner_crystal.analysis import an_exact_diag
        

        xx = an_exact_diag.an_dynamics.an_random_spin
        xx.outline()
        db = xx(alpha=1.0, V=100.0, spin_up=0.3, dt=100)
        print_vars(vars(),  ['db'])
        
        #def get_t_vs_field(self, field_name, sh, tlim=None, ts=None, integrate=False,  kac_rescale=None,  force=False):
        tt,  mag = db.get_t_vs_field('magnetization', (12, 'max'))
        mag = db.get_magnetization((12, 'max'), tlist=[100, 200, 300] )
        print_vars(vars(),  ['mag'])
        raise  
        vv=[0.0, 0.5, 1.0, 4.0, 16.0]
        
        i=0
        v = 0.0
        
        db=xx(Jzz=v, dt=0.5, surfix='fix_err_1em12')
        sh=(64, 'max')
        tt=np.arange(1., 180, 1.0)
        ii=range(sh[0])
        tt, data = db.get_magnetization(sh, tt, return_t=1,)
        print_vars(vars(),  ['db.dir_name', 'db.model_param'])
        raise  
        print_vars(vars(),  ['data.shape'])
        

            
       
        tt, a = get_ball_center(tt, data, 2, 8, 1, 32)
        print_vars(vars(),  ['len(a)'])
        print_vars(vars(),  ['a'])
        #print_vars(vars(),  ['np.sum(data[1] + 1)'])
        
        
        raise  
       
        
    def test_insert_and_fetch(self): 
        a = [1, 2, 3, 4]
        dic=ResultDB.make_nested_dict(a, val='love')
        a = [1, 2, 3, 5]
        dic=ResultDB.make_nested_dict(a, dic, val='hate')
        self.assertEqual(dic[1][2][3][4], 'love')
        self.assertEqual(dic[1][2][3][5], 'hate')
        
        db = ResultDB('/tmp')
        #db.insert( 'xxx', (0, 4), -1, 'love', [1])
        #db.insert( 'xxx', (0, 4), -1, 'hate', [2])
        db.insert( 'xxx', (0, 4), 'love', [1])
        db.insert( 'xxx', (0, 4), 'hate', [2])
        print(db['xxx']) 
        rec = db.fetch_easy('xxx', (0, 4), [1])
        self.assertTrue( rec=='love')
        rec = db.fetch_easy('xxx', (0, 4), [2])
        self.assertTrue( rec=='hate')
        print(' ========================')
        rec = db.fetch_easy('xxx', (0, 4), [1000])  # 1000 is a wrong key 
        self.assertTrue(rec is None)
        
    
    def test_idmrg_get_EE(self): 
        parpath =  '/home/zhli/mps_backup_folder/run-heisbg-long/idmrg/alpha=2.0'
        db = ResultDB_idmrg(parpath)
        
        print(db.get_EE(2))
   
    def test_add_key_list(self): 
        xx = [1, 2, 3, 4]
        val = 'xxx'
        self.db.add_key_list(xx, val='xxx')
        self.assertTrue(self.db.has_key_list(xx ))
        self.assertFalse(self.db.has_key_list(xx + [5]))
        

if __name__ == '__main__': 
    #warnings.filterwarnings('ignore')
    warnings.filterwarnings('once')
    if 1: 
        FIELD_NAME_LIST = [
            'energy', 
        'central_charge', 'scaling_dim', 'correlation', 'correlation_extra',  
        'entanglement_entropy', 'entanglement_spectrum', 'magnetization', 
        'entanglement_brute_force', 'entanglement_brute_force_6', 
        # these are too slow and not useful,  disable them
        #'entanglement_brute_force_9',   
        #'entanglement_brute_force_aver', 
        #'entanglement_brute_force_6_aver',
        #'entanglement_brute_force_9_aver' , 
        ]
           
        which_list = FIELD_NAME_LIST + ['structure_factor', 'fit_correlation']
        temp  =  list(ResultDB.__dict__.keys())
        which_list = [x for x in temp if 'plot' in x or 'show' in x or 'get' in x]
        which_list += ['fit_correlation'] 
        #import types 
        #temp = filter(lambda x: isinstance(temp[x], types.FunctionType) , temp.keys())
        
        parser = argparse.ArgumentParser(description='dont know')
        #parser.add_argument('-dir', '--dir', default='./')
        parser.add_argument('dir', nargs='?', default='./')
        parser.add_argument('-direct', '--direct', default='zz', choices=['pm', 'zz', 'xx', 'yy'])
        parser.add_argument('-w', '--which', choices=which_list,  default=None)
        parser.add_argument('-sw', '--switch', default=None)
        #parser.add_argument('-sh', '--mera_shape_list', nargs='*', default='all')
        parser.add_argument('-sh', '--mera_shape_list', nargs=2, type=int,  default=None)

        #parser.add_argument('-show', '--show', action='store_true', default=False) 
        parser.add_argument('-show', '--show', default=False) 
        argcomplete.autocomplete(parser)
        args = parser.parse_args()
        args = vars(args)
        if args['mera_shape_list'] is not None : 
            args['mera_shape_list'] = tuple(args['mera_shape_list'])

        
        sh = args['mera_shape_list']
        which = args['which']
        if which is not None: 
            
            args.pop('which')
            res=db.__getattribute__(which)(sh=sh, show_fig=1, show_val=1, **args)
            print(res)
            #db.upgrade()
            #db.backup()

    #if len(sys.argv)<=1 :   
    if 0: 
        args= {}
        parpath = '/home/zhli/Documents/mera_backup_tensor/run-long-better/alpha=0.3'
        #parpath = '/tmp/'
        sh = (4, 4);  
        args['dir'] = parpath
        args['mera_shape_list'] = sh
        args['draw'] = 1
        #args['show_fig'] = 1
        args['which'] = 'plot_structure_factor' 
        args['switch'] = 'diff'
        args['fit'] = 'power'
        args['plot'] = 1
        args['plot_fit'] = 1
        #args['xlim'] = (0, 10000)
        
        parpath = args['dir'] 
        db = ResultDB(parpath)
    
        if 0: 
            sh = args['mera_shape_list']
            #sh[0] = eval(sh[0]); sh[1]=eval(sh[1])
            which = args['which']
            
            #if which in which_list: 
            if 1: 
                args.pop('which')
                res=db.__getattribute__(which)(sh=sh, show_fig=1, show_val=1, **args)
                print(res)
            else: 
                pass
                #db.upgrade()
                #db.backup()

    
    if len(sys.argv)<=1 :   
      
        if 0:
            TestResultDB.test_temp=unittest.skip("skip test_temp")(TestResultDB.test_temp) 
            unittest.main()
        else: 
            suite = unittest.TestSuite()
            add_list = [
               #'test_insert_and_fetch', 
               #'test_add_key_list', 
               #'test_idmrg_get_EE', 
               'test_temp', 
  
            ]
            for a in add_list: 
                suite.addTest(TestResultDB(a))
            unittest.TextTestRunner().run(suite)
           
        
        
    
