from merapy.measure_and_analysis.correlation import correlation, eta_fit, curve_fit, correlation_extra
from merapy.measure_and_analysis.central_charge import *#central_charge, entanglement_entropy, entanglement_spectrum
from merapy.measure_and_analysis.scaling_dimension import calc_scaling_dim_2site as scaling_dim
from merapy.measure_and_analysis.magnetization import magnetization
from merapy.measure_and_analysis.result_db import ResultDB
#from merapy.measure_and_analysis.measurement import Measurement, measure_all
from merapy.measure_and_analysis.tabulate import tabulate


def energy(S) : 
    if isinstance(S, dict): 
        return S['energy']
    else: 
        return S.energy 


