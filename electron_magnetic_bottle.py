# %%
"""
This scripts loads magnetic bottle data. It does a few things:
    1. Separates foreground and background shots (determined by "BACKGROUND_PERIOD" in the .h5 files)
    2. Calibration for the electron magnetic bottle
"""

# %%  Getting electron magnetic bottle data
import numpy as np
import h5py
import matplotlib.pyplot as plt

filename = "./test_data/Run_001/rawdata/Run_001_0000000.h5"
with h5py.File(filename, 'r') as f:
    bunches = f['bunches'][()]  # each FEL laser pulse (i.e. bunch) is numbered as a "bunch number"
    background_period: int = f['BACKGROUND_PERIOD'][()]
    eon_data = f["digitizer/channel3"][()]  # groupname for electron magnetic bottle
    
# separate foregorund from background
background_shots = bunches% background_period == 0
foreground_shots = ~background_shots
eon_foreground = eon_data[foreground_shots]
eon_background = eon_data[background_shots]

avg_eon_foreground = np.averager(eon_data, axis=0)
avg_eon_background = np.averager(eon_data, axis=0)
avg_eon_subtracted = avg_eon_foreground - avg_eon_background

# %%  Calibration of electron magnetic bottle
import lmfit

def residuals(params, model, x, y):
    return y - model(params, x)

def ke_fit_model(params, tof):
    """ functional form for the calibration; different people use different conventions! e.g. common to see C^2 instead o C """
    C = params['C']  # proportionality constant
    T0 = params['T0']  # time-zero
    KE0 = params['KE0']  # minimum kinetic energy (retardation voltage)
    ke = 1 / (C * (tof - T0)**2) + KE0
    return ke

with h5py.File(filename, 'r') as f:
    try:
        voltage_1 = f['endstation/MagneticBottle/voltage_ch1'][()]
        voltage_2 = f['endstation/MagneticBottle/voltage_ch2'][()]
        voltage_3 = f['endstation/MagneticBottle/voltage_ch3'][()]
        voltage_1_on = f['endstation/MagneticBottle/ch1_is_enabled'][()]
        voltage_2_on = f['endstation/MagneticBottle/ch2_is_enabled'][()]
        voltage_3_on = f['endstation/MagneticBottle/ch3_is_enabled'][()]
        retardation_voltage = voltage_1*voltage_1_on + voltage_2*voltage_2_on - voltage_3*voltage_3_on
    except Exception as e:
        print("Cannot determine retardation voltage from files. Following exception occured:")
        print(e)
        retardation_voltage = 0

initial_params = lmfit.Parameters()
initial_params.add_many(
        # name, initial value, vary, min, max
        ('C', 1, True, 0, None),
        ('T0', 0, True, 0, None),
        ('KE0', retardation_voltage, False, 0, None),
        )

# for time-zero, use an arbitrary large KE value
tof_ns = np.array([5020, 6300])
ke_ev = np.array([1000, 32])  

fit_results = lmfit.minimize(residuals, initial_params,
                             args=(ke_fit_model, tof_ns, ke_ev))
fit_params = fit_results.params
propconst = fit_params['C'].value
timezero = fit_params['T0'].value
ke0 = fit_params['KE0'].value

# %%% defining transformation functions here
def ke_coordinate_func(tof):
    ke_coordinate = tof*np.nan  #  helps avoid unphysical values, and errors from square-rooting negative numbers
    ke_coordinate[tof>timezero] = 1/propconst / (tof[tof>timezero] - timezero+0.)**2 + ke0
    return ke_coordinate

def ke_jacobian_func(tof):
    """ jacobian = |dT/dKE| """
    jacobian = tof*np.nan
    jacobian[tof>timezero] = (tof[tof>timezero] - timezero+0.)**3 / (-2 * 1/propconst)
    return jacobian

def tof_coordinate_func(ke):
    tof = ke*np.nan
    tof[ke<0] = np.sqrt(1/propconst / (ke[ke<0] - ke0+0.)) + timezero
    return tof

def tof_jacobian_func(ke):
    """ jacobian = |dKE/dT| """
    jacobian = (-2 * 1/propconst) / (tof_coordinate_func(ke) - timezero+0.)**3
    return jacobian
