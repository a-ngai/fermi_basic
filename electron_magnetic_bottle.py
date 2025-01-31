# %%
"""
This scripts loads magnetic bottle data. It does a few things:
    1. Separates foreground and background shots (determined by "BACKGROUND_PERIOD" in the .h5 files)
    2. Calibration for the electron magnetic bottle
    3. Rebinning of data
"""

# %%  Getting electron magnetic bottle data
import h5py
import numpy as np
import matplotlib.pyplot as plt

filename = "./test_data/Run_001/rawdata/Run_001_012345678.h5"
with h5py.File(filename, 'r') as f:
    bunches = f['bunches'][()]  # each FEL laser pulse (i.e. bunch) is numbered as a "bunch number"
    background_period: int = f['Background_Period'][()]
    eon_tof_all = f["digitizer/channel3"][()]  # groupname for electron magnetic bottle
    
# separate foregorund from background
background_shots = bunches% background_period == 0
foreground_shots = ~background_shots
eon_tof_fore = eon_tof_all[foreground_shots]
eon_tof_back = eon_tof_all[background_shots]

avg_eon_tof_fore = np.average(eon_tof_fore, axis=0)
avg_eon_tof_back = np.average(eon_tof_back, axis=0)
avg_eon_tof = avg_eon_tof_fore - avg_eon_tof_back
avg_eon_tof /= -1  # signals are recorded as voltage drops, so we flip them

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

timezero = 5000  # typically find the photopeak by averaging the TOF, and looking for a small feature

initial_params = lmfit.Parameters()
initial_params.add_many(
        # name, initial value, vary, min, max
        ('C', 1.2e-6, True, 0, None),
        ('T0', timezero, False, timezero, None),
        ('KE0', retardation_voltage, False, 0, None),
        )

# for time-zero, use an arbitrary large KE value
tof_peaks_ns = np.array([5800, 6100])
eke_peaks_ev = np.array([20, 10])  

fit_results = lmfit.minimize(residuals, initial_params,
                             args=(ke_fit_model, tof_peaks_ns, eke_peaks_ev))
fit_params = fit_results.params
propconst = fit_params['C'].value
timezero = fit_params['T0'].value
ke0 = fit_params['KE0'].value

# %%% Defining helper transformation functions here
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

# transform TOF to eKE spectrum w/ Jacobian correction
tof_ns = np.arange(len(avg_eon_tof))
eke_ev = ke_coordinate_func(tof_ns)
avg_eon_eke = avg_eon_tof * ke_jacobian_func(tof_ns)

not_nan = np.isfinite(eke_ev)
eke_ev = eke_ev[not_nan][::-1]
avg_eon_eke = -avg_eon_eke[not_nan][::-1]

# transform backwards (for completeness)
_tof_ns = tof_coordinate_func(eke_ev)
_avg_eon_tof = avg_eon_eke * tof_jacobian_func(eke_ev)


fig, ax = plt.subplots(1, 1, figsize=(3, 3))
_model_tof = np.linspace(0, np.max(tof_peaks_ns) + 1000, num=100)
_model_eke = ke_coordinate_func(_model_tof)
ax.plot(tof_peaks_ns, eke_peaks_ev, marker="o", linestyle="", color="black", label="calibration points")
ax.plot(_model_tof, _model_eke, label="calibration result")
ax.set_ylim(0, 50)
ax.set_xlabel("TOF (ns)")
ax.set_ylabel("eKE (eV)")
ax.set_title("show calibration points")
plt.show()

# %% Example plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3))
ax1.plot(tof_ns, avg_eon_tof)
ax1.set_xlabel("TOF (ns)")
ax1.set_ylabel("Signal (arb.u.)")
ax1.set_title("Example electron magnetic bottle spectrum\nTOF-coordinate")

show_eke_window = eke_ev < 50
ax2.plot(eke_ev[show_eke_window], avg_eon_eke[show_eke_window])
ax2.set_xlabel("eKE (eV)")
ax2.set_ylabel("Signal (arb.u.)")
ax2.set_title("\neKE-coordinate")

ax1.set_xlim(0, 10000)
ax2.set_ylim(0, 10000)
fig.tight_layout()

plt.show()

# %% Simple rebinning and linearization of eKE coordinate

# Note: a large negative signal near 0-0.1 eV is expected; this comes from a
# small negative baseline at large TOF. Re-binning integrates this small
# negative baseline over a large TOF range (low eKE encompasses many TOF
# points) gives a big negative signal!

from scipy.interpolate import interp1d
from scipy.integrate import cumulative_simpson

new_eke_ev = np.linspace(np.min(eke_ev), 50, num=1000)
not_nan = np.isfinite(avg_eon_eke)
eke_spec_interpolation = interp1d(eke_ev, cumulative_simpson(avg_eon_eke[not_nan], x=eke_ev[not_nan], initial=0))
new_avg_eon_eke = np.gradient(eke_spec_interpolation(new_eke_ev)) / np.gradient(new_eke_ev)

fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.plot(new_eke_ev, new_avg_eon_eke, color='black', alpha=0)
ylim = ax.set_ylim(0, None)
ax.plot(eke_ev[show_eke_window], avg_eon_eke[show_eke_window], label="original eKE spec.")
ax.plot(new_eke_ev, new_avg_eon_eke, label="rebinned spec.")
ax.set_xlabel("eKE (eV)")
ax.set_ylabel("Signal (arb.u.)")
ax.set_title("Example electron magnetic bottle spectrum")
ax.legend()
fig.tight_layout()