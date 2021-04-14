#%%

import numpy as np 
from scipy.io.wavfile import read
import matplotlib.pyplot as plt
import warnings


def diff(x: np.array, 
         tau_max: int, 
         start: int,
         w: int) -> np.array:
    """
    B. Step 2: Difference Function. Calculate the difference function of one frame
    in Eq (6). 
    This is a sequential implementation, to be optimised

    Args:
        x (np.array): The input signal
        tau_max (int): An upper limit for the lag parameter tau (exclusive)
        start (int): The starting index in the signal
        w (int): integration window size

    Returns:
        np.array: A list of results calculated with different values of tau
    """
    
    # w < len(x)
    # tau_max - 1 + start + max_j (w) < len(x) => tau_max < len(x) + 1 - start - w
    assert tau_max < len(x) + 1 - start - w
    difference: list = []
    x = np.array(x, dtype=np.int32)
    for t in range(tau_max): 
        df: float = 0
        df2 = 0
        for j in range(1, w + 1): 
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    df += (x[start + j] - x[start + j + t]) ** 2
                except Exception:
                    print(f'{x[start + j]} - {x[start + j + t]} gave a warning')
        difference.append(df)
    
    return np.array(difference, dtype=np.float64)

def cmndiff(difference: np.array) -> np.array:
    """C. Step 3: Cumulative Mean Normalized Difference Function. Calculate the cmndiff of 
    one frame in Eq (8)

    Args:
        difference (np.array): the difference function calculated from diff

    Returns:
        np.array: a list of results
    """
    
    # add a 0 at the beginning to make the index align 
    diff_cumulative = difference.cumsum()
    diff_zero_tau = 1
    diff_positive_tau = difference[1:] / ((1 / np.array([i for i in range(1, len(difference))])) * difference.cumsum()[1:])
    
    return np.concatenate((np.array([diff_zero_tau]), diff_positive_tau))

def abs_threshold(cmndiff: np.array, threshold: float = 0.1) -> int:
    """D. Step 4: Absolute Threshold. Choose an absolute threshold and take the minimum value of tau
    that gives the minimum of d' deeper than that threshold. If none is found, choose the global minimum instead.
    A minimum of d' deeper than the threshold is equivalent to choosing one threshold and find the minimum tau
    that gives the value that falls below this new threshold. So we will use threashold as (the threshold stated in paper - d')

    Args:
        cmndiff (np.array): The cumulative mean normalized difference function calculated from cmndiff.
        threshold (float, optional): The threshold set to prevent "octave error". Defaults to 0.1.

    Returns:
        int: returns the value of tau chosen
    """
    tau_global_mean = len(cmndiff) - 1
    for t in range(len(cmndiff)):
        if cmndiff[t] <= threshold:
            return t # return the minimum value of tau that gives the value below the threshold
        
        if cmndiff[t] < cmndiff[tau_global_mean]:
            tau_global_mean = t
    
    return tau_global_mean

def parabolic_interpolation(tau_selected: int, cmndiff: np.array):
    """E. Step 5: Parabolic Interpolation. Perform parabolic interpolation on the difference function calculated.

    Args:
        tau_selected (int): the tau value selected for minimum difference
        cmndiff (np.array): the result calculated from cmndiff function.
    """
    assert tau_selected > 0 and tau_selected < len(cmndiff) - 1
    ordinates = np.array(cmndiff[tau_selected - 1: tau_selected + 2]) # get the y coordinates
    abscissae = np.array([tau_selected - 1, tau_selected, tau_selected + 1])
    
    coeffs = np.polyfit(abscissae, ordinates, 2)
    p = np.poly1d(coeffs)
    
    critical_pts = p.deriv().r
    real_critical_pts = critical_pts[critical_pts.imag==0].real
    critical_pt = real_critical_pts[0] # take the critical point check if it's between the first and third points
    # if it is, then use it as the result, if not (should be impossible, implicitly)
    
    if critical_pt > tau_selected - 1 and real_critical_pts < tau_selected + 1:
        return critical_pt
    else:
        return min(abscissae, key=lambda abscissa: cmndiff[abscissa]) # return the minimum of the three points
     
    
#%%
# test on sine waves
fs, data = read('audio/flute-alto-C-corrected.wav')

f1 = 400
f2 = 800
f3 = 1200
x = np.arange(10000)
# data = np.sin(2 * np.pi * f1 * x / fs) + np.sin(2 * np.pi * f2 * x / fs) + np.sin(2 * np.pi * f3 * x / fs)
diff_signal = diff(data, 3000, 10000, 200)
cmndiff_signal = cmndiff(diff_signal)
tau = abs_threshold(cmndiff_signal, threshold=0.1)
tau_interpolated = parabolic_interpolation(tau, cmndiff_signal)
print(f'tau (not interpolated): {tau}')
print(f'tau chosen: {tau_interpolated}')
detected_freq = 1 / (tau_chosen / fs)
print(f'detected frequency: {detected_freq}')

plt.figure(1)
plt.plot(data)
plt.figure(2)
plt.plot(diff_signal)
plt.figure(3)
plt.plot(cmndiff_signal)
print(len(data))


# %%
