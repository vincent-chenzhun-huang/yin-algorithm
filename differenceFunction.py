import numpy as np

def differenceFunction(x, N, tau_max):
    """
    Compute difference function of data x. This corresponds to equation (6) in [1]
    Fastest implementation. Use the same approach than differenceFunction_scipy.
    This solution is implemented directly with Numpy fft.
    :param x: audio data
    :param N: length of data
    :param tau_max: integration window size
    :return: difference function
    :rtype: list
    """

    x = np.array(x, np.float64)
    w = x.size
    tau_max = min(tau_max, w)
    x_cumsum = np.concatenate((np.array([0.]), (x * x).cumsum()))
    size = w + tau_max
    p2 = (size // 32).bit_length()
    nice_numbers = (16, 18, 20, 24, 25, 27, 30, 32)
    size_pad = min(x * 2 ** p2 for x in nice_numbers if x * 2 ** p2 >= size)
    fc = np.fft.rfft(x, size_pad)
    conv = np.fft.irfft(fc * fc.conjugate())[:tau_max]
    return x_cumsum[w:w - tau_max:-1] + x_cumsum[w] - x_cumsum[:tau_max] - 2 * conv

x = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 3, 4, 5, 1, 2, 3, 4]
print(differenceFunction(x, len(x), 3))


