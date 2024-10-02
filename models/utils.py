import numpy as np
import numpy.random

from scipy.signal import lfilter


def getdbm(x):
    return 10 * np.log10(np.mean(abs(x) ** 2, axis=1) / 50 / 1e-3)


def setdbm(x, dnew):
    if getdbm(x) > -100:
        out = x * np.tile(10 ** ((dnew - getdbm(x)) / 20), np.shape(x), 1)
    else:
        out = x

    return out


def limiter(x, limithi=None, limit=1):
    y = x
    if limithi is None:
        i = (abs(x) > limit)
        y[i] = limit * x[i] / abs(x[i])
    else:
        i = (x < limit)
        y[i] = limit
        i = (x > limithi)
        y[i] = limithi

    p = len(i) / len(x)

    return y, p


def softlimiter(x, p: int = 1):
    # y = (atan(abs(x). ^ p)*2/pi). ^ (1/p).*exp(1i*angle(x));
    return (np.tanh(abs(x) ** p)) ** (1/p) * np.exp(1j * np.angle(x))


def rc(t, beta):
    # i2=find(1-4*beta^2*t.^2<1e-13);t(i2)=0.11;
    i2 = ((4 * beta ** 2 * t) ** 2) == 1
    t[i2] = 0.11

    pulse = (np.sinc(t) * np.cos(np.pi * beta * t) /
             (1 - 4 * beta ** 2 * t ** 2))

    pulse[i2] = 0
    return pulse


def delay(x, n: list = [1], filter_length: int = 512):
    """ delay(x, n) Delays the vector x by n steps. If filter_length is
    % omitted a length of 512 is assumed for the filter. If the
    % signal x is oversampled by a factor of 2 or more, the
    % filter_length can be set to 5, making the function
    % considerably faster.
    % If n is omitted it is set to 1. If n is a vector, the output
    % y is a matrix where the columns are delayed x, with the
    % different delays in n.
    % x can also be a matrix, in which case y is a matrix with
    % size(x, 2)*length(n) columns. The first size(x, 2) columns is
    % then x delayed by the first element of n, and so on.
    % Updated 2017-07-06: bug found. Changed from t+frac to to t-frac
    % Updated 2018-12-23: allowing for vector of delays n
    % Updated 2020-05-12: allowing for matrix x
    % Updated 2021-03-22: bug found, in vector n operation.
    % Updated 2023-08-15: The filter length can now be specified as
    % an input parameter.
    % (c) Thomas Eriksson 2017-2023
    % Time: The function can delay ~ 7e6 complex samples per
    % second, with FilterLength = 512.
    """

    cols = len(x[0])

    y = np.zeros((len(x), len(n) * cols))
    for ind in n:
        nint = round(ind)
        y1 = np.roll(x, nint)
        if ind >= 0:
            y1[0:nint+1, :] = np.zeros((nint, cols))
        else:
            y1[len(y1)+nint:, :] = np.zeros((abs(nint), cols))

        frac = ind - nint
        # Check if integer
        if frac != 0:
            # Simple version if there is some oversampling
            if filter_length < 10:
                h = rc(np.transpose(
                    range(-filter_length, filter_length)) - frac, 0.5)
            else:
                h = np.sinc(np.transpose(
                    range(-filter_length, filter_length) - frac) * np.hanning(2 * filter_length + 1))

            y1 = [[y1], [np.zeros((filter_length, cols))]]
            y1 = lfilter(h, 1, y1)
            y1 = y1[len(y1)-len(x):, :]

        y[:, (ind-1)*cols+1:ind*cols] = y1


def randn_c(rows: int = 1, cols: int = 1, threshold: int = 0):
    """Complex normal random numbers, neg"""
    if threshold <= 0:
        x = np.sqrt(threshold ** 2 - 2 * np.log(numpy.random.normal(rows, cols))
                    ) * np.exp(1j * 2 * np.pi * np.random.normal(rows, cols))
    else:
        x = np.sqrt(-2 * np.log(1 - (1 - np.exp(-threshold ** 2 / 2)) * numpy.random.normal(
            rows, cols))) * np.exp(1j * 2 * np.pi * numpy.random.normal(rows, cols))
