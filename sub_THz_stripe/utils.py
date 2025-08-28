import numpy as np
import numpy.random

from scipy.signal import lfilter

global_seed = None


def getdbm(x):
    return 10 * np.log10(np.mean(abs(x) ** 2) / 50 / 1e-3)


def setdbm(x, dnew):
    if getdbm(x) > -100:
        out = x * np.tile(10 ** ((dnew - getdbm(x)) / 20), len(x))
    else:
        out = x

    return out


def db_to_power(db_value: float | int):
    """Convert a value in decibels (dB) to a power value.

    The power/decibel relationship is as follows:
    db_value = 10 * log10(power_value)
    """
    return 10.0 ** (db_value / 10.0)


def power_to_db(power_value: float | int):
    """Convert a power value to decibels (dB).

    The power/decibel relationship is as follows:
    db_value = 10 * log10(power_value)
    """
    return 10 * np.log10(power_value)


def db_to_magnitude(db_value: float | int):
    """Convert a value in decibels (dB) to a magnitude value.

    The power/decibel relationship is as follows:
    db_value = 20 * log10(power_value)
    """
    return 10.0 ** (db_value / 20.0)


def limiter(x, limit=1, limithi=None):
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
    """Delays the vector x by n steps.

    If filter_length is omitted a length of 512 is assumed for the filter.
    If x is oversampled by a factor of 2 or more, the filter_length can be set to 5,
    making the function considerably faster. Omitting sets it to 1.

    When n is a vector, the output is a matrix where the rows are delayed x, with
    the different delays in n.
    When x is a matrix the output is a matrix with size(x, 2)*length(n) columns. The
    first size(x, 2) columns is then x delayed by the first element of n, and so on.

    (c) Thomas Eriksson 2017-2023

    Time: The function can delay ~ 7e6 complex samples per second, with FilterLength = 512.
    """

    # In case the input is a matrix the output matrix has the same amount of columns but (x.rows * len(n)) rows.
    # Every x.rows of the output matrix forms one delayed version of the input.
    shape = (x.shape[0] * len(n), x.shape[1])
    rows = x.shape[0]

    # Zero array/matrix used to build up the output.
    y = np.zeros(shape=shape, dtype=np.complex128)

    # Iterate over all the given delays.
    for i, delayn in enumerate(n):
        # Round the delay so we can use it as index.
        delay_int = round(delayn)
        # Roll the input with the given delay.
        y1 = np.roll(x, delay_int)

        # If the delay is greater than zero add zeros to the beginning of the array.
        # When it is smaller than zero add zeros at the end of the array.
        if delayn >= 0:
            y1[:, :delay_int] = np.zeros((rows, delay_int))
        else:
            y1[:, y1.shape[1] +
                delay_int:] = np.zeros((rows, abs(delay_int)))

        frac = delayn - delay_int
        # Check if integer
        if frac != 0:
            # Simple version if there is some oversampling
            if filter_length < 10:
                h = rc(np.arange(-filter_length, filter_length + 1) - frac, 0.5)
            else:
                h = (np.sinc(np.arange(-filter_length, filter_length +
                     1) - frac) * np.hanning(2 * filter_length + 1))

            y1 = np.concat([y1, np.zeros((rows, filter_length))], axis=1)
            y1 = lfilter(h, 1, y1)
            y1 = y1[:, y1.shape[1]-x.shape[1]:]

        y[(i*rows):((i+1) * rows), :] = y1

    return y


def randn_c(rows: int = 1, cols: int = 1, threshold: int = 0):
    """Complex normal random numbers, neg"""
    if rows == 1:
        rand1 = numpy.random.uniform(size=cols)
        rand2 = numpy.random.uniform(size=cols)
    else:
        rand1 = numpy.random.uniform(size=(rows, cols))
        rand2 = numpy.random.uniform(size=(rows, cols))
    if threshold <= 0:
        x = (np.sqrt(threshold ** 2 - 2 * np.log(rand1))
             * np.exp(1j * 2 * np.pi * rand2))
    else:
        x = (np.sqrt(-2 * np.log(1 - (1 - np.exp(-threshold ** 2 / 2))
             * rand1)) * np.exp(1j * 2 * np.pi * rand2))

    return x


def makespiralconstellation(m, f):
    m = np.arange(1, m+1)
    ang = np.sqrt((4 * np.pi * m) ** 2 * f / 2 +
                  np.sqrt(((4 * np.pi * m) ** 2 * f / 2) ** 2 + (4 * np.pi * m) ** 2))

    return ang * np.exp(1j * ang)


def randconst(rows, cols, m=16, type: str = 'QAM'):
    """RANDCONST Complex constellation."""
    match type:
        case 'QAM':
            match m:
                case 2:
                    c = np.array([-1, 1])
                case 8:
                    c = np.array([[-3, -1, 1, 3], [-3, -1, 1, 3]] +
                                 1j * np.array(list(np.ones((1, 4))) +
                                               list(-np.ones((1, 4)))))
                    c = c.flatten('F')
                case 32:
                    xpoints = np.arange(-5, 6, 2)
                    ypoints = np.arange(-3, 4, 2)
                    x, y = np.meshgrid(xpoints, ypoints)
                    x = x.flatten('F')
                    y = y.flatten('F')
                    c = x + 1j * y
                    arr = np.arange(-3, 4, 2)
                    c = np.array(list(arr - 1j * 5) +
                                 list(c) + list(arr + 1j * 5))
                case 128:
                    xpoints = np.arange(-11, 12, 2)
                    ypoints = np.arange(-7, 8, 2)
                    x, y = np.meshgrid(xpoints, ypoints)
                    x = x.flatten('F')
                    y = y.flatten('F')
                    c = x + 1j * y
                    arr = np.arange(-7, 8, 2)
                    c = np.array(list(arr - 1j * 11) + list(arr - 1j * 9) +
                                 list(c) + list(arr + 1j * 9) + list(arr + 1j * 11))
                case 512:
                    xpoints = np.arange(-15, 16, 2)
                    ypoints = np.arange(-23, 24, 2)
                    x, y = np.meshgrid(xpoints, ypoints)
                    x = x.flatten('F')
                    y = y.flatten('F')
                    xpoints = np.arange(-3, 4, 2)
                    ypoints = np.arange(-15, 16, 2)
                    x2, y2 = np.meshgrid(xpoints, ypoints)
                    x2 = x2.flatten('F')
                    y2 = y2.flatten('F')
                    c = np.array(list(x + 1j * y) + list(x2 + 1j *
                                                         y2 - 20) + list(x2 + 1j * y2 + 20))
                case _:
                    q = np.log2(m)
                    if (q % 2) != 0:
                        raise ValueError('Bad constellation size.')
                    q = round(np.sqrt(m))
                    r = np.arange(1, q + 1) - (q + 1) / 2
                    c = np.reshape(np.tile(np.transpose(
                        [r]), (1, q)) + 1j * np.tile(r, (q, 1)), (q ** 2,))
        case 'PSK':
            c = np.exp(1j * 2 * np.pi * np.arange(1, m+1) / m)
        case 'SPIRAL':
            c = makespiralconstellation(m, 0)
        case _:
            raise ValueError(
                'Wrong constellation type. Choose QAM, PSK or SPIRAL')

    c = np.sqrt(2) / np.std(c) * c

    rng = numpy.random.default_rng(seed=global_seed)
    i = rng.integers(len(c), size=(rows, cols))
    x = c[i]

    return x, c


def rrc(t, beta):
    i1 = (t == 0)
    i2 = ((4 * beta * t) ** 2 == 1)

    t[i1] = 0.11  # value 0.11 not used
    t[i2] = 0.11
    pulse = (np.sin(np.pi * t * (1 - beta)) + 4 * beta * t *
             np.cos(np.pi * t * (1 + beta))) / (np.pi * t * (1 - (4 * beta * t) ** 2))

    pulse[i1] = (1 - beta + 4 * beta / np.pi)
    pulse[i2] = beta / np.sqrt(2) * ((1 + 2 / np.pi) * np.sin(
        np.pi / 4 / beta) + (1 - 2/np.pi) * np.cos(np.pi / 4 / beta))

    return pulse


def pulseshape(x, oversampling=5, beta=0.07, fl: int = 25):
    """
    % pulseshape(x,oversampling,beta)
    % Root-raised-cosine pulse shaping.
    % x is a vecor or matrix where the columns are samples to pulseshape.
    % beta is the rootraisedcosine factor.
    % FL is the half filterlength, in symbols. Default 25.
    % h is pulse shape. Default rrc.
    %
    % The function oversamples and filters with rrc(beta) filter.
    % The resulting vector contains extra samples for filling and emptying the
    % filters in the pulseshape and matched filter functions. Remove
    % 2*oversampling*FL samples from the beginning of the vector
    """

    if fl < round(1 / beta):
        raise ValueError('Warning: FilterLength too small to guarantee 35 dB.')

    if type(oversampling) is int:
        pulse_filter = rrc(
            np.arange(0, fl + (1/oversampling), 1/oversampling), beta)
        flipped_filter = np.flip(pulse_filter[1:])
        # this makes sure that t = 0 is represented
        pulse_filter = np.concatenate((flipped_filter, pulse_filter))

        # append zeros for filter delay
        shape = (x.shape[0], fl)
        xz = np.zeros(shape=shape, dtype=np.dtype('complex128'))
        xz = np.concatenate((x, xz), axis=1)

        # multiphase filtering. Same as zeropadding and filtering.
        shape = (xz.shape[0], xz.shape[1] * oversampling)
        xps = np.zeros(shape=shape, dtype=np.dtype('complex128'))
        for d in range(oversampling):
            l2 = lfilter(pulse_filter[d::oversampling], 1, xz)
            xps[:, d::oversampling] = l2
    else:
        zs = np.zeros(shape=(1, fl), dtype=np.dtype('complex128'))
        x = np.concatenate((zs, x), axis=1)
        # t = linspace(1, length(x)+1-1/oversampling-0.000001, round(length(x)*oversampling))
        t = np.arange(1, x.shape[1] + 1 - 0.001, 1/oversampling)
        x = np.concatenate(
            (np.zeros((1, fl)), x, np.zeros((1, fl + 1))), axis=1)

        t = t + fl
        xps = np.zeros(shape=(x.shape[0], len(t)),
                       dtype=np.dtype('complex128'))

        p1 = np.round(t - fl).astype(int)
        p2 = np.round(t + fl).astype(int)

        for k in range(len(t)):
            t2 = np.arange(p1[k], p2[k] + 1) - t[k] + np.finfo(float).eps
            a = np.array([np.sin(np.pi * t2 * (1 - beta)) + 4 *
                          beta * t2 * np.cos(np.pi * t2 * (1 + beta))])
            b = np.array([np.pi * t2 * (1 - (4 * beta * t2) ** 2)])
            c = np.transpose(x[:, p1[k]-1:p2[k]])
            xps[:, k] = np.matmul((a / b), c)

    return xps
