import numpy as np

def generate_gaussian_symbols(k=1, ndata=5000, p=1):
    """Doc
    :param ndata: number of symbols to generate
    :param p: signal variance
    :return: k x ndata: symbols sampled from a complex gaussian with variance p

    Note that this generates a variable which is drawn from a complex gaussian distribution with variance p
    which is equivalent to a + bj with a and b sampled from a gaussian distribution with variance p/2.
    Here we first sample a and b from a gaussian with mean 0 and variance 1, by multiplying with sqrt(p)/sqrt(2)
    we obtain variance p/2 for both a and b, given that var(constant * X) = constant^2 var(X)
    """
    s = np.sqrt(p) / np.sqrt(2) * (np.random.randn(k, ndata) +
                                    1j * np.random.randn(k, ndata))
    return s.astype(np.complex64)
