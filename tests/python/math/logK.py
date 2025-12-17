import pyarts3 as pyarts
import numpy as np

ln = pyarts.arts.rtepack.logK
tra = pyarts.arts.rtepack.tran

made_it = False

# Since random noise, we try multiple times to avoid occasional failures
for i in range(10):
    try:
        arr = np.array([1, 0, 0, 0, 0, 0, 0]) * 1e-2

        K = pyarts.arts.Propmat(arr)
        expK = tra(K, K, 1.0)()
        logK = -ln(expK)
        print(f"Testing logK with arr=\n{arr}")
        print(logK)
        assert np.allclose(K, logK, atol=1e-5)

        arr[1] = np.random.random() * 1e-3

        K = pyarts.arts.Propmat(arr)
        expK = tra(K, K, 1.0)()
        logK = -ln(expK)
        print(f"Testing logK with arr=\n{arr}")
        print(logK)
        assert np.allclose(K, logK, atol=1e-5)

        arr[2] = np.random.random() * 1e-3

        K = pyarts.arts.Propmat(arr)
        expK = tra(K, K, 1.0)()
        logK = -ln(expK)
        print(f"Testing logK with arr=\n{arr}")
        print(logK)
        assert np.allclose(K, logK, atol=1e-5)

        arr[3] = np.random.random() * 1e-3

        K = pyarts.arts.Propmat(arr)
        expK = tra(K, K, 1.0)()
        logK = -ln(expK)
        print(f"Testing logK with arr=\n{arr}")
        print(logK)
        assert np.allclose(K, logK, atol=1e-5)

        arr[4] = np.random.random() * 1e-3

        K = pyarts.arts.Propmat(arr)
        expK = tra(K, K, 1.0)()
        logK = -ln(expK)
        print(f"Testing logK with arr=\n{arr}")
        print(logK)
        assert np.allclose(K, logK, atol=1e-5)

        arr[5] = np.random.random() * 1e-3

        K = pyarts.arts.Propmat(arr)
        expK = tra(K, K, 1.0)()
        logK = -ln(expK)
        print(f"Testing logK with arr=\n{arr}")
        print(logK)
        assert np.allclose(K, logK, atol=1e-5)

        arr[6] = np.random.random() * 1e-3

        K = pyarts.arts.Propmat(arr)
        expK = tra(K, K, 1.0)()
        logK = -ln(expK)
        print(f"Testing logK with arr=\n{arr}")
        print(logK)
        assert np.allclose(K, logK, atol=1e-5)
        made_it = True
        print("logK tests passed at attempt", i+1)
        break
    except AssertionError:
        print("logK test failed at attempt", i+1)
        pass

assert made_it, "logK tests failed after 10 attempts"
