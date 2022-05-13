import ctypes
myctest = ctypes.CDLL('./test.so')
myctest.myprint()
