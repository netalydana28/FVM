import ctypes 
lib = ctypes.CDLL(".\mathoper.dll")
print(lib.sum(2, 6))