def BusidincodeMapFData(BusNameinCode,BusNameinData):
  import numpy as np
  BusNameinData=BusNameinData-1
  busid=np.array(BusNameinCode).searchsorted(BusNameinData)
  return busid