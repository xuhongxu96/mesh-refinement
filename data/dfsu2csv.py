import mikeio
import numpy

ds = mikeio.read('data/velocity.dfsu')

last_time = ds.isel([8], axis=0)
numpy.savetxt('data/velocity.txt', last_time.data[5], delimiter='\n')