from b26_toolkit.src.instruments.microwave_generator import MicrowaveGenerator
import timeit
import matplotlib.pyplot as plt
import numpy as np

def switch_frequency(mw_generator):
    mw_generator.update({'frequency': 2.7e9})
    mw_generator.update({'frequency': 2.7e9})


times = []
mw_generator = MicrowaveGenerator()
for i in range(0,10):
    print(i)
    times.append(timeit.timeit(lambda: switch_frequency(mw_generator), number = 500))
plt.plot(times, '.')

print(times)

print(np.std(np.array(times)))

plt.title('Average: {0:.2f}, Stddev: {1:.2f}'.format(np.mean(np.array(times)), np.std(np.array(times))))
plt.xlabel('Trial Number')
plt.ylabel('Execution time (s)')

plt.show()