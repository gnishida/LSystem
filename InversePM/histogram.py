import numpy as np
import matplotlib.pyplot as plt

data = []
f = open('pm_linearity.txt', 'r')
for line in f:
	data.append(float(line))

f.close()
plt.hist(data, bins=30)
plt.savefig('pm_histogram.png')
plt.show()