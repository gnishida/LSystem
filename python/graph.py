import numpy as np
import matplotlib.pyplot as plt

data = []
f = open('angle_effect.txt', 'r')
for line in f:
	rec = []
	for v in line.split(' '):
		rec.append(v)
	data.append(rec)
data = np.array(data);
f.close()

print data

plt.figure(figsize=(20, 10))
for r in range(5):
	for c in range(5):
		plt.subplot(5, 5, r * 5 + c + 1)
		#plt.xlabel('angle')
		#plt.ylabel('#branches')
		plt.plot(data[:,0], data[:, (4 - r) * 5 + c + 1])

plt.savefig('angle_effect.png')
plt.show()