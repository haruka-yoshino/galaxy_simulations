import pickle
import matplotlib.pyplot as plt

footprint = "/Users/harukayoshino/Documents/uvic/masters1/ngvs catalogue"

data = {}

with open(footprint+"/NGVSbounds.txt", "r") as f:
	header = f.readline().strip().split(",")
	
	for h in header:
		data[h] = []



	for line in f:
		values = line.strip().split(",")
		
		for i in range(len(header)):
			data[header[i]].append(float(values[i]))




with open(footprint+f"/footprint_dictionary.pkl", "wb") as f:
	pickle.dump(data, f)


plt.plot(data["ra"], data["dec"])
plt.show()
plt.close()
