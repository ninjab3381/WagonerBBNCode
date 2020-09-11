import math
import numpy as np
import matplotlib.pyplot as plt

filename1 = 'new123_top_only.dat'
filename2 = 'new123_bottom_only.dat'

h2_h_abundance_list = []
h3_h_abundance_list = []
he4_abundance_list = []
h_abundance_list = []
heh_h_abundance_list = []
temp_list = []
time_list = []

i = 0
with open(filename1, 'r') as f:
    next(f)
    next(f)
    for line in f:
        field = line.split()
        h2_h_abundance = field[3].strip()
        h2_h_abundance_list.append(float(h2_h_abundance))
        h3_h_abundance = field[4].strip()
        h3_h_abundance_list.append(float(h3_h_abundance))
        he4_abundance = field[6].strip()
        he4_abundance_list.append(float(he4_abundance))

with open(filename2, 'r') as f:
    next(f)
    next(f)
    for line in f:
        fields = line.split()
        temp = fields[0].strip()
        temp_exp = int(temp[-3:])
        time = fields[1].strip()
        time_exp = int(time[-3:])
        h_abundance = fields[9].strip()
        heh_h_abundance = fields[10].strip()
        temp_list.append(float(temp_exp))
        time_list.append(float(time_exp))
        h_abundance_list.append(float(h_abundance))
        heh_h_abundance_list.append(float(heh_h_abundance))

print(h2_h_abundance_list)
print(h3_h_abundance_list)
print(he4_abundance_list)
print(heh_h_abundance_list)
print(temp_list)

heh_abundance_list = [heh_h * h for heh_h, h in zip(heh_h_abundance_list, h_abundance_list)]

heh_he4_ratio_list = [math.log(heh/he4) for heh, he4 in zip(heh_abundance_list, he4_abundance_list)]
h3_h2_ratio_list = [math.log(h3_h/h2_h) for h2_h, h3_h in zip(h2_h_abundance_list, h3_h_abundance_list)]
heh_h2_ratio_list = [math.log(heh_h/h2_h) for heh_h, h2_h in zip(heh_h_abundance_list, h2_h_abundance_list)]

print(heh_he4_ratio_list)
print(h3_h2_ratio_list)
print(heh_h2_ratio_list)

fig, axs = plt.subplots(2,2)
fig.set_figheight(12)
fig.set_figwidth(10)

axs[0,0].plot(temp_list, heh_he4_ratio_list)
axs[0,0].set_title('ln(HeH/He4) abundance ratio vs Temp(T9)')
axs[0,0].set_xlabel('temp_exp')
axs[0,0].set_ylabel('ln(HeH/He4 abundance ratio)')
axs[0,0].invert_xaxis()

axs[0,1].plot(temp_list, h3_h2_ratio_list)
axs[0,1].set_title('H3/H2 abundance ratio vs Temp(T9)')
axs[0,1].set_xlabel('temp_exp')
axs[0,1].set_ylabel('ln(H3/H2 abundance ratio)')
axs[0,1].invert_xaxis()

axs[1,0].plot(temp_list, heh_h2_ratio_list)
axs[1,0].set_title('HeH/H2 abundance ratio vs Temp(T9)')
axs[1,0].set_xlabel('temp_exp')
axs[1,0].set_ylabel('ln(HeH/H2 abundance ratio)')
axs[1,0].invert_xaxis()

fig.subplots_adjust(hspace=.5)

plt.show()
