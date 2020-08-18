import math
import matplotlib.pyplot as plt

filename1 = 'new123_top_only.dat'
filename2 = 'new123_bottom_only.dat'

h3_abundance_list = []
heh_abundance_list = []
temp_list = []
time_list = []

i = 0
with open(filename1, 'r') as f:
    next(f)
    next(f)
    for line in f:
        field = line.split()
        h3_abundance = field[4].strip()
        h3_abundance_list.append(math.log(float(h3_abundance)))

with open(filename2, 'r') as f:
    next(f)
    next(f)
    for line in f:
        fields = line.split()
        temp = fields[0].strip()
        temp_exp = int(temp[-3:])
        time = fields[1].strip()
        time_exp = int(time[-3:])
        heh_abundance = fields[10].strip()
        temp_list.append(float(temp_exp))
        time_list.append(float(time_exp))
        heh_abundance_list.append(math.log(float(heh_abundance)))

print(h3_abundance_list)
print(heh_abundance_list)
print(temp_list)

fig, axs = plt.subplots(2,2)
fig.set_figheight(12)
fig.set_figwidth(10)

axs[0,0].plot(temp_list, heh_abundance_list)
axs[0,0].set_title('HeH abundance vs Temp(T9)')
axs[0,0].set_xlabel('temp_exp')
axs[0,0].set_ylabel('log(HeH abundance)')
axs[0,0].invert_xaxis()

axs[0,1].plot(temp_list, h3_abundance_list)
axs[0,1].set_title('H3 abundance vs Temp(T9)')
axs[0,1].set_xlabel('temp_exp')
axs[0,1].set_ylabel('log(H3 abundance)')
axs[0,1].invert_xaxis()

axs[1,0].plot(time_list, heh_abundance_list)
axs[1,0].set_title('HeH abundance vs Time')
axs[1,0].set_xlabel('time_exp')
axs[1,0].set_ylabel('log(HeH abundance)')

axs[1,1].plot(time_list, h3_abundance_list)
axs[1,1].set_title('H3 abundance vs Time')
axs[1,1].set_xlabel('time_exp')
axs[1,1].set_ylabel('log(H3 abundance)')

fig.subplots_adjust(hspace=.5)

plt.show()
