import math
import matplotlib.pyplot as plt

choice_dict = {1:'H3', 2: 'HeH'}
filename1 = 'new123_top_only.dat'
filename2 = 'new123_bottom_only.dat'

time_list = []
h3_abundance_list = []
heh_abundance_list = []
time_list = []
temp_list = []
z_list = []

with open(filename1, 'r') as f:
    next(f)
    next(f)
    for line in f:
        fields = line.split()
        temp = fields[0].strip
        temp_list.append(int(temp[-3:]))
        h3_abundance = fields[4].strip()
        h3_abundance_list.append(math.log(float(h3_abundance)))

with open(filename2, 'r') as f:
    next(f)
    next(f)
    for line in f:
        fields = line.split()
        time = fields[1].strip()
        time_exp = int(time[-3:])
        heh_abundance = fields[10].strip()
        time_list.append(math.log(float(time)))
        heh_abundance_list.append(math.log(float(heh_abundance)))
        redshift = math.sqrt((((88.2 * math.pow(10, 16))) / math.pow(10, time_exp)) - 1) - 1
        z_list.append(math.log(redshift))

print(h3_abundance_list)
print(heh_abundance_list)
print(time_list)
print(temp_list)
print(z_list)

fig, axs = plt.subplots(2, 2)
fig.set_figheight(12)
fig.set_figwidth(10)

axs[0, 0].plot(temp_list, heh_abundance_list)
axs[0, 0].set_title('HeH abundance vs Temp (T9)')
axs[0, 0].set_xlabel('log(temp)')
axs[0, 0].set_ylabel('log(HeH abundance)')

axs[0, 1].plot(z_list, heh_abundance_list)
axs[0, 1].invert_xaxis()
axs[0, 1].set_title('HeH abundance vs Redshift')
axs[0, 1].set_xlabel('log(Redshift)')
axs[0, 1].set_ylabel('log(HeH abundance)')

axs[1, 0].plot(time_list, h3_abundance_list)
axs[1, 0].set_title('H3 abundance vs Time')
axs[1, 0].set_xlabel('log(time)')
axs[1, 0].set_ylabel('log(H3 abundance)')

axs[1, 1].plot(z_list, h3_abundance_list)
axs[1, 1].invert_xaxis()
axs[1, 1].set_title('H3 abundance vs Redshift')
axs[1, 1].set_xlabel('log(Redshift)')
axs[1, 1].set_ylabel('log(H3 abundance)')

fig.subplots_adjust(hspace=.5)

plt.show()
