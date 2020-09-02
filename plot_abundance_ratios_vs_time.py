import math
import matplotlib.pyplot as plt

filename1 = 'new123_top_only.dat'
filename2 = 'new123_bottom_only.dat'

h2_abundance_list = []
h3_abundance_list = []
he4_abundance_list = []
heh_abundance_list = []
temp_list = []
time_list = []

i = 0
with open(filename1, 'r') as f:
    next(f)
    next(f)
    for line in f:
        field = line.split()
        h2_abundance = field[3].strip()
        h2_abundance_list.append(math.log(float(h2_abundance)))
        h3_abundance = field[4].strip()
        h3_abundance_list.append(math.log(float(h3_abundance)))
        he4_abundance = field[6].strip()
        he4_abundance_list.append(math.log(float(he4_abundance)))

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

print(h2_abundance_list)
print(h3_abundance_list)
print(he4_abundance_list)
print(heh_abundance_list)
print(temp_list)

heh_he4_ratio_list = [heh/he4 for heh, he4 in zip(heh_abundance_list, he4_abundance_list)]
h3_h2_ratio_list = [h3/h2 for h2, h3 in zip(h2_abundance_list, h3_abundance_list)]

print(heh_he4_ratio_list)
print(h3_h2_ratio_list)
print(len(heh_he4_ratio_list))
print(len(temp_list))

fig, axs = plt.subplots(1,2)
fig.set_figheight(12)
fig.set_figwidth(10)

axs[0].plot(temp_list, heh_he4_ratio_list)
axs[0].set_title('HeH to He4 abundance ratios vs Temp(T9)')
axs[0].set_xlabel('temp_exp')
axs[0].set_ylabel('HeH to He4 abundance ratio')
axs[0].invert_xaxis()

axs[1].plot(temp_list, h3_h2_ratio_list)
axs[1].set_title('H3 to H2 abundance ratios vs Temp(T9)')
axs[1].set_xlabel('temp_exp')
axs[1].set_ylabel('H3 to H2 abundance ratio')
axs[1].invert_xaxis()

fig.subplots_adjust(hspace=.5)

plt.show()
