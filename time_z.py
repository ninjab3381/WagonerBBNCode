import math

time_exp = 16
redshift = math.sqrt((((88.2 * math.pow(10, 16))) / math.pow(10, time_exp)) - 1) - 1
print(redshift)
print(math.log(redshift))