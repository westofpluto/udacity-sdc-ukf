##########################################################################
# plotnis.py
# Script to plot NIS values 
##########################################################################
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def buckets(nisarray):
    bins=[0.0,0.35,7.82,1000.0]
    nb=len(bins)
    cnt=[0]*nb
    np=0
    for p in nisarray:
        np+=1
        for k in range(1,nb):
            if p >= bins[k-1] and p < bins[k]:
                cnt[k-1]+=1
                break
  
    if np>0: 
        for k in range(nb):
            cnt[k]/=float(np) 
    return cnt

dataRadar = np.loadtxt( "NIS_radar.txt", usecols=[0], skiprows=1 )
fracs=buckets(dataRadar)
print "Fraction between 0.0 and 0.35 = %s " % fracs[0]
print "Fraction between 0.35 and 7.82 = %s " % fracs[1]
print "Fraction above 7.82 = %s " % fracs[2]

dataLaser = np.loadtxt( "NIS_laser.txt", usecols=[0], skiprows=1 )

nisRadar = np.transpose(dataRadar)
nisLaser = np.transpose(dataLaser)
nis = [nisRadar, nisLaser]

#
# 95% confidence levels depend on degrees of freedom
#
confRadar = 7.82
confLaser = 5.99

confidences = np.array([confRadar, confLaser])

fig = plt.figure(figsize=(16,8))
axes = []
plotnum = 1
for data, confidence in zip(nis, confidences):
    subplot = 120 + plotnum
    axis = fig.add_subplot(subplot)
    axes.append( axis )
    axis.tick_params(which='both',direction='in')
    axis.set_xlabel( 'Step', labelpad=10 )
    axis.plot( np.arange(0,len(data)), data, 'r-', label="NIS" )
    axis.axhline( y=confidence, color='b',linestyle='-' ,label="95% confidence threshold" )
    axis.legend( prop={'size':20} )
    plotnum += 1

axes[0].set_title('RADAR', fontsize=20)
axes[1].set_title('LASER', fontsize=20)

plt.tight_layout()
plt.savefig( "NIS.png", bbox_inches = 'tight', dpi = 300 )

plt.show()
