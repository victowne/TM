import matplotlib as mpl
import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib import style
import os
import pandas as pd
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 

# os.chdir('c:/')
# plt.switch_backend('agg')  
# plt.style.use('ggplot') 
#rcParams
# mpl.rcParams['xtick.direction'] = 'in'
# mpl.rcParams['ytick.direction'] = 'in'
# mpl.rcParams['xtick.major.size'] = 0
# mpl.rcParams['ytick.major.size'] = 0
# xmajorLocator   = MultipleLocator(20) #将x主刻度标签设置为20的倍数  
# xmajorFormatter = FormatStrFormatter('%5.1e') #设置x轴标签文本的格式  
# xminorLocator   = MultipleLocator(5) #将x轴次刻度标签设置为5的倍数  
# ymajorLocator   = MultipleLocator(0.5) #将y轴主刻度标签设置为0.5的倍数  
# ymajorFormatter = FormatStrFormatter('%1.0e') #设置y轴标签文本的格式  
# yminorLocator   = MultipleLocator(0.1) #将此y轴次刻度标签设置为0.1的倍数 
mpl.rcParams['mathtext.default'] = 'regular'
#
font = {'family': 'Times New Roman',
         'style': 'italic',
        'weight': 'normal',
         'color': 'black', 
          'size': 18,
        }


x = np.linspace(0,100,1001)
y = np.sin(x)
z = np.cos(x)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x,y,label = 'sin(x)')
ax1.set_xlabel('where are you going',fontdict = font)
ax1.set_ylabel('fine thank you', fontdict = font)
ax1.set_title(r'This is a contour of $\psi_i$',fontdict = font)
ax2 = ax1.twinx()
ax2.plot(x,z,'r',label = 'cos(x)')
ax1.legend(loc='upper left')
ax2.legend(loc='best')
# ax1.set_xlim(0,50)
# ax1.set_xticks((1,10,15))
# ax1.set_xticklabels(('a','b','c'),fontsize=15)
# plt.autoscale(tight=True)
# 设置坐标刻度值的大小以及刻度值的字体
plt.tick_params(labelsize=23)
ax2.spines['left'].set_color('red')
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['top'].set_linewidth(2)
ax1.spines['left'].set_linewidth(2)
ax1.spines['right'].set_linewidth(2)
#设置主刻度标签的位置,标签文本的格式  
# ax.xaxis.set_major_locator(xmajorLocator)  
# ax.xaxis.set_major_formatter(xmajorFormatter)    
# ax.yaxis.set_major_locator(ymajorLocator)  
# ax.yaxis.set_major_formatter(ymajorFormatter)  
#显示次刻度标签的位置,没有标签文本  
# ax.xaxis.set_minor_locator(xminorLocator)  
# ax.yaxis.set_minor_locator(yminorLocator)  
# ax.xaxis.grid(True, which='major') #x坐标轴的网格使用主刻度  
# ax.yaxis.grid(True,linestyle= '--', which='major') #y坐标轴的网格使用次刻度 
plt.show()