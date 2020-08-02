# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 19:04:14 2020

@author: JingCe
"""

import pandas as pd
import datetime
import time
import numpy as np
import matplotlib.pyplot as plt
import enu_tools as tools


def UTC2GPST(date):
    if date.month <= 2:
        date.year = date.year - 1
        date.month = date.month + 12
    julday = int(365.25 * date.year) + int(30.6001 * (date.month + 1)) + date.day + 1720981.5 + date.hour / 24.0 + \
                 date.minute/1440.0 + date.second/86400.0
    
    gps_week = int((julday - 2444244.5)/7)
    day_of_week = int((julday - 2444244.5) % 7)
    second_of_week = 24 * 60 * 60 * day_of_week + date.hour * 60 * 60 + date.minute * 60 + date.second
    # print(gps_week, day_of_week, second_of_week)
    return [gps_week, second_of_week ,day_of_week]

#-------------------导入定位结果文件(ECEF)-------------------#
# data_frame=pd.read_csv('rtk1_CSZ1.csv')
data_frame=pd.read_csv('cusv35501.csv')

# print(len(data_frame))

#-------------------对时间字段格式化-------------------#
for ii in range(0,len(data_frame)): 
    date_str=data_frame.at[ii,'Calendar']+' '+data_frame.at[ii,'GPST']
    date_t=datetime.datetime.strptime(date_str,'%Y/%m/%d %H:%M:%S')
    date_gpst=UTC2GPST(date_t)
    data_frame.at[ii,'GPS WEEK']=date_gpst[0]
    data_frame.at[ii,'SOW']=date_gpst[1]
    data_frame.at[ii,'SOWW']=date_gpst[1]

#-------------------构造连续时间序列-------------------#
t_beg=data_frame.at[0,'SOW']
t_end=data_frame.at[len(data_frame)-1,'SOW']
t_list= np.arange(t_beg,t_end,30.0)

#-------------------修改Dataframe索引字段-------------------#
data_frame.set_index(['SOW'],inplace=True)
data_frame.rename(columns={'SOWW':'SOW'},inplace=True)

#-------------------提取XYZ坐标字段-------------------#
x_mean=data_frame['x-ecef(m)'].mean()
y_mean=data_frame['y-ecef(m)'].mean()
z_mean=data_frame['z-ecef(m)'].mean()
xyz_mean_array=np.array([x_mean,y_mean,z_mean])

llh_mean=tools.ecef2llh([x_mean,y_mean,z_mean])
lat_mean_deg=llh_mean[0]/np.pi*180.0
lon_mean_deg=llh_mean[1]/np.pi*180.0

#-------------------提取中断时间的历元赋Nan-------------------#
for ii in t_list:
    if (ii in data_frame.index):
        xyz_array=np.array([
            data_frame.at[ii,'x-ecef(m)'],
            data_frame.at[ii,'y-ecef(m)'],
            data_frame.at[ii,'z-ecef(m)']
                           ])
        
        if(data_frame.at[ii,'Q']==5):
            enu_t=tools.ecef2enu(xyz_array,xyz_mean_array,lat_mean_deg,lon_mean_deg)
            data_frame.at[ii,'E']=enu_t[0]
            data_frame.at[ii,'N']=enu_t[1]
            data_frame.at[ii,'U']=enu_t[2]
    else:
        data_frame.at[ii,'Q']=0


#-------------------重新构造enu Dataframe-------------------#
enu_frame=data_frame[['SOW','E','N','U']]

#-------------------pandas画图-------------------#
# enu_frame.plot(subplots=True,kind='line',x='SOW',y='E',figsize=(10, 6),color='red',grid=True)
# enu_frame.plot(subplots=True,kind='line',x='SOW',y='N',figsize=(10, 6))

# enu_frame.plot.scatter(x='SOWW',y='E')
# enu_frame.plot.scatter(x='SOWW',y='N')
# enu_frame.plot.scatter(x='SOWW',y='U')

#-------------------matplotlib画图-------------------#
fig=plt.figure(figsize=(10,6),dpi=80)

ax=fig.add_subplot(3,1,1)
ax.scatter(enu_frame['SOW'],enu_frame['E'],s=6)
# ax.plot(enu_frame['E'])
ax=fig.add_subplot(3,1,2)
ax.scatter(enu_frame['SOW'],enu_frame['N'],s=6)
ax=fig.add_subplot(3,1,3)
ax.scatter(enu_frame['SOW'],enu_frame['U'],s=6)

