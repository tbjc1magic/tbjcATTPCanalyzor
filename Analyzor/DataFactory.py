from sqlalchemy import *
import numpy as np
import pandas as pd
import cv2
import time
import re

class DataFactory(object):

    def __init__(self,data_path,map_path):

        ##### set single thread computation #########
        cv2.setNumThreads(1)
        self.fname = data_path.split('/')[-1][:-3]

        print self.fname+" initialization will take some time"
        start_time = time.time()
        engine = create_engine('sqlite+pysqlite:///'+data_path)
        self.ADCdf = pd.io.sql.read_sql("SELECT * FROM ADC limit 10000", engine)
        end_time = time.time()

        engine = create_engine('sqlite+pysqlite:///'+map_path)
        self.ProtoMapdf = pd.io.sql.read_sql("SELECT * FROM ProtoMap", engine)
        print  self.fname+"loading finished"

    def InitT3(self):

        print self.fname+" begin processing"
        start_time = time.time()
        ADCdfn = self.ADCdf.copy()
        ProtoMapdf = self.ProtoMapdf

        ADCdfn.columns = [np.uint16(_[1:]) if re.match('t\d+',_) is not None else str(_) for _ in ADCdfn.columns ]
        ADCdfn['max'] = ADCdfn.iloc[:,3+50:-50].max(axis=1)
        mask= (ADCdfn.iloc[:,3:-1]>20) & (ADCdfn.iloc[:,3:-1].gt(ADCdfn['max']*0.2,axis=0))
        ADCdfn.iloc[:,3:-1] = ADCdfn.iloc[:,3:-1][mask].fillna(0)
        end_time = time.time()

        start_time = time.time()
        n1 = [_ for _ in ADCdfn.columns if type(_) is np.uint16]
        n2 = [_ for _ in ADCdfn.columns if type(_) is not np.uint16]
        t2 = pd.melt(ADCdfn.iloc[:],id_vars=n2,value_vars=n1).drop(['ID'],axis=1)
        t2.columns = ['EventID','PadNum','max','time','charge']
        self.t3 = pd.merge(t2[t2['charge']>20],ProtoMapdf,on='PadNum')[['EventID','PadNum','time','PadX','PadY','charge']]

        print self.fname+" processing finished"

        self.ADCdfn = ADCdfn

    def ConstructImage(self,EID):

        t3 = self.t3

        tmp = t3[(t3['EventID']==EID)&(t3['charge']>3)].copy()
        p = (tmp['PadNum']-1)%63+1

        tmp['PadPos'] = ((p-9)*(p>9)+p)*(tmp['PadNum']!=0)

        Q = [tmp[(tmp['PadX']>0)&(tmp['PadY']>0)] ,tmp[(tmp['PadX']<0)&(tmp['PadY']>0)],
            tmp[(tmp['PadX']<0)&(tmp['PadY']<0)], tmp[(tmp['PadX']>0)&(tmp['PadY']<0)]]

        image1 = np.zeros([300,600])

        image1[(-Q[0]['PadPos'].values+151).astype(np.int), Q[0]['time'].values.astype(np.int)] =255
        image1[(Q[2]['PadPos'].values+150).astype(np.int), Q[2]['time'].values.astype(np.int)] =255

        image2 = np.zeros([300,600])

        image2[(-Q[1]['PadPos'].values+151).astype(np.int), Q[1]['time'].values.astype(np.int)] =255
        image2[(Q[3]['PadPos'].values+150).astype(np.int), Q[3]['time'].values.astype(np.int)] =255

        width,height = 300,600
        a = np.concatenate([np.arange(width/2,0,-1),np.arange(1,width/2+1,1)])
        weights = np.tile(a,height).reshape(height,width).T/150.0

        image1[:50,:] = 0
        image2[:50,:] = 0

        image1[-50:,:] = 0
        image2[-50:,:] = 0

        if np.sum(image1*image2*weights**4)>np.sum(image1*image2[::-1]*weights**4):
            image = image1+image2#[::-1]
        else:
            image = image1+image2[::-1]

        image = np.where(image>100,255,0).astype(np.uint8)

        gray = cv2.GaussianBlur(image, (3, 3), 0)
        ret,im = cv2.threshold(gray.astype(np.uint8), 10, 255, cv2.THRESH_BINARY)

        thresh = im.astype(np.uint8)
        for _ in range(2):
            thresh = cv2.erode(thresh, None, iterations=1)
        for _ in range(2):
            thresh = cv2.dilate(thresh, None, iterations=1)

        return thresh

    def InitMesh(self):
        def Process(Event):
            mesh = Event.iloc[:,3:][Event.iloc[:,3:]>20].sum(axis=0)
            return mesh
        self.mesh_df = self.ADCdf.groupby('EventID').apply(Process)

    def plotMap(self):
        ProtoMapdf = self.ProtoMapdf
        for row in ProtoMapdf.iloc[:252].iterrows():
            plt.scatter(row[1]['PadX'],row[1]['PadY'],marker='${}$'.format(row[1]['PadNum']),s=200)
