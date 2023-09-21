import math
from run_method import RunMethod
from config import Configuration

RATEC = 100/math.sqrt(3)
class PowerSummation(RunMethod):
    def __init__(self,config:Configuration):
        super().__init__(config=config)
        self.busGroup = []# cac bus tuong ung o cac slack khac nhau
        for bs1 in self.param.busSlack:
            r1 = self.config.__findBusConnected__(bs1,self.config.busC,self.config.lineC)
            #
            self.busGroup.append(r1)
    #
    def __lineDirection__(self):
        ba = self.param.busSlack[:]
        lset = set()
        for ii in range(20000):
            ba2 = []
            for b1 in ba:
                for l1 in self.config.setLineHnd:
                    if l1 not in lset:
                        if b1==self.config.lineC[l1][1]:
                            d = self.config.lineC[l1][0]
                            self.config.lineC[l1][0] = self.config.lineC[l1][1]
                            self.config.lineC[l1][1] = d
                            lset.add(l1)
                            ba2.append(d)
                        elif b1==self.config.lineC[l1][0]:
                            lset.add(l1)
                            ba2.append(self.config.lineC[l1][1])
            if len(ba2)==0:
                break
            ba= ba2.copy()
    #
    def __ordCompute__(self):
        busC = dict() # connect [LineUp,[LineDown]]
        for h1 in self.config.setBusHnd:
            busC[h1] = [0,set()]
        #
        for h1,l1 in self.config.lineC.items():
            busC[l1[1]][0]= h1     # frombus
            busC[l1[0]][1].add(h1) # tobus
        #
        self.ordc,self.ordv = [],[]
        for bs1 in self.busGroup:
            busC1 = {k:v for k,v in busC.items() if k in bs1}
            #bus already
            balr = {h1:True for h1 in bs1}
            #set order
            sord = set()
            #order compute
            ordc1 = []
            for k,v in busC1.items():
                if len(v[1])==0:
                    if v[0]!=0:
                        ordc1.append(v[0])
                        sord.add(v[0])
                        balr[k]=False
            #
            for ii in range(500):
                for k,v in busC1.items():
                    if balr[k]:
                        if len(v[1]-sord)==0:
                            if k in self.param.setSlack:
                                break
                            #
                            if v[0]!=0:
                                ordc1.append(v[0])
                            sord.add(v[0])
                            balr[k]=False
            ordv1 = [ordc1[-i-1]  for i in range(len(ordc1))]
            self.ordc.append(ordc1)
            self.ordv.append(ordv1)
    #
    def __run1config__(self,fo=''):
        """
        - result (dict): {'FLAG':,'RateMax%', 'Umax[pu]','Umin[pu]','DeltaA','RateMax%'}
        - FLAG (str): 'CONVERGENCE' or 'DIVERGENCE' or 'LOOP' or 'ISLAND'
        - DeltaA: MWH
        """
        # ok run PSM
        self.__lineDirection__()
        #print(self.lineC)
        self.__ordCompute__()
        #print(self.ordc)
        #
        res = {'FLAG':'CONVERGENCE','RateMax[%]':0, 'Umax[pu]':0,'Umin[pu]':100,'DeltaA':0,'cosP':0,'cosN':0}
        # B of Line
        BUSb = {}
        for bri,v in self.param.LINEb.items():
            if bri not in self.config.lineOff:
                bfrom = self.config.lineC[bri][0]
                bto = self.config.lineC[bri][1]
                #
                if bfrom in BUSb.keys():
                    BUSb[bfrom]+=v
                else:
                    BUSb[bfrom]=v
                #
                if bto in BUSb.keys():
                    BUSb[bto]+=v
                else:
                    BUSb[bto]=v

        # Shunt
        for k1,v1 in self.param.BUSbs.items():
            if k1 not in self.config.shuntOff:
                if k1 in BUSb.keys():
                    BUSb[k1]+=v1
                else:
                    BUSb[k1]=v1
        #
        if fo:
            rB,rG,rL = super().__initcsv__(fo)
        #
        va,ra,cosP,cosN = [],[],[1],[-1]
        #
        for pi in self.param.profileID:
            res['DeltaA'] -= self.param.loadAll[pi].real
            if self.param.dgAll: 
                res['DeltaA'] += self.param.dgAll[pi].real
            sa1,va1,dia1 = dict(),dict(),dict()# for 1 profile
            for i1 in range(self.param.nSlack):# with each slack bus
                bs1 = self.param.busSlack[i1]
                ordc1 = self.ordc[i1]
                ordv1 = self.ordv[i1]
                setBusHnd1 = self.busGroup[i1]
                vbus = {h1:complex(self.param.Ubase,0) for h1 in setBusHnd1}
                vbus[bs1] = complex(self.param.genProfile[pi][bs1],0)
                #
                s0 = 0
                du,di = dict(),dict()
                for ii in range(self.param.iterMax+1):
                    sbus = {k:v for k,v in self.param.loadProfile[pi].items() if k in setBusHnd1}
                    # 
                    for k in self.config.dgOn:
                        if k in setBusHnd1:
                            sbus[k] -= self.param.dgProfile[pi][k]
                    # B of Line + Shunt
                    for k1,v1 in BUSb.items():
                        if k1 in setBusHnd1:
                            vv = abs(vbus[k1])
                            sbus[k1] += complex(0, -vv*vv*v1)
                    # cal cong suat nguoc
                    for bri in ordc1:
                        bfrom = self.config.lineC[bri][0]
                        bto = self.config.lineC[bri][1]
                        rx = self.param.LINE[bri][2]
                        #
                        du[bri] = sbus[bto].conjugate()/vbus[bto].conjugate()*rx
                        ib = abs(sbus[bto]/vbus[bto])
                        di[bri] = ib
                        ds1 = ib*ib*rx
                        #
                        if ds1.real>0.2 and ds1.real>sbus[bto].real:# if delta S is greater than S
                            return {'FLAG':'DIVERGENCE'}
                        #
                        sbus[bfrom]+=ds1+sbus[bto]
                    # cal dien ap xuoi
                    for bri in ordv1:
                        bfrom = self.config.lineC[bri][0]
                        bto = self.config.lineC[bri][1]
                        vbus[bto]=vbus[bfrom]-du[bri]
                    #
                    if abs(s0-sbus[bs1])<self.param.epsilon:
                        break
                    else:
                        s0 = sbus[bs1]
                    #
                    if ii==self.param.iterMax:
                        return {'FLAG':'DIVERGENCE'}
                # finish
                # loss P
                res['DeltaA']+=sbus[bs1].real
                # Umax[pu]/Umin[pu]
                va.extend( [abs(v) for v in vbus.values()] )
                #
                try:
                    if sbus[bs1].imag>=0:
                        cosP.append(sbus[bs1].real/abs(sbus[bs1]))
                    else:
                        cosN.append(-sbus[bs1].real/abs(sbus[bs1]))
                except:
                    pass
                # RateMax
                for bri in ordc1:
                    ra.append( di[bri]/self.param.LINE[bri][3]*RATEC )
                #
                if fo:
                    va1.update(vbus)
                    dia1.update(di)
                    sa1.update(sbus)
            #
            if fo:
                rB,rL,rG = super().__update1profile__(pi,rB,rL,rG,sa1,va1,dia1)
        #
        res = super().__update_result__(res,va,ra,cosP,cosN)
        #
        if fo:
            super().__export_profiles__(fo,rB,rL,rG,res)
        return res
