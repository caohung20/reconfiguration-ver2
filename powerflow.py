import math
import time
from powersummation import PowerSummation
from newton_raphson import NewtonRaphson
from config import Configuration
from sparse_matrix_lfnewton import SparseNewtonRaphson
class PowerFlow:
    def __init__(self,config: Configuration):
        self.config = config
        self.param = self.config.param
    
    def run1Config_WithObjective(self,option=None,fo=''):
        #
        v1 = self.run1Config(fo)
        if v1['FLAG']!='CONVERGENCE':
            obj = math.inf
        else: #RateMax[%]    Umax[pu]    Umin[pu]    Algo_PF    option_PF
            obj = v1['DeltaA']
            # constraint
            obj+=self.param.setting['RateMax[%]'][1]*max(0, v1['RateMax[%]']-self.param.setting['RateMax[%]'][0])
            obj+=self.param.setting['Umax[pu]'][1]*max(0, v1['Umax[pu]']-self.param.setting['Umax[pu]'][0])
            obj+=self.param.setting['Umin[pu]'][1]*max(0,-v1['Umin[pu]']+self.param.setting['Umin[pu]'][0])
            #cosP ycau cosP>0.9
            obj+=self.param.setting['cosPhiP'][1]*max(0,-v1['cosP']+self.param.setting['cosPhiP'][0])
            #cosN ycau cosN<-0.95
            obj+=self.param.setting['cosPhiN'][1]*max(0,v1['cosN']-self.param.setting['cosPhiN'][0])
        #
        #self.config.lineOff.sort()
        #self.config.shuntOff.sort()
        v1['Objective'] = obj
        v1['LineOff'] = self.config.lineOff
        v1['ShuntOff'] = self.config.shuntOff
        return v1
    
    def run1Config(self,fo=''):
        # check loop island and return if loop or island appears 
        t0 = time.time()
        #
        """ run PF 1 config """
        if self.param.setting['Algo_PF']=='PSM':
            return PowerSummation(self.config).__run1config__(fo)  
        elif self.param.setting['Algo_PF']=='N-R':
            return NewtonRaphson(self.config).__run1config__(fo) 
        elif self.param.setting['Algo_PF'] == 'SNR':
            return SparseNewtonRaphson(self.config).__run1config__(fo)
        return None

