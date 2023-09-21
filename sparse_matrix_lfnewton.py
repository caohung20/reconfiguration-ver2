import numpy as np
import math
from ybus import YbusMatrix
from run_method import RunMethod
from config import Configuration
class SparseNewtonRaphson(RunMethod):
    def __init__(self,config:Configuration):
        super().__init__(config=config)
    def __run1config__(self, fo=''):
        # ready to run Newton raphson
        Ybus = YbusMatrix(self.config).__get_sparse_Ybus__()
        #
        # Measure memory usage of the variable
        #memory_usage = sys.getsizeof(Ybus)
        # Print the memory usage
        #print(f"Memory usage of the variable: {memory_usage} bytes")
        #
        res = {'FLAG':'CONVERGENCE','RateMax[%]':0, 'Umax[pu]':0,'Umin[pu]':100,'DeltaA':0,'cosP':0,'cosN':0}
        #
        if fo:
            rB,rG,rL = super().__initcsv__(fo)
        Ym, theta = dict(),dict()
        for k,v in Ybus.items():
            Ym[k] = abs(v)
            theta[k] = np.angle(v,deg=False)
        #
        countSlack = 0
        countPV = 0
        slackCounted = []
        PVcounted = []
        for bi in self.config.setBusHnd:
            if self.param.BUS[bi][3] == 3:
                countSlack+=1 
            elif self.param.BUS[bi][3] == 2:
                countPV += 1 
            slackCounted.append(countSlack)
            PVcounted.append(countPV)
        no_jacobi_equation = 2 * self.param.nBus - countPV - 2 * self.param.nSlack         
        #
        DC = np.zeros(no_jacobi_equation) 
        #
        va,ra,cosP,cosN = [],[],[1],[-1]
        for pi in self.param.profileID:
            # initialize voltage magnitude of all buses
            # iterate in setBUShnd in case some buses are off 
            Vm = {bi:float(self.param.Ubase) for bi in self.config.setBusHnd}
            # initialize voltage angle 
            delta = {bi:0 for bi in self.config.setBusHnd}
            P = {bi:(-self.param.loadProfile[pi][bi]).real for bi in self.config.setBusHnd}    #+self.Pgen
            Q = {bi:(-self.param.loadProfile[pi][bi]).imag for bi in self.config.setBusHnd}    #+self.Qgen
            #update P and Q of distributed generation in buses
            for bi in self.config.dgOn:
                P[bi] += self.param.dgProfile[pi][bi].real
                Q[bi] += self.param.dgProfile[pi][bi].imag
            # update Vm of slack buses
            for bs in self.param.setSlack:
                Vm[bs] = self.param.genProfile[pi][bs]
            for ii in range(self.param.iterMax+1):
                # Initialize Jacobian Matrix
                A = np.zeros((no_jacobi_equation,no_jacobi_equation))
                for b1 in self.config.setBusHnd:
                    #
                    J1_row_offdiag_idx = J2_row_offdiag_idx = int((b1 - slackCounted[b1-1])-1)
                    J1_diag_idx = J2_row_diag_idx = J3_col_diag_idx = J1_row_offdiag_idx 
                    #
                    J3_row_offdiag_idx = J4_row_offdiag_idx = int((self.param.nBus+b1-slackCounted[b1-1]-PVcounted[b1-1]-self.param.nSlack)-1)
                    J4_diag_idx = J2_col_diag_idx = J3_row_diag_idx = J3_row_offdiag_idx 
                    #
                    J11 = 0
                    J22 = 0
                    J33 = 0
                    J44 = 0
                    #
                    for li in self.config.busC[b1]:
                        for b2 in self.config.lineC[li]:
                            if b1 == b2:
                                pass
                            else:
                                line = (b1,b2)
                                # diagonal elements of J1
                                J11 += Vm[b1] * Vm[b2] * Ym[line] * math.sin((theta[line] - delta[b1] + delta[b2]).real)
                                # diagonal elements of J3          
                                J33 += Vm[b1] * Vm[b2] * Ym[line] * math.cos((theta[line] - delta[b1] + delta[b2]).real)
                                if self.param.BUS[b1][3] != 3:
                                    #slackbus doesn't exist in J2 and J4
                                    #diagonal elements of J2
                                    J22 += Vm[b2] * Ym[line] * math.cos((theta[line] - delta[b1] + delta[b2]).real)
                                    # diagonal elements of J4               
                                    J44 += Vm[b2] * Ym[line] * math.sin((theta[line] - delta[b1] + delta[b2]).real)
                                if self.param.BUS[b1][3] != 3 and self.param.BUS[b2][3] != 3:
                                    J1_col_offdiag_idx = J3_col_offdiag_idx = int(b2 - slackCounted[b2-1] - 1)
                                    J2_col_offdiag_idx = J4_col_offdiag_idx = int(self.param.nBus + b2 - PVcounted[b2-1] - slackCounted[b2-1] - self.param.nSlack - 1)
                                    # off diagonal elements of J1
                                    A[J1_row_offdiag_idx][J1_col_offdiag_idx] = -Vm[b1] * Vm[b2] * Ym[line] * math.sin((theta[line] - delta[b1] + delta[b2]).real)
                                
                                    if self.param.BUS[b2][3] == 1:
                                        # off diagonal elements of J2
                                        A[J2_row_offdiag_idx][J2_col_offdiag_idx] = Vm[b1] * Ym[line] * math.cos((theta[line] - delta[b1] + delta[b2]).real)
                                
                                    if self.param.BUS[b1][3] == 1:
                                        # off diagonal elements of J3
                                        A[J3_row_offdiag_idx][J3_col_offdiag_idx] = -Vm[b1] * Vm[b2] * Ym[line] * math.cos((theta[line] - delta[b1] + delta[b2]).real)
                                
                                    if self.param.BUS[b1][3] == 1 and self.param.BUS[b2][3] == 1:
                                        # off diagonal elements of J4
                                        A[J4_row_offdiag_idx][J4_col_offdiag_idx] = -Vm[b1] * Ym[line] * math.sin((theta[line] - delta[b1] + delta[b2]).real)
                    #   
                    Pk = Vm[b1]**2 * Ym[b1] * math.cos((theta[b1])) + J33
                    Qk = -Vm[b1]**2 * Ym[b1] * math.sin((theta[b1])) - J11
                    if self.param.BUS[b1][3] == 3:
                        # swing bus
                        P[b1] = Pk
                        Q[b1] = Qk
                    if self.param.BUS[b1][3] == 2:
                        Q[b1] = Qk
                        """
                        if Qmax[b1-1] != 0:
                            Qgc = Q[b1-1] + Qd[n-1] - Qsh[n-1]
                            if ii <= 7:                   #Between the 2th & 6th iterations
                                if ii > 2:                #the Mvar of generator buses are
                                    if Qgc < Qmin[n-1]:   #tested. If not within limits Vm(n)
                                        Vm[n-1] += 0.01   #is changed in steps of 0.01 pu to
                                    elif Qgc > Qmax[n-1]: #bring the generator Mvar within
                                        Vm[n-1] -= 0.01   #the specified limits.            
                        """
                    if self.param.BUS[b1][3] != 3:
                        # diagonal elements of J1
                        A[J1_diag_idx][J1_diag_idx] = J11        
                        DC[J1_diag_idx] = P[b1] - Pk
                    if self.param.BUS[b1][3] == 1:
                        # diagonal elements of J2
                        A[J2_row_diag_idx][J2_col_diag_idx] = 2 * Vm[b1] * Ym[b1] * math.cos(theta[b1]) + J22    
                        # diagonal elements of J3
                        A[J3_row_diag_idx][J3_col_diag_idx] = J33
                        # diagonal elements of J4         
                        A[J4_diag_idx][J4_diag_idx] = -2 * Vm[b1] * Ym[b1] * math.sin(theta[b1]) - J44   
                        DC[J4_diag_idx] = Q[b1] - Qk
                #matrix A left division Matrix DC
                DX = np.linalg.solve(A, DC.T)
                for bi in self.config.setBusHnd:
                    del_update_DXidx = int(bi - slackCounted[bi-1]-1)
                    Vm_update_DXidx = int(self.param.nBus + bi - PVcounted[bi-1] - slackCounted[bi-1] - self.param.nSlack-1)
                    if self.param.BUS[bi][3] != 3:
                        delta[bi] +=  DX[del_update_DXidx]
                    if self.param.BUS[bi][3] == 1:
                        Vm[bi] += DX[Vm_update_DXidx]
                maxerror = max(abs(DC))
                if maxerror < self.param.epsilon:
                    break
                if ii==self.param.iterMax:
                    return {'FLAG':'DIVERGENCE'}
            # finish NR
            # store all voltage magnitude value in all profile
            va.extend(Vm.values())
            vbus = {bi:(Vm[bi]* np.cos(delta[bi])+Vm[bi]*1j*np.sin(delta[bi])) for bi in self.config.setBusHnd }
            sbus = {bi:(P[bi] + Q[bi]*1j) for bi in self.config.setBusHnd}
            for bi in self.param.busSlack:
                sbus[bi] += self.param.loadProfile[pi][bi]
            """
            for bi in self.busPV:
                # update Q load in PV bus
                sbus[bi] += self.loadProfile[pi][bi].imag * 1j 
            """
            slt = 0
            Il = dict()
            for li,bi in self.config.lineC.items():
                line = (bi[0],bi[1])
                Ib1 = (vbus[bi[0]]-vbus[bi[1]])*(-Ybus[line])
                Ib2 = (vbus[bi[1]]-vbus[bi[0]])*(-Ybus[line])
                #
                if abs(Ib1) >= abs(Ib2):
                    Il[li] = abs(Ib1)
                else: 
                    Il[li] = abs(Ib2)   
                rate = ((Il[li])/self.param.LINE[li][3])*RATEC
                ra.append(rate)
                Snk = vbus[bi[0]]*np.conj(Ib1)
                Skn = vbus[bi[1]]*np.conj(Ib2)
                slt += Snk +Skn
            for bs in self.param.busSlack:
                if sbus[bs].imag:
                    cosP.append(sbus[bs].real/abs(sbus[bs]))
                else:
                    cosN.append(-sbus[bs].real/abs(sbus[bs])) 
            if fo:
                sa1,va1,dia1 = dict(),dict(),dict()# for 1 profile
                va1.update(vbus)
                dia1.update(Il)
                sa1.update(sbus)
                rB,rL,rG =super().__update1profile__(pi,rB,rL,rG,sa1,va1,dia1)
            res['DeltaA'] += slt.real
        res = super().__update_result__(res,va,ra,cosP,cosN)
        if fo:
            super().__export_profiles__(fo,rB,rL,rG,res)
        return res