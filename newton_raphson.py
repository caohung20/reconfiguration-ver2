import sys
import numpy as np
from run_method import RunMethod
from config import Configuration
from ybus import YbusMatrix
class NewtonRaphson(RunMethod):
    def __init__(self,config:Configuration):
        super().__init__(config=config)
    def __run1config__(self, fo=''):
        # ready to run Newton raphson
        Ybus = YbusMatrix(self.config).__getYbus__()
        #
        # Measure memory usage of the variable
        memory_usage = sys.getsizeof(Ybus)
        # Print the memory usage
        print(f"Memory usage of the variable: {memory_usage} bytes")
        res = {'FLAG':'CONVERGENCE','RateMax[%]':0, 'Umax[pu]':0,'Umin[pu]':100,'DeltaA':0,'cosP':0,'cosN':0}
        #
        if fo:
            rB,rG,rL = super().__initcsv__(fo)
        # return Y bus magnitude
        Ym = np.abs(Ybus)
        # return phase angle of Ybus        
        theta = np.angle(Ybus,deg=False)
        countSlack = 0
        countPV = 0
        slackCounted = []
        PVcounted = []
        for bi in self.param.BUS.keys():
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
            Vm = [float(self.param.Ubase) for _ in self.config.setBusHnd]
            P = [(-self.param.loadProfile[pi][bi]).real for bi in self.config.setBusHnd]    #+self.Pgen
            Q = [(-self.param.loadProfile[pi][bi]).imag for bi in self.config.setBusHnd]    #+self.Qgen
            # include distributed generation to Pbus
            for bi in self.config.dgOn:
                P[bi-1] += self.param.dgProfile[pi][bi].real
                Q[bi-1] += self.param.dgProfile[pi][bi].imag
            # update Vm of slack buses
            for bs in self.param.setSlack:
                Vm[bs-1] = self.param.genProfile[pi][bs]
            # initialize voltage angle 
            delta = [0 for _ in self.config.setBusHnd]
            sa1,dia1,va1 = dict(),dict(),dict()# for 1 profile
            for ii in range(self.param.iterMax+1):
                A = np.zeros((no_jacobi_equation,no_jacobi_equation))
                for b1 in self.config.setBusHnd:
                    J1_row_offdiag_idx = J2_row_offdiag_idx = int((b1 - slackCounted[b1-1])-1)
                    J1_diag_idx = J2_row_diag_idx = J3_col_diag_idx = J1_row_offdiag_idx 

                    J3_row_offdiag_idx = J4_row_offdiag_idx = int((self.param.nBus+b1-slackCounted[b1-1]-PVcounted[b1-1]-self.param.nSlack)-1)
                    J4_diag_idx = J2_col_diag_idx = J3_row_diag_idx = J3_row_offdiag_idx 

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
                                # diagonal elements of J1
                                J11 += Vm[b1-1] * Vm[b2-1] * Ym[b1-1][b2-1] * math.sin((theta[b1-1][b2-1] - delta[b1-1] + delta[b2-1]).real)
                                # diagonal elements of J3          
                                J33 += Vm[b1-1] * Vm[b2-1] * Ym[b1-1][b2-1] * math.cos((theta[b1-1][b2-1] - delta[b1-1] + delta[b2-1]).real)
                                if self.param.BUS[b1][3] != 3:
                                    #slackbus doesn't exist in J2 and J4
                                    #diagonal elements of J2
                                    J22 += Vm[b2-1] * Ym[b1-1][b2-1] * math.cos((theta[b1-1][b2-1] - delta[b1-1] + delta[b2-1]).real)
                                    #diagonal elements of J4               
                                    J44 += Vm[b2-1] * Ym[b1-1][b2-1] * math.sin((theta[b1-1][b2-1] - delta[b1-1] + delta[b2-1]).real)

                                if self.param.BUS[b1][3] != 3 and self.param.BUS[b2][3] != 3:
                                    J1_col_offdiag_idx = J3_col_offdiag_idx = int(b2 - slackCounted[b2-1] - 1)
                                    J2_col_offdiag_idx = J4_col_offdiag_idx = int(self.param.nBus + b2 - PVcounted[b2-1] - slackCounted[b2-1] - self.param.nSlack - 1)
                                    # off diagonal elements of J1
                                    A[J1_row_offdiag_idx][J1_col_offdiag_idx] = -Vm[b1-1] * Vm[b2-1] * Ym[b1-1][b2-1] * math.sin((theta[b1-1][b2-1] - delta[b1-1] + delta[b2-1]).real)
                                
                                    if self.param.BUS[b2][3] == 1:
                                        # off diagonal elements of J2
                                        A[J2_row_offdiag_idx][J2_col_offdiag_idx] = Vm[b1-1] * Ym[b1-1][b2-1] * math.cos((theta[b1-1][b2-1] - delta[b1-1] + delta[b2-1]).real)
                                
                                    if self.param.BUS[b1][3] == 1:
                                        # off diagonal elements of J3
                                        A[J3_row_offdiag_idx][J3_col_offdiag_idx] = -Vm[b1-1] * Vm[b2-1] * Ym[b1-1][b2-1] * math.cos((theta[b1-1][b2-1] - delta[b1-1] + delta[b2-1]).real)
                                
                                    if self.param.BUS[b1][3] == 1 and self.param.BUS[b2][3] == 1:
                                        # off diagonal elements of J4
                                        A[J4_row_offdiag_idx][J4_col_offdiag_idx] = -Vm[b1-1] * Ym[b1-1][b2-1] * math.sin((theta[b1-1][b2-1] - delta[b1-1] + delta[b2-1]).real)
                    #   
                    Pk = Vm[b1-1]**2 * Ym[b1-1][b1-1] * math.cos((theta[b1-1][b1-1])) + J33
                    Qk = -Vm[b1-1]**2 * Ym[b1-1][b1-1] * math.sin((theta[b1-1][b1-1])) - J11
                    if self.param.BUS[b1][3] == 3:
                        # swing bus
                        P[b1-1] = Pk
                        Q[b1-1] = Qk
                    if self.param.BUS[b1][3] == 2:
                        Q[b1-1] = Qk 
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
                        #diagonal elements of J1
                        A[J1_diag_idx][J1_diag_idx] = J11        
                        DC[J1_diag_idx] = P[b1-1] - Pk
                    if self.param.BUS[b1][3] == 1:
                        #diagonal elements of J2
                        A[J2_row_diag_idx][J2_col_diag_idx] = 2 * Vm[b1-1] * Ym[b1-1][b1-1] * math.cos(theta[b1-1][b1-1]) + J22    
                        #diagonal elements of J3
                        A[J3_row_diag_idx][J3_col_diag_idx] = J33
                        #diagonal elements of J4         
                        A[J4_diag_idx][J4_diag_idx] = -2 * Vm[b1-1] * Ym[b1-1][b1-1] * math.sin(theta[b1-1][b1-1]) - J44   
                        DC[J4_diag_idx] = Q[b1-1] - Qk
                #matrix A left division Matrix DC
                DX = np.linalg.solve(A, DC.T)
                for bi in self.config.setBusHnd:
                    del_update_DXidx = int(bi - slackCounted[bi-1]-1)
                    Vm_update_DXidx = int(self.param.nBus + bi - PVcounted[bi-1] - slackCounted[bi-1] - self.param.nSlack-1)
                    if self.param.BUS[bi][3] != 3:
                        delta[bi-1] +=  DX[del_update_DXidx]
                    if self.param.BUS[bi][3] == 1:
                        Vm[bi-1] += DX[Vm_update_DXidx]
                maxerror = max(abs(DC))
                if maxerror < self.param.epsilon:
                    break
                if ii==self.param.iterMax:
                    return {'FLAG':'DIVERGENCE'}
            #finish NR
            # store all voltage magnitude value in all profile
            va.extend(Vm)
            V = [Vm[bi-1]* np.cos(delta[bi-1])+Vm[bi-1]*1j*np.sin(delta[bi-1]) for bi in self.config.busC ]
            Il = dict()
            slt = 0
            vbus,sbus = dict(),dict()
            for i in range(len(Vm)):
                vbus[i+1] = V[i]
                sbus[i+1] = P[i] + Q[i]*1j
            for bi in self.param.busSlack:
                sbus[bi] += self.param.loadProfile[pi][bi]
            for li,bi in self.config.lineC.items():
                Ib1 = (vbus[bi[0]]-vbus[bi[1]])*(-Ybus[bi[0]-1][bi[1]-1])
                Ib2 = (vbus[bi[1]]-vbus[bi[0]])*(-Ybus[bi[0]-1][bi[1]-1])
                #
                if abs(Ib1) >= abs(Ib2):
                    Il[li] = abs(Ib1)
                else: 
                    Il[li] = abs(Ib2)   
                rate = ((Il[li])/self.param.LINE[li][3])*RATEC
                ra.append(rate)
                Snk = V[bi[0]-1]*np.conj(Ib1)
                Skn = V[bi[1]-1]*np.conj(Ib2)
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
