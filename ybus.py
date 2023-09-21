from config import Configuration
class YbusMatrix:
    
    def __init__(self, config:Configuration):
        self.config = config
        self.param = config.param
        self.setLineHnd = config.setLineHnd
        self.setLinebHnd = config.setLinebHnd
        self.shuntOff = config.shuntOff

    def __get_sparse_Ybus__(self):
        sparse_ybus = dict()
        #initialize branch admittance
        y = {k:(1+0j) for k in self.setLineHnd}
        #formation of the off diagonal elements
        for k,v in self.config.lineC.items():
            y[k] = y[k]/self.param.LINE[k][2]
            sparse_ybus[(v[0],v[1])] = -y[k]
            sparse_ybus[(v[1],v[0])] = sparse_ybus[(v[0],v[1])]
        #calculate yline with lineb: 
        for lbi in self.setLinebHnd:
            y[lbi] += self.param.LINEb[lbi]*1j
        # Shunt 
        for k1,v1 in self.param.BUSbs.items():
            if k1 not in self.shuntOff:
                sparse_ybus[k1] = v1*1j
        #formation of the diagonal elements
        for bi in self.config.busC.keys():
            for li in self.config.busC[bi]:
                if bi not in sparse_ybus:
                    sparse_ybus[bi] =  y[li]
                else:
                    sparse_ybus[bi] +=  y[li]
        # case if bus slack does not connect to any bus
        for bsi in self.param.busSlack:
            if not self.config.busC[bsi]:
                sparse_ybus[bsi] = 0
        return sparse_ybus
    #
    def __getYbus__(self):
        #initialize Ybus
        Ybus = []
        for _ in range(self.param.nBus):
            row = [0] * self.param.nBus
            Ybus.append(row)
        #initialize branch admittance
        y = {k:(1+0j) for k in self.config.setLineHnd}
        #formation of the off diagonal elements
        for k,v in self.config.lineC.items():
                y[k] = y[k]/self.param.LINE[k][2]
                Ybus[v[0]-1][v[1]-1] = -y[k]
                Ybus[v[1]-1][v[0]-1] = Ybus[v[0]-1][v[1]-1]
        #calculate yline with lineb: 
        for lbi in self.setLinebHnd:
            y[lbi] += self.param.LINEb[lbi]*1j
        #formation of the diagonal elements
        for bi in self.config.busC.keys():
            for li in self.config.busC[bi]:
                Ybus[bi-1][bi-1] +=  y[li]
        # Shunt 
        for k1,v1 in self.param.BUSbs.items():
            if k1 not in self.shuntOff:
                Ybus[k1-1][k1-1] += v1*1j
        return Ybus

    def __get_sparse_ybus_list__(self):
        row,col,data = [],[],[]
        #initialize branch admittance
        y = {k:(1+0j) for k in self.setLineHnd}
        #formation of the off diagonal elements
        for k,v in self.config.lineC.items():
            y[k] = y[k]/self.param.LINE[k][2]
            ib1 = v[0] - 1
            ib2 = v[1] - 1
            row.append(ib1)
            col.append(ib2)
            row.append(ib2)
            col.append(ib1)
            data.append(-y[k])
            data.append(-y[k])
        #calculate yline with lineb: 
        for lbi in self.setLinebHnd:
            y[lbi] += self.param.LINEb[lbi]*1j
        diag = [(0+0j) for _ in self.config.setBusHnd]
        for k1,v1 in self.param.BUSbs.items():
            if k1 not in self.shuntOff:
                diag[k1-1] += v1*1j
        #formation of the diagonal elements
        for bi in self.config.busC.keys():
            for li in self.config.busC[bi]:
                diag[bi-1] +=  y[li]
        for i in range(len(self.config.setBusHnd)):
            col.append(i)
            row.append(i)
        data.extend(diag)
        sparse_Ybus = csr_matrix((data, (row, col)), shape=(self.param.nBus,self.param.nBus))
        return sparse_Ybus

