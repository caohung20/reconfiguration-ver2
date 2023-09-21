import copy
from parameter import Parameter
class Configuration:
    def __init__(self,param: Parameter,lineOff=[],shuntOff=[],varFlag=None):
        self.param = param
        if varFlag is not None:
            if len(varFlag)!=self.param.nVar:
                raise Exception('Error size of varFlag')
            lineOff = self.getLineOff(varFlag[:self.param.nL])
            shuntOff = self.getShuntOff(varFlag[self.param.nL:])
        #
        self.setBusHnd =self.param.setBusHnd
        self.setLineHnd = self.param.setLineHndAll - set(lineOff)
        self.setLinebHnd = self.param.LINEb.keys() - set(lineOff)
        #
        self.shuntOff = set(shuntOff)
        self.lineOff = set(lineOff)
        #
        self.lineSureISL,busc1 = self.__getLineISL__(self.param.BUSC)
        self.bus0ISL = set(busc1.keys())
        self.busISL = self.setBusHnd - self.bus0ISL
        #
        self.lineC = copy.deepcopy(self.param.AllLine2Bus.copy())
        #
        self.busC = copy.deepcopy(self.param.BUSC)

        for li in self.lineOff:
            for bi in self.lineC[li]:
                self.busC[bi].remove(li)
        #
        self.bus2bus = copy.deepcopy(self.param.bus2bus)
        for li in self.lineOff:
            bi = self.lineC[li]
            self.bus2bus[bi[0]].remove(bi[1])
            self.bus2bus[bi[1]].remove(bi[0])
        #
        for li in self.lineOff:
            self.lineC.pop(li) 
        #
        #Dg On
        self.dgOn = self.param.dgOn

        
        
    #
    def getLineFlag3(self):
        """ cac Branch co the dong mo """
        return self.param.lineFLAG3
    #
    def getShuntFlag3(self):
        """ cac Shunt co the dong mo """
        return self.param.shuntFLAG3
    #
    def getLineOff(self,lineFlag): # 0: inservice, 1 off service
        lineOff = []
        for i in range(len(self.param.lineFLAG3)):
            if lineFlag[i]:
                lineOff.append(self.param.lineFLAG3[i])
        return lineOff
    #
    def getShuntOff(self,shuntFlag): # 0: inservice, 1 off service
        shuntOff = []
        for i in range(len(self.param.shuntFLAG3)):
            if shuntFlag[i]:
                shuntOff.append(self.param.shuntFLAG3[i])
        return shuntOff
    #
 
    def getVarFlag(self,lineOff,shuntOff,dgOff):
        varFlag = [0]*self.param.nVar
        for i in range(self.param.nL):
            if self.param.lineFLAG3[i] in lineOff:
                varFlag[i] = 1
        for i in range(self.param.nSht):
            if self.param.shuntFLAG3[i] in shuntOff:
                varFlag[self.param.nL+i] = 1
        return varFlag
    def __checkLoopIsland__(self,lineOff):
        # check island/loop multi slack ----------------------------------------
        if lineOff.intersection(self.lineSureISL):
            return 'ISLAND'
        #
        self.setLineHnd = self.param.setLineHndAll - lineOff
        
        r11 = self.setBusHnd.copy()
        self.busGroup = []# cac bus tuong ung o cac slack khac nhau
        for bs1 in self.param.busSlack:
            r1 = self.__findBusConnected__(bs1,self.busC,self.lineC)
            if len(r1.intersection(self.param.setSlack))>1:
                return 'LOOP MULTI SLACK'
            #
            self.busGroup.append(r1)
            r11.difference_update(r1)
        #
        if r11:
            return 'ISLAND'
        #
        # LOOP
        if self.__checkLoop__(self.bus0ISL,self.busC,self.lineC):
            return 'LOOP'
        #
        return ''
        
    def __findBusConnected__(self,bi1,bset,lset):
        ## find all bus connected to Bus b1
        res = {bi1}
        ba = {bi1}
        while True:
            ba2 = set()
            for b1 in ba:
                for li in bset[b1]:
                    for bi in lset[li]:
                        if bi not in res:
                            ba2.add(bi)
            if ba2:
                res.update(ba2)
                ba=ba2
            else:
                break
        return res
    def __checkLoop__(self,busHnd,bus,br):
        setBusChecked = set()
        setBrChecked = set()
        for o1 in busHnd:
            if o1 not in setBusChecked:
                setBusChecked.add(o1)
                tb1 = {o1}
                #
                for i in range(20000):
                    if i==19999:
                        raise Exception('Error in checkLoop()')
                    tb2 = set()
                    for b1 in tb1:
                        for l1 in bus[b1]:
                            if l1 not in setBrChecked:
                                setBrChecked.add(l1)
                                for bi in br[l1]:
                                    if bi!=b1:
                                        if bi in setBusChecked or bi in tb2:
                                            return bi
                                        tb2.add(bi)
                    if len(tb2)==0:
                        break #ok finish no loop for this group
                    setBusChecked.update(tb2)
                    tb1=tb2.copy()
        return None
    def __getLineISL__(self,busC0):
        lineISL = set() # line can not be off => ISLAND
        busC = busC0.copy()
        while True:
            n1 = len(lineISL)
            for k,v in busC.items():
                if len(v)==1:
                    lineISL.update(v)
            if n1==len(lineISL):
                break
            busc1 = dict()
            for k,v in busC.items():
                if len(v)!=1:
                    busc1[k]=v-lineISL
            busC = busc1.copy()
        return lineISL,busC
    def dfs(self, node, parent, path):
        if node in path:
            start_idx = path.index(node)
            loop_found = path[start_idx:]
            set_loop_found = set (loop_found)
            if set_loop_found not in self.setloops.values():
                l = len(self.loops)
                self.loops[l+1] = loop_found
                self.setloops[l+1] = set_loop_found
            return 
        path.append(node)
        for neighbor in self.bus2bus[node]:
            if neighbor != parent:
                self.dfs(neighbor, node, path.copy())
        path.pop()
    def getloop(self,node):
        self.loops = dict()
        self.setloops = dict()
        
        self.dfs(node, None, [])
        return self.loops
    def get_lines_from_loop(self,lineOn):
        node = self.lineC[lineOn][0]
        loop = self.getloop(node)
        self.line_in_loop = set()
        for i in range(len(loop[1])):
            line = {loop[1][i-1],loop[1][i]}
            position = self.param.busesin1line.index(line)
            self.line_in_loop.add(self.param.listline[position])
        return self.line_in_loop
    
    def get_list_lines_in_loop(self,lineOn):
        node =self.lineC[lineOn][0]
        loop = self.getloop(node)
        self.line_in_loop = []
        for i in range(len(loop[1])):
            line = {loop[1][i-1],loop[1][i]}
            position = self.param.busesin1line.index(line)
            self.line_in_loop.append(self.param.listline[position])
        
        return self.line_in_loop


if __name__ == '__main__':
    fi = 'tromvia200percentSmax.xlsx'
    param = Parameter(fi)
    lineoff = [33,35,36,37]
    config = Configuration(param,lineoff)
    lineinloop = config.get_lines_from_loop(2)
    print(lineinloop)