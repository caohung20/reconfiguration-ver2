import os
import csv
import time
from config import Configuration
def toString(v,nRound=5):
    """ convert object/value to String """
    if v is None:
        return 'None'
    t = type(v)
    if t==str:
        if "'" in v:
            return ('"'+v+'"').replace('\n',' ')
        return ("'"+v+"'").replace('\n',' ')
    if t==int:
        return str(v)
    if t==float:
        if v>1.0:
            s1 = str(round(v,nRound))
            return s1[:-2] if s1.endswith('.0') else s1
        elif abs(v)<1e-8:
            return '0'
        s1 ='%.'+str(nRound)+'g'
        return s1 % v
    if t==complex:
        if v.imag>=0:
            return '('+ toString(v.real,nRound)+' +' + toString(v.imag,nRound)+'j)'
        return '('+ toString(v.real,nRound) +' '+ toString(v.imag,nRound)+'j)'
    try:
        return v.toString()
    except:
        pass
    if t in {list,tuple,set}:
        s1=''
        for v1 in v:
            s1+=toString(v1,nRound)+','
        if v:
            s1 = s1[:-1]
        if t==list:
            return '['+s1+']'
        elif t==tuple:
            return '('+s1+')'
        else:
            return '{'+s1+'}'
    if t==dict:
        s1=''
        for k1,v1 in v.items():
            s1+=toString(k1)+':'
            s1+=toString(v1,nRound)+','
        if s1:
            s1 = s1[:-1]
        return '{'+s1+'}'
    return str(v)

def add2CSV(nameFile,ares,delim):
    """
    append array String to a file CSV
    """
    pathdir = os.path.split(os.path.abspath(nameFile))[0]
    try:
        os.mkdir(pathdir)
    except:
        pass
    #
    if not os.path.isfile(nameFile):
        with open(nameFile, mode='w') as f:
            ew = csv.writer(f, delimiter=delim, quotechar='"',lineterminator="\n")
            for a1 in ares:
                ew.writerow(a1)
            f.close()
    else:
        with open(nameFile, mode='a') as f:
            ew = csv.writer(f, delimiter=delim, quotechar='"',lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
            for a1 in ares:
                ew.writerow(a1)

class RunMethod:
    def __init__(self,config:Configuration):
        self.config = config
        self.param = config.param
        self.lineOff = self.config.lineOff
        self.shuntOff = self.config.shuntOff
    def __initcsv__(self,fo):
        add2CSV(fo,[[],[time.ctime()],['PF 1Profile','lineOff',str(list(self.lineOff)),'shuntOff',str(list(self.shuntOff))]],',')
        #
        rB = [[],['BUS/Profile']]
        rB[1].extend([bi for bi in self.param.lstBusHnd])
        #
        rL = [[],['LINE/Profile']]
        rL[1].extend([bi for bi in self.param.lstLineHnd])
        #
        rG = [[],['GEN/Profile']]
        for bi in self.param.busSlack:
            rG[1].append(str(bi)+'_P')
            rG[1].append(str(bi)+'_Q')
            rG[1].append(str(bi)+'_cosPhi')
        return rB,rG,rL
    def __update1profile__(self,pi,rB,rL,rG,sa1,va1,dia1):
        #sa1,dia1,va1 = dict(),dict(),dict() # for 1 profile
        rb1 = [pi]
        rl1 = [pi]
        rg1 = [pi]
        for bi1 in self.param.lstBusHnd:
            rb1.append(toString(abs(va1[bi1])/self.param.Ubase))
        #
        for bri in self.param.lstLineHnd:
            try:
                r1 = abs(dia1[bri])/self.param.LINE[bri][3]*RATEC
                rl1.append( toString(r1,2) )
            except:
                rl1.append('0')
        #
        for bs1 in self.param.busSlack:
            rg1.append(toString(sa1[bs1].real))
            rg1.append(toString(sa1[bs1].imag))
            if sa1[bs1].imag>=0:
                rg1.append(toString(sa1[bs1].real/abs(sa1[bs1]),3))
            else:
                rg1.append(toString(-sa1[bs1].real/abs(sa1[bs1]),3))
        #
        rB.append(rb1)
        rL.append(rl1)
        rG.append(rg1)
        return rB,rL,rG
    
    def __update_result__(self,res,va,ra,cosP,cosN):
        res['Umax[pu]'] = max(va)/self.param.Ubase
        res['Umin[pu]'] = min(va)/self.param.Ubase
        res['RateMax[%]'] = max(ra)
        res['cosP'] = min(cosP)
        res['cosN'] = max(cosN)
        return res

    def __export_profiles__(self,fo,rB,rL,rG,res):
        rB.append(['','Umax[pu]',toString(res['Umax[pu]']),'Umin[pu]',toString(res['Umin[pu]']) ])
        add2CSV(fo,rB,',')
        #
        rL.append(['','RateMax[%]',toString(res['RateMax[%]'],2)])
        add2CSV(fo,rL,',')
        #
        rG.append(['','cosPmin',toString(res['cosP'],3),'cosNMax',toString(res['cosN'],3)])
        add2CSV(fo,rG,',')
        return  