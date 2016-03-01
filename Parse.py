#!/usr/bin/python

import xml.etree.ElementTree as ET
import sys,os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import SolveEq
import collections
from scipy.stats.stats import pearsonr
#TODO: The scirp is merely functional. Need be cleaned For Lixing Fri 29 Jan 2016 04:08:17 PM EST.

DEBUG  = int(sys.argv[1])
PLOT   = DEBUG
FORMAL = not(DEBUG)
CELLULAR = False
GAP = False


def countAggregation(interrtt,window,thred):
    pct=[]
    pdt=[]
    start=0
    for i in range(0,len(interrtt)-window):
        pcount = 0
        sum = 0
        psum = 0
        for j in interrtt[i:i+window]:
            if j > thred:
                pcount += 1
                psum += j
            sum += abs(j)
        pct.append(pcount/float(window))
        pdt.append(psum/float(sum))

    #for i in range(window):
    #    pct.append(pcount/float(window))

    return pct,pdt


def statInterrtt(list):
    pSum=[]
    nSum=[]
    for i,x in enumerate(list):
        if x > 0:
            pSum.append(x)
        elif x< 0:
            nSum.append(x)

    if len(nSum) == 0 and len(pSum) != 0:
        pnratio = float('inf')
        pct = len(pSum)/float(len(pSum)+len(nSum))
    elif len(pSum) == 0 and len(nSum) == 0:
        pnratio = 1
        pct = 1
    else:
        pnratio=sum(pSum)/float(-sum(nSum))
        pct = len(pSum)/float(len(pSum)+len(nSum))

    return pnratio,pct



def getProbeRate(intersent,packetSize):
    return [(int(packetSize))*8/float(sent*1000) for sent in intersent]



def getEstResult_Aggre(interrtt,interack,intersent,packetSizes,N,Min,Max,Step,Up_Down):
    def getAckRate(interack,start,end):
        real_start = start
        real_end = end
        thred = 0.4#np.std(interack)
        for ii in range(0,len(interack),1):
            # print "start",interack[start + ii],interack[start - ii]
            if start > -1 and start<len(interack):
                if start + ii < len(interack) and (interack[start + ii] > thred) :
                    real_start = ii + start
                    break

                if start - ii > 1 and (interack[start - ii] > thred) :
                    real_start = start - ii
                    break
        end = min(len(interack)-1,end)
        for ii in range(0,len(interack),1):
            if end > -1 and end<len(interack):
                if end + ii > real_start and end + ii < len(interack) and (interack[end + ii] > thred) :
                    real_end = ii + end-1
                    break
                if end - ii > real_start and (interack[end - ii] > thred) :
                    real_end = end - ii-1
                    break

        # for ii in range(start,0,-1):
        #     if ii < len(interack) and ii>0 and (interack[ii] > thred) :
        #         real_start = ii
        #         break
        # for ii in range(end,len(interack),1):
        #     if ii < len(interack) and ii>0 and interack[ii] > thred:
        #         real_end = ii - 1
        #         break
        # print start,end,":",real_start,real_end
        return real_start,real_end

    pnratio_st,pct_all = statInterrtt(interrtt[:N-1])
    pct_all = sum([ inter<0.4 for inter in interack[:N-1]])/float(len(interack[:N-1]))
    pcount_thre = 0
    pcount_len = 10
    step=int(Step)
    Lsubtrain=int(N)/step
    proberate=[float(Min) + (float(Max)-float(Min))/float(N)*(i/Lsubtrain)*Lsubtrain for i in range(int(N))]
    if FORMAL:
        print "%2.3f"%pnratio_st,"%2.3f"%pct_all,
    elif DEBUG:
        print " PN %2.3f"%pnratio_st,"PCT %2.3f"%pct_all,
    start = int(0)
    bundle_thre = 0.2
    aggre_thre = 0.4#np.std(interack)
    pct,pdt=countAggregation(interack[:N],5,aggre_thre)

    pcount = 0
    p_b = 0
    b_g = 0
    unbundled={}
    for ix,x in enumerate(pct):
        if ( x > bundle_thre  ) :#and pct[ix-1] > bundle_thre
            b_g = ix - p_b
            # print "Unaggre:",ix,x,b_g
            if b_g > 5:
                start = p_b - pcount+1
                # print "Breaks",start,pcount
                # pcount += 1
                while pcount > pcount_len:
                    unbundled[start] = pcount_len
                    pcount -= pcount_len
                    start += pcount_len
                if pcount > 3:
                    unbundled[start] = pcount
                pcount = 1
            else:
                if pcount == 0:
                    start = ix
                pcount += max(int( ix - p_b ),1)
            p_b = ix
        if ix == len(pct) - 1 and pcount > 3:
            start =  p_b - pcount+1
            # pcount += 1
            while pcount > pcount_len:
                unbundled[start] = pcount_len
                pcount -= pcount_len
                start += pcount_len
            if pcount > 3:
                unbundled[start] = pcount
            pcount = 0



    start = 0
    start_b=start
    start_e=start
    unbundled = collections.OrderedDict(sorted(unbundled.items()))
    if DEBUG:
        print "Estimation Aggre -----------------------------------"
        print "Get starts:",unbundled
    pre_true=True
    c_gap = 0
    last_b=-1
    last_t = -1
    accu_false = 0
    first_false = True
    dict_accu_false = {}
    b_gap = {}
    dict_accu_false[0]=0
    b_gap[0] = 0
    rates= {}
    rates[0] = 1
    for p_start in unbundled:
        count=unbundled[p_start]
        # pp_start = int( p_start/10 )*10

        #rate = (sum(packetSizes[pp_start:pp_start+10])*8)/float(1000*sum(interack[pp_start:pp_start+10]))
        #s_rate = (sum(packetSizes[pp_start:pp_start+10])*8)/float(1000*sum(intersent[pp_start:pp_start+10]))
        r_start = 0
        if p_start == 0:
            r_start = 0
        else:
            r_start = p_start - 1

        # real_start,real_end = getAckRate(interack,r_start,r_start + count - 1)
        real_start,real_end = getAckRate(interack,p_start,p_start + count - 1)

        try:
            rate = (sum(packetSizes[real_start:real_end+1])*8)/float(1000*sum(interack[real_start:real_end+1]))
            s_rate = (sum(packetSizes[real_start:real_end+1])*8)/float(1000*sum(intersent[real_start:real_end+1]))
        except:
            rate = 1
            s_rate = 1
        #if s_rate > proberate[p_start + count/2]:
        #    s_rate = proberate[p_start + count/2]
        delta = min(0.1*s_rate,0.5)
        if rate  < s_rate + 2 and rate+delta > s_rate and rate<float(Max) :
            if pre_true:
                start_b=p_start + count

                #first_false = True
            #if pre_true:
            #    start = p_start + count
            #    start_e = start
            #    rates[start] = rate
            #else:
            #    pre_true = True

            start = p_start + count
            start_e = start
            rates[start] = rate
            pre_true = True

            last_t = p_start
            dict_accu_false[p_start+count] = accu_false
        else:
            #if pre_true == False:
            #    break
            pre_true = False
            accu_false += count
            if accu_false > 19:
                break
            #if first_false == True:
            #    start_b = p_start + count
            #    rates[start_b] = rate
            #    first_false = False
            dict_accu_false[p_start+count] = accu_false

        if last_b == -1:
            c_gap = p_start
        else:
            c_gap = p_start - (last_b + unbundled[last_b])
        #if c_gap > 19:
        #    break
        #if c_gap > 10:
        #    break
        if  last_b > -1 and b_gap[(last_b + unbundled[last_b])] < c_gap:
            b_gap[p_start + count] = c_gap
        elif last_b == -1:
            b_gap[p_start + count] = c_gap
        else:
            b_gap[p_start + count] = b_gap[last_b + unbundled[last_b]]
        last_b = p_start
        if DEBUG:
            print "Verify start: ",p_start, rate, s_rate, \
                rate  < s_rate + 2 and rate+delta > s_rate and rate<float(Max),\
                "Last_True",last_b,"Gap",c_gap,b_gap[p_start + count],\
                "False:",accu_false

    #if start == 0 and pnratio_st < 1.2 :
    #    if last_t != -1:
    #        start = last_t + unbundled[last_t]
    #        print "guess",
    #    else:
    #        print "bad",

    #final evaluate
    real_start,real_end = getAckRate(interack,start,N-1)
    real_end = N - 2
    rate = (sum(packetSizes[real_start:real_end+1])*8)/float(1000*sum(interack[real_start:real_end+1]))
    s_rate = (sum(packetSizes[real_start:real_end+1])*8)/float(1000*sum(intersent[real_start:real_end+1]))
    if DEBUG:
        print "Final Verify start: ",start,rate,s_rate,rate+0.1  > s_rate
        print "Conf.",int(pct_all>0.2 or Up_Down>0.5 ),dict_accu_false[start],b_gap[start]
    elif FORMAL:
        print "Conf.",int(pct_all>0.2 or Up_Down>0.5 ),dict_accu_false[start],b_gap[start],
    if rate + 0.1  > s_rate and pct_all>0.17 and start > 0:
        dict_accu_false[int( 0.9*int(N) )]=dict_accu_false[start]
        b_gap[int( 0.9*int(N) )]=b_gap[start]
        start = int( 0.9*int(N) )
        rates[start] = rate
    pnratio,pct_ = statInterrtt(interrtt[start:-2])
    #while start < len(pct)-1 and (proberate[start]>40):#pnratio<1.15 or
    #    pnratio,pct_ = statInterrtt(interrtt[start:-2])
    #    if DEBUG:
    #        print "END Verify start:",start,pnratio,pct_
    #    start += 1
    #################################################
    ## Check the achivable Thoughtput ###############
    #################################################
    #real_start,real_end = getAckRate(interack,90,99)
    #try:
    #    rate = (sum(packetSizes[real_start:real_end+1])*8)/float(1000*sum(interack[real_start:real_end+1]))
    #except:
    #    rate = 30
    #dd = 0
    #if DEBUG:
    #    print "Cal AT:",real_start,rate
    #while dd < 90 and ( rate > 12.5 or real_start>90 ):
    #    dd += 10
    #    real_start,real_end = getAckRate(interack,90-dd,99)
    #    try:
    #        rate = (sum(packetSizes[real_start:real_end+1])*8)/float(1000*sum(interack[real_start:real_end+1]))
    #    except:
    #        rate = 30
    #    if DEBUG:
    #        print "Cal AT:",real_start,rate
    #s_rate = (sum(packetSizes[real_start:real_end])*8)/float(1000*sum(intersent[real_start:real_end]))
    #rate = (sum(packetSizes[real_start:real_end])*8)/float(1000*sum(interack[real_start:real_end]))
    ##print "AT: %2.2f/%2.2f"%(rate,s_rate),
    aggre_thre = 0.4#np.std(interack)
    if FORMAL or DEBUG:
        ei = []
        for ij,i_interack in enumerate(interack):
            if i_interack > aggre_thre and ij > N -1 :
                ei.append(i_interack)
        print "EI: %2.2f"%(np.mean(ei)),
    #print "C: %2.2f"%(rate*(s_rate-rates[start_e])/float(s_rate - rate)),



    result=[]
    if start == int( 0.9*int(N) ) or ( rates[start_b] + rates[start_e] )/2.0 >12:
        result.append('Nan')
        result.append('Nan')
        #return (proberate[start],None)
    else:
        result.append(  rates[start_b] )
        result.append(  rates[start_e] )
        #return (None, proberate[start])
    return result




def detectUp(list,list1,N):
    # base=None
    # base1=None
    if N >= len(list):
        N = len(list) + 1
    # for i in range(N-1):
    #     if list[i] > list[0] and list1[i] > list1[0]:
    #         base = list[i]-list[0]
    #         base1 = list1[i]-list1[0]
    #         break

    final=[]
    final1=[]
    # if base is None:
    #     base = 1
    # if base1 is None:
    #     base1 = 1
    for i in range(N-2):
        t_deviation = ( list[i+1]-list[i] )/float(list[N-2]-list[0])
        final.append(t_deviation)

        t_deviation = ( list1[i+1]-list1[i] )/float(list1[N-2]-list1[0])
        # t_deviation = ( list1[i+1]-list1[i] )/float(base1)

        final1.append(t_deviation)

    return final,final1

def getTime(str):
    temp=str.split('.')
    return float(temp[0]) + float(temp[1])/1000000
def loadXML_Parse(filename):
    df = pd.DataFrame({ 'LinkType' : [],
                       'Power' :[] ,
                       'Time_PktSent' :[] ,
                       'Time_PktRcvd' :[] ,
                       'Time_AckSent' :[] ,
                       'Size_PktSent' :[] })

    tree = ET.parse(filename)
    root = tree.getroot()
    request=tree.getroot().get('Request')
    #if request is None or "Martins" not in request:
    #    return None,None,None,None
    params=request.split('?')
    if 'train' in params[0]:
        global CELLULAR,GAP
        CELLULAR = int(params[1])
        if len(params) > 8:
            if params[2].split('=')[0] == 'PacketGap':
                GAP = True
            packetSize=params[2].split('=')[1]
            packetGap=params[3].split('=')[1]
            packetN=params[4].split('=')[1]
            Rmin=params[5].split('=')[1]
            Rmax=params[6].split('=')[1]
            #In Lab test there is one more params Step
            if "Step" in params[7]:
                Step=params[7].split('=')[1]
                Capacity=params[8].split('=')[1]
                Util=params[9].split('=')[1]
            else:
                Step=5
                Capacity=params[7].split('=')[1]
                Capacity=Capacity.replace("Mbps",'')
                if len(params)>10 and params[9] == 'Android':
                    Util=params[11].split('=')[1]
                else:
                    Util=params[8].split('=')[1]

        else:
            if params[2].split('=')[0] == 'PacketGap':
                GAP = True
            packetSize=1400
            packetGap=1200
            packetN=100
            Rmin=1
            Rmax=15
            Capacity="NA"
            Util="NA"
    elif 'custom' in params[0]:
        packetSize=params[6].split('=')[1]
        packetGap=params[5].split('=')[1]
        packetN=params[1].split('=')[1]
        Rmin="00"
        Rmax="00"
        Capacity="NA"
        Util="NA"



    if request is None:
        return None,None,None,None

    linktype=tree.getroot().get('Request')[11:]
    power=request[18:]
    recvlist=tree.getroot().find('MeasureRcvd')
    # recvlist=tree.getroot().findall(".//MeasureRcvd/*[@Meta]")
    #print len( testlist )
    sentlist=tree.getroot().find('MeasureSent')
    idx=0
    if  len(sentlist) - len(recvlist) <1 and\
            len(sentlist) - len(recvlist) >-3:
        for i,item in enumerate(sentlist):
            df=df.append(pd.DataFrame({'LinkType': [linktype] ,
                                       'Power':[ power ],
                                       'Time_PktSent':[ getTime( item.get('Time') ) ],
                                       'Time_PktRcvd':[ getTime( recvlist[i].get('Time') ) ],
                                       'Time_AckSent':[  int( recvlist[i].get('TsVal'))   ],
                                       'Size_PktSent':[ int( item.get('IPLength') ) ],
                                       },index=[idx]))
            idx += 1

    elif  len(sentlist) - len(recvlist) >0 and \
            len(sentlist) - len(recvlist) < 50:
        for i,item in enumerate(recvlist):
            df=df.append(pd.DataFrame({'LinkType': [linktype] ,
                                       'Power':[ power ],
                                       'Time_PktSent':[ getTime( sentlist[i].get('Time') ) ],
                                       'Time_PktRcvd':[ getTime( item.get('Time') ) ],
                                       'Time_AckSent':[  int( item.get('TsVal'))   ],
                                       'Size_PktSent':[ int( sentlist[i].get('IPLength') ) ],
                                       },index=[idx]))
            idx += 1

    if len(df) == 0:
        if DEBUG:
            print len(recvlist),len(sentlist)
            print "Error, len(df) == 0"
        return None,None,None,None
    df['RTT']=1000*(df['Time_PktRcvd']-df['Time_PktSent'])
    segments=[len(df)-2]#[14,45,69]#,20,20,20,20,19]#[14,15,15,15,15,14]
    start=0
    interack_all=[]
    interrtt_all=[]
    intersent_all=[]
    realRate=[]
    seg_start=[0,2,2]
    if FORMAL or DEBUG:
        print filename,(filename.split('-')[-1]).replace('.xml',''),
        #1000*( df.at[len(df)-1,'Time_PktRcvd']-df.at[0,'Time_PktSent'] )
    if FORMAL:
        print int( CELLULAR ),packetN,packetSize,packetGap,Rmax,Capacity,Util,\
            "%2.2f"%(len(df)/float(int(packetN) + 40)),\
            "%2.2f"%(np.mean(df.ix[0:int(packetN),['RTT']])),
    if DEBUG:
        print "********Bandwidth*******************"
    for iseg,seg in enumerate(segments):
        end = start + seg - 1

        interack=[]
        intersent=[]
        interrtt=[]
        up_interrtt=[]
        down_interrtt=[]


        for i in range(start+seg_start[iseg],end):
            interack.append(1000* (df.ix[i+1,'Time_PktRcvd']-df.ix[i,'Time_PktRcvd']) )
            intersent.append(1000* (df.ix[i+1,'Time_PktSent']-df.ix[i,'Time_PktSent']) )
            interrtt.append(1* (df.ix[i+1,'RTT']-df.ix[i,'RTT']) )


        interrtt_all += interrtt
        interack_all += interack
        intersent_all += intersent

        interack_mean=np.mean(interack)

        intersent_mean=np.mean(intersent)

        df=df[:-1]

        start = start + seg


        up1,up2 = detectUp(df['Time_AckSent'].tolist(),df['Time_PktRcvd'].tolist(),int( packetN ))
        up_down=(pearsonr(up1,up2)[0])
        estimation_Aggre=getEstResult_Aggre(interrtt,interack,intersent,df['Size_PktSent'].tolist(),int( packetN ),Rmin,Rmax,Step,up_down)
        if FORMAL or DEBUG:
            print "Est:", \
                "Aggre_ %2.2f"%float(estimation_Aggre[0]),\
                "%2.2f"%float(estimation_Aggre[1]),Step,





    return df,interack_all,interrtt_all,intersent_all
def main (argv):
    if PLOT:
        f,(plt1,plt2,plt3,plt5)=plt.subplots(4,sharex=True)
    for filename in argv:
        df,interack,interrtt,intersent = loadXML_Parse(filename) # the filename
        if df is not None and PLOT:
            plt1.plot(df['RTT'], linewidth=2.0,label="RTT"+filename)
            plt2.plot(interack,linewidth=2.0,label="InterAck"+filename)
            pct,pdt=countAggregation(interack[:100],5,0.4)


            plt3.plot(pct,marker='x',label="PCT"+filename)
            up1,up2 = detectUp(df['Time_AckSent'].tolist(),df['Time_PktRcvd'].tolist(),100)
            plt5.plot(up1, '-o',linewidth=2.0,label="ACK_SENT"+filename)
            plt5.plot(up2, '-+',linewidth=2.0,label="Pkt_Rcvd"+filename)
        if df is not None:
            up1,up2 = detectUp(df['Time_AckSent'].tolist(),df['Time_PktRcvd'].tolist(),100)
            if FORMAL or DEBUG:
                print "Cor. %2.3f"%(pearsonr(up1,up2)[0])

    if PLOT:
        plt1.set_xlabel("Packet Unit")
        plt2.set_xlabel("Packet Unit")
        plt3.set_xlabel("Packet Unit")
        plt5.set_xlabel("Packet Unit")
        plt1.set_ylabel("RTT (ms)")
        plt2.set_ylabel("Inter-Arrival of ACKs (ms)")
        plt3.set_ylabel("Local PCT")
        plt5.set_ylabel("InterAck vs InterSent")
        plt1.legend(loc='upper left')
        plt2.legend(loc='upper left')
        plt3.legend(loc='lower right')
        plt5.legend(loc='lower right')
        plt.show()

if __name__ == "__main__":
    main(sys.argv[2:])
