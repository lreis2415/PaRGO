# coding:utf-8

import geopandas as gpd
import queue
import math
import time
import os
from streamRasToShp import streamSHP
import hydroAnalysis as hydro
cur_dir = os.getcwd()

workdir = "cd wkflow & "

nrank = 8  # 进程数
mpi = r"mpiexec -n " + str(nrank)
#文件夹与中间结果名称定义
middata = r"data\\mid\\"
finaldata = r"data\\final\\"
nbr = r"data\\neighbor\\moore.nbr "

demres = middata + r"res_dem.tif "
demresfel = middata + r"res_demfel.tif "
resdirfile = middata + r"res_dir.tif "
resscafile = middata + r"res_sca.tif "
weifile = middata + r"weight.tif "
ssafile = middata + r"w_sca.tif "
resnet = middata + r"res_net.tif "
res_so = middata + r"res_str.tif "
res_ws = middata + r"res_watershed.tif "
netshp = middata + r"res_rivernet.shp"
resd8file = middata + r"res_d8.tif "
demfel = finaldata + r"demfel.tif "
dirfile = finaldata + r"dir.tif "
scafile = finaldata + r"sca.tif "
netfile = finaldata + r"net.tif "
streamfile = finaldata + r"streamOrder.tif "
wsfile = finaldata + r"watershed.tif "








def largearea(_g, buf, watershed, dir, sca):
    exe = r" apps\Largearead8.exe "
    cmd = workdir + mpi + exe + res_ws + dirfile + nbr + scafile + str(_g) + " " + str(buf) + " "
    net = gpd.read_file(netshp)
    l = len(net['Watershed'])
    m = net['Watershed']
    '''
    while (q.qsize() < l):
        for it in net.iloc:
            if (it.USLINKNO1 == -1 and it.USLINKNO2 == -1 and it.LINKNO != -1):
                q.put(it.LINKNO)
                index1 = net.loc[net.USLINKNO1 == it.LINKNO].index.to_list()
                index2 = net.loc[net.USLINKNO2 == it.LINKNO].index.to_list()
                iin = net.loc[net.LINKNO == it.LINKNO].index.to_list()
                net.loc[iin[0], 'LINKNO'] = -1
                if (len(index1)):
                    net.loc[index1[0], 'USLINKNO1'] = -1
                # print('up1:',net.loc[index1[0],'USLINKNO1'])
                if (len(index2)):
                    net.loc[index2[0], 'USLINKNO2'] = -1
'''
    c = cmd + str(1)  # init
    g = os.system(c)
    t = 0
    start = time.time()
    for i in range(0, l):
        id = str(m[i])
        c = cmd + id
        print("working process: large_aread8-", id)
        g = os.system(c)
        #print(g)
        t = t + g
    end = time.time()
    print("process time:", end - start)
    return t / 10000

def largePitremove(_g, buf, watershed, dem, demfel):
    exe = r" apps\Largepitremove.exe "

    cmd = workdir + mpi + exe + watershed + dem + nbr + demfel + str(_g) + " " + str(buf) + " "
    print(cmd)
    net = gpd.read_file(netshp)
    l = len(net['Watershed'])
    m = net['Watershed']

    '''
    while (q.qsize() < l):
        for it in net.iloc:
            if (it.USLINKNO1 == -1 and it.USLINKNO2 == -1 and it.LINKNO != -1):
                q.put(it.LINKNO)
                index1 = net.loc[net.USLINKNO1 == it.LINKNO].index.to_list()
                index2 = net.loc[net.USLINKNO2 == it.LINKNO].index.to_list()
                iin = net.loc[net.LINKNO == it.LINKNO].index.to_list()
                net.loc[iin[0], 'LINKNO'] = -1
                if (len(index1)):
                    net.loc[index1[0], 'USLINKNO1'] = -1
                # print('up1:',net.loc[index1[0],'USLINKNO1'])
                if (len(index2)):
                    net.loc[index2[0], 'USLINKNO2'] = -1
    '''
    print("number of watersheds:", l)
    c = cmd + str(1)
    g = os.system(c)
    start = time.time()
    for i in range(0, l):
        id = str(m[i])
        print(id)
        c = cmd + id
        print("working process: large_pitremove-", id)
        g = os.system(c)
        print(g)
    end = time.time()
    print("process time:", end - start)



def hydroAnalysis(_g, buf, nrank,  threshold = 5000):
    hydro.nrank=nrank
    hydro.resample(dem, demres, _g)
    hydro.pitremove(demres, demresfel)

    #hydro.flowdirdinf(demresfel, resdirfile)
    hydro.flowdird8(demresfel, resd8file)
    hydro.scad8(0, resd8file, resscafile)
    #hydro.scadinf(resdirfile, resscafile)

    #hydro.PeukerDouglas(demresfel, weifile)
    #hydro.scad8(0, resdirfile, ssafile)
    #thresh = hydro.dropanalysis(10, 2, 500, demresfel, resdirfile)
    #print(thresh)
    hydro.extractnet(resscafile, resnet, threshold)
    hydro.strahler(resnet, resd8file, res_so)
    hydro.watershed_py(res_so, resd8file, res_ws)
    streamSHP(res_so, netshp, resd8file, res_ws)



def LargeHydroAnalysis(_g, buf, nrank, threshold = 500000):
    hydro.nrank = nrank
    #largePitremove(_g, buf, res_ws, dem, demfel)
    hydro.flowdird8(demfel, dirfile)
    largearea(_g, buf, res_ws, dirfile, scafile)
    hydro.threshold(threshold, scafile, netfile)
    hydro.strahler(netfile, dirfile, streamfile)
    hydro.watershed_py(streamfile, dirfile, wsfile)
    #hydro.streamSHP(streamfile, netshp, dirfile, wsfile)

if __name__ == "__main__":
    '''
    Large are hydroAnalysis based on Pargo:
    First round: 将dem重采样至较粗分辨率后进行第一轮快速分割
    Second round: 基于第一轮的粗分割子流域，对每个子流域进行流域分割
    '''

    _g = 5  # 重采样倍数
    buf = 2  # 缓冲区大小
    nrank= 16  # 进程数
    mpi = r"mpiexec -n " + str(nrank)
    print(workdir)
    print("rank number:", nrank)
    dem = r"data\\dem.tif "  # 输入DEM

    threshold = 20000 #自定义第一轮粗流域分割河网提取阈值
    threshold2 = 200000

    #hydroAnalysis(_g, buf, nrank, threshold) #第一轮粗分割

    #两轮最好分开运行，调整阈值确定第一轮的粗分割子流域数量，粗分割子流域数量由数据量大小与计算能力共同决定，50个左右为宜
    LargeHydroAnalysis(_g, buf, nrank, threshold2)#第二轮细分割






