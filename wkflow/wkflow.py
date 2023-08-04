# coding:utf-8

import geopandas as gpd
import queue
import math
import time
import os
from streamRasToShp import streamSHP
import hydroAnalysis as hydro
cur_dir = os.getcwd()

#workdir = r"G: & cd egc\\pargo-dev\\wkflow & "
workdir = "cd wkflow & "

work = r"G:\\egc\\pargo-dev\\wkflow\\"  # 文件存放位置
nrank = 8  # 进程数
mpi = r"mpiexec -n " + str(nrank)
middata = r"data\\process\\"
finaldata = r"data\\final\\"
nbr = r"data\\neighbor\\moore.nbr "
dem = r"data\\dem.tif "

demres = middata + r"res_dem.tif "
demresfel = middata + r"res_demfel.tif "
resdirfile = middata + r"res_dir.tif "
resscafile = middata + r"res_sca.tif "
weifile = middata + r"weight.tif "
ssafile = middata + r"w_sca.tif "
resnet = middata + r"res_net.tif "
res_so = middata + r"res_stream.tif "
res_ws = middata + r"res_watershed.tif "
netshp = middata + r"res_rivernet.shp"

demfel = finaldata + r"demfel.tif "
dirfile = finaldata + r"dir.tif "
scafile = finaldata + r"sca.tif "
netfile = finaldata + r"net.tif "
streamfile = finaldata + r"streamOrder.tif "









def largearea(_g, buf, watershed, dir, sca):
    exe = r" D:\EGC\PaRGO\build\apps\hydrology\Release\Largearead8.exe "
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
    exe = r" D:\EGC\PaRGO\build\apps\hydrology\Release\Largepitremove.exe "

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



def hydroAnalysis(_g, buf):
    hydro.resample(dem, demres, _g)
    hydro.pitremove(demres, demresfel)
    hydro.flowdird8(demresfel, resdirfile)
    hydro.scad8(0, resdirfile, resscafile)
    hydro.PeukerDouglas(demresfel, weifile)
    hydro.scad8(1, resdirfile, ssafile, weifile)
    thresh = hydro.dropanalysis(10, 2, 500, demresfel, resdirfile)
    #print(thresh)
    hydro.threshold(thresh * 1000, resscafile, resnet)
    hydro.stream(resnet, resdirfile, res_so)
    hydro.watershed(res_so, resdirfile, res_ws)
    #largePitremove(_g, buf)
    #flowdird8(demfel, dirfile)
    #t = largearea(_g, buf)
    #print("process time:", t)
    #threshold(5000, scafile, netfile)
    #stream(netfile, dirfile, streamfile)


if __name__ == "__main__":
    _g = 10  # 重采样倍数
    buf = 2  # 缓冲区大小
    nrank= 12  # 进程数
    mpi = r"mpiexec -n " + str(nrank)
    print(workdir)
    print("rank number:", nrank)
    dem = r"data\\dem.tif "  # 输入DEM
    #hydroAnalysis(_g, buf)
    #streamSHP(res_so, netshp, resdirfile, res_ws)
    #largePitremove(_g, buf, res_ws, dem, demfel)
    hydro.flowdird8(demfel, dirfile)
    largearea(_g, buf, res_ws, dirfile, scafile)






