import math
import sys
import time
import os
import rasterio
import numpy as np
from streamRasToShp import streamSHP
cur_dir = os.getcwd()
currentlimt=sys.getrecursionlimit()
newlim=50000
sys.setrecursionlimit(newlim)
#workdir = r"G: & cd egc\\pargo-dev\\wkflow & "
workdir = "cd wkflow & "

work = r"G:\\egc\\pargo-dev\\wkflow\\"  # 文件存放位置
nrank = 8 # 进程数
mpi = r"mpiexec -n " + str(nrank)
nbr = r"data\\neighbor\\moore.nbr "
middata = r"data\\process\\"
finaldata = r"data\\final\\"
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

def resample(dem, resdem, _g):
    exe = r" apps\resample.exe "
    cmd = workdir + mpi + exe + "-input " + dem + "-output " + resdem + "-g " + str(_g)
    print("working process: resample(g=", _g)
    start = time.time()
    g = os.system(cmd)
    end = time.time()
    print("process time:", end - start)

def pitremove(dem, demfel):
    exe = r" apps\pitremove.exe "
    cmd = workdir + mpi + exe + dem + nbr + demfel
    print(cmd)
    print("working process: PitRemove of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def flowdird8(dem, dir):
    exe = r" apps\flowdird8.exe "
    cmd = workdir + mpi + exe + dem + nbr + dir + str(8)
    print("working process: calculate d8 of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def flowdirdinf(dem, dinf):
    exe = r" apps\flowdirdinf.exe "
    cmd = workdir + mpi + exe + dem + nbr + dinf
    print("working process: calculate d-inf of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def scadinf(dinf, scadinf):
    exe = r" apps\scadinf.exe "
    cmd = workdir + mpi + exe + dinf + nbr + scadinf
    print("working process: calculate sca of ", dinf)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def scad8(w, dir, sca, wei=weifile):
    exe = r" apps\scad8.exe "
    cmd = workdir + mpi + exe + dir + nbr + sca
    wcmd = workdir + mpi + exe + dir + nbr + wei + sca
    start = time.time()
    if (w):
        print("working process: calculate weighted_sca of ", dir)
        start = time.time()
        g = os.system(wcmd)
    else:
        print("working process: calculate sca of ", dir)
        start = time.time()
        g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def PeukerDouglas(dem, weifile):
    exe = r" apps\PeukerDouglas.exe "
    cmd = workdir + mpi + exe + dem + nbr + weifile
    print("working process: catch possible river net of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def dropanalysis(nthresh, threshmin, threshmax, dem, dir):
    exe = r" apps\dropanalysis.exe "
    cmd = workdir + mpi + exe + dem + dir + resscafile + ssafile + nbr
    print("working process: drop analysis of ", dem)
    start = time.time()
    t = threshmin
    for th in range(0, nthresh):
        r = math.exp((math.log(threshmax) - math.log(threshmin)) / (nthresh - 1))
        thresh = threshmin * pow(r, th)
        c = cmd + str(thresh)
        g = os.system(c)
        print(g)
        if (g == 0):
            t = thresh
            break
    end = time.time()
    print("process time:", end - start)
    print("thresh:",t)
    return t


def threshold(thresh, sca, net):
    exe = r" apps\threshold.exe "
    cmd = workdir + mpi + exe + sca + net + str(thresh)
    print("working process: extract river net with threshold:", thresh)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def stream(net, dir, sm):
    exe = r" apps\stream.exe "
    cmd = workdir + mpi + exe + net + dir + nbr + sm
    print("working process: calculate stream order of ", net)
    print(cmd)
    g = os.system(cmd)
    print(g)

def watershed(stream, dir, ws):
    exe = r" apps\watershed.exe "
    cmd = workdir + mpi + exe + stream + dir + nbr + ws
    print(cmd)
    print("working process: calculate stream order of ", stream)
    g = os.system(cmd)
    print(g)

def dinf2d8(dinffile, d8file):
    with rasterio.open(dinffile) as src:
        dinf = src.read(1)
        _noData = src.nodata
        profile = src.profile
    net = np.full(dinf.shape, np.nan)
    Eps = 0.01
    iNeighborCells = 1

    for iRow in range(iNeighborCells, dinf.shape[0] - iNeighborCells):
        for iCol in range(iNeighborCells, dinf.shape[1] - iNeighborCells):
            if not math.isclose(sca[iRow][iCol], _noData, abs_tol=Eps):
                if sca[iRow][iCol] >= th:
                    net[iRow][iCol] = 1
                else:
                    net[iRow][iCol] = 0
            else:
                net[iRow][iCol] = _noData

    # 更新nodata值
    profile.update(nodata=_noData)
        # 写入新的tif文件
    with rasterio.open(netfile, 'w', **profile) as dst:
        dst.write(net, 1)


def extractnet(scafile, netfile, th):
    with rasterio.open(scafile) as src:
        sca = src.read(1)
        _noData = src.nodata
        profile = src.profile

    # 创建一个空数组来存储分级结果
    net = np.full(sca.shape, np.nan)
    Eps = 0.01
    iNeighborCells=1

    for iRow in range(iNeighborCells, net.shape[0] - iNeighborCells):
        for iCol in range(iNeighborCells, net.shape[1] - iNeighborCells):
            if not math.isclose(sca[iRow][iCol], _noData, abs_tol=Eps):
                if sca[iRow][iCol] >= th:
                    net[iRow][iCol] = 1
                else:
                    net[iRow][iCol] = 0
            else:
                net[iRow][iCol] = _noData

    # 更新nodata值
    profile.update(nodata=_noData)
        # 写入新的tif文件
    with rasterio.open(netfile, 'w', **profile) as dst:
        dst.write(net, 1)

def strahler(netfile, dirfile, sm):
    with rasterio.open(netfile) as src:
        net = src.read(1)
        _noData = src.nodata
        profile = src.profile
    with rasterio.open(dirfile) as src:
        dir = src.read(1)
    # 创建一个空数组来存储分级结果
    orderout = np.full(net.shape, np.nan)
    d1 = [-1, -1, -1, 0, 0, 0, 1, 1, 1]
    d2 = [-1, 0, 1, -1, 0, 1, -1, 0, 1]
    Eps = 0.01
    iNeighborCells=1
    contribs = np.zeros_like(net)
    #orderout = np.zeros_like(net)
    for iRow in range(iNeighborCells, net.shape[0] - iNeighborCells):
        for iCol in range(iNeighborCells, net.shape[1] - iNeighborCells):
            if abs(net[iRow, iCol] - _noData) > Eps and net[iRow, iCol] == 1:

                # 遍历所有邻居
                for tRow in range(iRow - iNeighborCells, iRow + 2 * iNeighborCells):
                    for tCol in range(iCol - iNeighborCells, iCol + 2 * iNeighborCells):
                        if abs(dir[tRow, tCol] - _noData) >= Eps and dir[tRow, tCol] >= 0:
                            d = int(dir[tRow, tCol])
                            if d1[d] + tRow == iRow and d2[d] + tCol == iCol and d + dir[iRow, iCol] != 8 and d != 4:
                                if net[tRow, tCol] == 1:
                                    contribs[iRow, iCol] += 1

            else:
                contribs[iRow, iCol] = -1
    for iRow in range(iNeighborCells, net.shape[0] - iNeighborCells):
        for iCol in range(iNeighborCells, net.shape[1] - iNeighborCells):
            if abs(net[iRow, iCol] - _noData) > Eps and net[iRow, iCol] == 1:
                orderdec(orderout, contribs,dir, iRow, iCol, _noData)
            else:
                orderout[iRow, iCol] = _noData
    # 更新nodata值
    profile.update(nodata=_noData)
        # 写入新的tif文件
    with rasterio.open(sm, 'w', **profile) as dst:
        dst.write(orderout, 1)

def newOrder(nOrder):
    junction = False
    source = True
    oOut = 1
    ordermax = 0
    count = 0
    for i in range(9):
        if nOrder[i] > 0:
            count += 1
            source = False
            if count == 1:
                oOut = nOrder[i]
                ordermax = nOrder[i]
            else:
                if nOrder[i] > oOut:
                    ordermax = nOrder[i]
                    oOut = nOrder[i]
                elif nOrder[i] == oOut:
                    oOut = ordermax + 1
                junction = True
    return oOut, junction, source

def orderdec(orderout, contribs, dir, iRow, iCol, _noData):
    d1 = [-1, -1, -1, 0, 0, 0, 1, 1, 1]
    d2 = [-1, 0, 1, -1, 0, 1, -1, 0, 1]
    Eps = 0.01
    if not contribs[iRow, iCol]:
        nOrder = np.zeros(9, dtype=np.int16)
        k = 0
        for tRow in range(iRow - 1, iRow + 1 + 1):
            for tCol in range(iCol - 1, iCol + 1 + 1):
                if dir[tRow, tCol] >= 0 and abs(dir[tRow, tCol] - _noData) > Eps:
                    d = int(dir[tRow, tCol])
                    if d1[d] + tRow == iRow and d2[d] + tCol == iCol and d + dir[
                        iRow, iCol] != 8 and d != 4 and abs(orderout[tRow, tCol] - _noData) > Eps:
                        nOrder[k] = orderout[tRow, tCol]
                k += 1
        oOut, junction, source = newOrder(nOrder)
        orderout[iRow, iCol] = oOut
        contribs[iRow, iCol] = -1
        if dir[iRow, iCol] >= 0 and abs(dir[iRow, iCol] - _noData) > Eps:
            d =int(dir[iRow, iCol])
            tRow = iRow + d1[d]
            tCol = iCol + d2[d]

            if contribs[tRow, tCol] > 0:
                contribs[tRow, tCol] -= 1
            if contribs[tRow, tCol] == 0:
                orderdec(orderout, contribs, dir, tRow, tCol, _noData)

def watershed_py(streamfile, dirfile, wsfile):
    with rasterio.open(streamfile) as src:
        str = src.read(1)
        _noData = src.nodata
        profile = src.profile
    with rasterio.open(dirfile) as src:
        dir = src.read(1)
    # 创建一个空数组来存储分级结果
    wsout = np.full(str.shape, np.nan)
    d1 = [-1, -1, -1, 0, 0, 0, 1, 1, 1]
    d2 = [-1, 0, 1, -1, 0, 1, -1, 0, 1]
    Eps = 0.01
    iNeighborCells = 1
    n1=0
    n2=0
    contribs = np.zeros_like(str)
    # orderout = np.zeros_like(net)
    for iRow in range(iNeighborCells, str.shape[0] - iNeighborCells):
        for iCol in range(iNeighborCells, str.shape[1] - iNeighborCells):
            if abs(dir[iRow, iCol] - _noData) > Eps and dir[iRow, iCol] >=0 :

                # 遍历所有邻居
                for tRow in range(iRow - iNeighborCells, iRow + 2 * iNeighborCells):
                    for tCol in range(iCol - iNeighborCells, iCol + 2 * iNeighborCells):
                        if abs(dir[tRow, tCol] - _noData) >= Eps and dir[tRow, tCol] >= 0:
                            d = int(dir[tRow, tCol])
                            if d1[d] + tRow == iRow and d2[d] + tCol == iCol and d + dir[iRow, iCol] != 8 and d != 4:
                                if abs(str[tRow, tCol] - _noData) >= Eps and str[tRow, tCol] >= 1:
                                    contribs[iRow, iCol] += 1

            else:
                contribs[iRow, iCol] = -2
    for iRow in range(iNeighborCells, str.shape[0] - iNeighborCells):
        for iCol in range(iNeighborCells, str.shape[1] - iNeighborCells):
            if dir[iRow, iCol] >= 0 and abs(dir[iRow, iCol] - _noData) > Eps:
                d = int(dir[iRow, iCol])
                #print(abs(dir[iRow, iCol] - _noData))
                tRow = iRow + d1[d]
                tCol = iCol + d2[d]
                if( contribs[iRow][iCol]>0 and contribs[tRow][tCol]==-2):
                    n2 = n2 + 1
                if (contribs[iRow][iCol] > 0 and contribs[tRow][tCol] >= 2):
                    n1 = n1 + 1
    #n = n2+n1
    n=0
    for iRow in range(iNeighborCells, str.shape[0] - iNeighborCells):
        for iCol in range(iNeighborCells, str.shape[1] - iNeighborCells):
            if abs(dir[iRow, iCol] - _noData) >= Eps and dir[iRow, iCol] >= 0 and contribs[iRow,iCol]>=0:
                d = int(dir[iRow, iCol])
                n = watersheddec(d, iRow, iCol, wsout, contribs, dir, n, _noData)
    # 更新nodata值
    profile.update(nodata=_noData)
        # 写入新的tif文件
    with rasterio.open(wsfile, 'w', **profile) as dst:
        dst.write(wsout, 1)

def watersheddec(d, iRow, iCol, ws, contribs, dir, n2, _noData):
    d1 = [-1, -1, -1, 0, 0, 0, 1, 1, 1]
    d2 = [-1, 0, 1, -1, 0, 1, -1, 0, 1]
    tRow = iRow + d1[d]
    tCol = iCol + d2[d]
    Eps = 0.01


    if (contribs[iRow][iCol] > 0 and (contribs[tRow][tCol] == -2 or contribs[tRow][tCol] == -3 or contribs[tRow][tCol] > 1)):
        ws[iRow][iCol]=n2+1
        n2=n2+1
        contribs[iRow][iCol]=-1

    elif (contribs[tRow][tCol] == -1 or contribs[tRow][tCol] == -3):
        ws[iRow][iCol]=ws[tRow][tCol]
        if (contribs[iRow][iCol] < 2):
            contribs[iRow][iCol]=-1
        else:
            contribs[iRow][iCol]=-3
    elif (contribs[iRow][iCol] == 0 and contribs[tRow][tCol] == -2):
        ws[iRow][iCol]=ws[tRow][tCol]
        contribs[iRow][iCol]=-2
    if abs(ws[iRow, iCol] - _noData) > Eps:
        for tRow in range(iRow - 1, iRow + 1 + 1):
            for tCol in range(iCol - 1, iCol + 1 + 1):
                if dir[tRow, tCol] >= 0 and abs(dir[tRow, tCol] - _noData) > Eps:
                    d = int(dir[tRow, tCol])
                    if d1[d] + tRow == iRow and d2[d] + tCol == iCol and d + dir[
                    iRow, iCol] != 8 and d != 4:
                        n2 = watersheddec(d, tRow, tCol, ws, contribs, dir, n2, _noData)
    return n2

def hydroAnalysis():

    pitremove(dem, demfel)
    flowdird8(demfel, dirfile)
    scad8(0, dirfile, scafile)

    threshold(800000, scafile, netfile)
    strahler(netfile, dirfile, strahlerfile)
    watershed_py(strahlerfile, dirfile, wsfile)
    streamSHP(strahlerfile, netshp, dirfile, wsfile)


if __name__ == "__main__":


    nrank = 1  # 进程数
    mpi = "mpiexec -n " + str(nrank)
    dem = "data\\XJdem.tif "  # 输入DEM
    scafile = "data\\sca.tif"
    netfile = "data\\net.tif"
    dirfile = "data\\d8.tif"
    strahlerfile = "data\\str.tif"
    wsfile = "data\\watershed.tif"
    #hydroAnalysis()

    threshold(800000, scafile, netfile)
    strahler(netfile, dirfile, strahlerfile)
    watershed_py(strahlerfile, dirfile, wsfile)
    streamSHP(strahlerfile, netshp, dirfile, wsfile)

