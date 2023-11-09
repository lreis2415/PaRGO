import math
import time
import os
import rasterio
import numpy as np
cur_dir = os.getcwd()

#workdir = r"G: & cd egc\\pargo-dev\\wkflow & "
workdir = "cd wkflow & "

work = r"G:\\egc\\pargo-dev\\wkflow\\"  # 文件存放位置
nrank = 16  # 进程数
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

def resample(dem, resdem, _g):
    exe = r" D:\EGC\PaRGO-fxc\build\apps\spatial\Release\resample.exe "
    cmd = workdir + mpi + exe + "-input " + dem + "-output " + resdem + "-g " + str(_g)
    print("working process: resample(g=", _g)
    start = time.time()
    g = os.system(cmd)
    end = time.time()
    print("process time:", end - start)

def pitremove(dem, demfel):
    exe = r" D:\EGC\PaRGO-dev-fxc\build\apps\hydrology\Release\pitremove.exe "
    cmd = workdir + mpi + exe + dem + nbr + demfel
    print("working process: PitRemove of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def flowdird8(dem, dir):
    exe = r" D:\EGC\PaRGO-dev-fxc\build\apps\hydrology\Release\flowdird8.exe "
    cmd = workdir + mpi + exe + dem + nbr + dir + str(8)
    print("working process: calculate d8 of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def flowdirdinf(dem, dinf):
    exe = r" D:\EGC\PaRGO-dev-fxc\build\apps\hydrology\Release\flowdirdinf.exe "
    cmd = workdir + mpi + exe + dem + nbr + dinf
    print("working process: calculate d-inf of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def scadinf(dinf, scadinf):
    exe = r" D:\EGC\PaRGO-dev-fxc\build\apps\hydrology\Release\scadinf.exe "
    cmd = workdir + mpi + exe + dinf + nbr + scadinf
    print("working process: calculate sca of ", dinf)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def scad8(w, dir, sca, wei=weifile):
    exe = r" D:\EGC\PaRGO-dev-fxc\build\apps\hydrology\Release\scad8.exe "
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
    exe = r" D:\EGC\PaRGO-dev-fxc\build\apps\hydrology\Release\PeukerDouglas.exe "
    cmd = workdir + mpi + exe + dem + nbr + weifile
    print("working process: catch possible river net of ", dem)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def dropanalysis(nthresh, threshmin, threshmax, dem, dir):
    exe = r" D:\EGC\PaRGO-dev-fxc\build\apps\hydrology\Release\dropanalysis.exe "
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
    exe = r" D:\EGC\PaRGO-fxc\build\apps\hydrology\Release\threshold.exe "
    cmd = workdir + mpi + exe + sca + net + str(thresh)
    print("working process: extract river net with threshold:", thresh)
    start = time.time()
    g = os.system(cmd)
    print(g)
    end = time.time()
    print("process time:", end - start)

def stream(net, dir, sm):
    exe = r" D:\EGC\PaRGO-fxc\build\apps\hydrology\Release\stream.exe "
    cmd = workdir + mpi + exe + net + dir + nbr + sm
    print("working process: calculate stream order of ", net)
    print(cmd)
    g = os.system(cmd)
    print(g)

def watershed(stream, dir, ws):
    exe = r" D:\EGC\PaRGO-fxc\build\apps\hydrology\Release\watershed.exe "
    cmd = workdir + mpi + exe + stream + dir + nbr + ws
    print(cmd)
    print("working process: calculate stream order of ", stream)
    g = os.system(cmd)
    print(g)

def strahler(netfile, dirfile, sm):
    with rasterio.open(netfile) as src:
        net = src.read(1)
        _noData = src.nodata
        profile = src.profile
    with rasterio.open(dirfile) as src:
        dir = src.read(1)
    # 创建一个空数组来存储分级结果
    orderout = np.full(net.shape, np.nan)
    d1 = {-1, -1, -1, 0, 0, 0, 1, 1, 1}
    d2 = {-1, 0, 1, -1, 0, 1, -1, 0, 1}
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
                            d = dir[tRow, tCol]
                            if d1[d] + tRow == iRow and d2[d] + tCol == iCol and d + dir[iRow, iCol] != 8 and d != 4:
                                if net[tRow, tCol] == 1:
                                    contribs[iRow, iCol] += 1

            else:
                contribs[iRow, iCol] = -1
    for iRow in range(iNeighborCells, net.shape[0] - iNeighborCells):
        for iCol in range(iNeighborCells, net.shape[1] - iNeighborCells):
            if abs(net[iRow, iCol] - _noData) > Eps and net[iRow, iCol] == 1:
                orderdec(orderout, contribs, iRow, iCol, _noData)
            else:
                orderout[iRow, iCol] = _noData
    # 更新nodata值
    profile.update(nodata=np.nan)
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

def orderdec(orderout, contribs, iRow, iCol, _noData):
    d1 = {-1, -1, -1, 0, 0, 0, 1, 1, 1}
    d2 = {-1, 0, 1, -1, 0, 1, -1, 0, 1}
    Eps = 0.01
    if not contribs[iRow, iCol]:
        nOrder = np.zeros(9, dtype=np.int16)
        k = 0
        for tRow in range(iRow - 1, iRow + 1 + 1):
            for tCol in range(iCol - 1, iCol + 1 + 1):
                if dir[tRow, tCol] >= 0 and abs(dir[tRow, tCol] - _noData) < Eps:
                    d = dir[tRow, tCol]
                    if d1[d] + tRow == iRow and d2[d] + tCol == iCol and d + dir[
                        iRow, iCol] != 8 and d != 4 and abs(orderout[tRow, tCol] - _noData) < Eps:
                        nOrder[k] = orderout[tRow, tCol]
                k += 1
        oOut, junction, source = newOrder(nOrder)
        orderout[iRow, iCol] = oOut
        contribs[iRow, iCol] = -1
        if dir[iRow, iCol] >= 0 and abs(dir[iRow, iCol] - _noData) < Eps:
            d = dir[iRow, iCol]
            tRow = iRow + d1[d]
            tCol = iCol + d2[d]
            if contribs[tRow, tCol] > 0:
                contribs[tRow, tCol] -= 1
            if contribs[tRow, tCol] == 0:
                orderdec(orderout, contribs, tRow, tCol, _noData)

if __name__ == "__main__":
    strahler("data\\", dirfile, sm)
