
from osgeo import gdal, ogr, osr
import os

def streamSHP(input_file, output_file, dir_file, watershed_file):
    """
    将栅格数据转换为线要素shp文件
    :param input_file: 输入栅格文件路径
    :param output_file: 输出shp文件路径
    """
    # 打开栅格数据集
    src_ds = gdal.Open(input_file)
    # 获取栅格波段
    srcband = src_ds.GetRasterBand(1)
    # 读取栅格数据为NumPy数组
    data = srcband.ReadAsArray()
    noData = srcband.GetNoDataValue()

    # 打开流向数据集
    dir_ds = gdal.Open(dir_file)
    # 获取流向波段
    dirband = dir_ds.GetRasterBand(1)
    # 读取流向数据为NumPy数组
    dir_data = dirband.ReadAsArray()
    # 打开流域数据集
    watershed_ds = gdal.Open(watershed_file)
    # 获取流域波段
    watershed_band = watershed_ds.GetRasterBand(1)
    # 读取流域数据为NumPy数组
    watershed_data = watershed_band.ReadAsArray()

    # 创建矢量数据源
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(output_file):
        driver.DeleteDataSource(output_file)
    dst_ds = driver.CreateDataSource(output_file)

    # 创建图层
    srs = osr.SpatialReference()
    srs.ImportFromWkt(src_ds.GetProjectionRef())
    dst_layer = dst_ds.CreateLayer("stream", srs=srs, geom_type=ogr.wkbLineString)

    # 创建字段
    fd = ogr.FieldDefn("Value", ogr.OFTInteger)
    dst_layer.CreateField(fd)
    fd = ogr.FieldDefn("Watershed", ogr.OFTInteger)
    dst_layer.CreateField(fd)
    # 获取栅格数据的仿射变换参数
    gt = src_ds.GetGeoTransform()

    def pixel2coord(x, y):
        """将像素坐标转换为地理坐标"""
        xp = gt[0] + x * gt[1] + y * gt[2]
        yp = gt[3] + x * gt[4] + y * gt[5]
        return xp, yp

    def get_upstream(x, y, value):
        """获取上游像素点坐标"""
        upstream = []
        for index, (nx, ny) in enumerate(
            [(x - 1, y - 1), (x, y - 1), (x + 1, y - 1), (x - 1, y), (x, y), (x + 1, y), (x - 1, y + 1), (x, y + 1),
             (x + 1, y + 1)]):
            if (0 <= nx < data.shape[1] and 0 <= ny < data.shape[0]) and (data[ny, nx] == value or data[ny, nx] == -1):
                d = dir_data[ny, nx]
                if index + d == 8 and index != 4 :
                    upstream.append((nx, ny))
        return upstream

    def get_downstream(x, y):
        d = dir_data[y, x]
        if d == 0:
            nx, ny = x - 1, y - 1
        elif d == 1:
            nx, ny = x, y - 1
        elif d == 2:
            nx, ny = x + 1, y - 1
        elif d == 3:
            nx, ny = x - 1, y
        elif d == 5:
            nx, ny = x + 1, y
        elif d == 6:
            nx, ny = x - 1, y + 1
        elif d == 7:
            nx, ny = x, y + 1
        elif d == 8:
            nx, ny = x + 1, y + 1
        return nx, ny

    def track_line(x, y, value):
        """追踪线要素"""
        line = []
        stack= [(x, y)]
        while stack:
            x, y = stack.pop()
            # 将当前像素点添加到线要素中
            line.append((x, y))
            up = get_upstream(x, y, value)
            if len(up) <= 1:
            # 将当前像素点设置为NoData值，避免重复追踪
                data[y, x] = -1
            #else:


            nx, ny = get_downstream(x, y)

            if 0 <= nx < data.shape[1] and 0 <= ny < data.shape[0] and len(up)<=1:
                if data[ny, nx] == value and len(get_upstream(nx, ny, value)) <= 1:
                    stack.append((nx, ny))
                if data[ny, nx] == value and len(get_upstream(nx, ny, value)) > 1:
                    line.append((nx, ny))
                    line.reverse()
                if data[ny, nx] > value:
                    line.append((nx, ny))
                    line.reverse()

            if len(up) == 1:
                nx, ny = up[0]
                if data[ny, nx] == value:
                    stack.append(up[0])
            else:
                ws = watershed_data[y, x]
                line.reverse()


        return ws, line

    def create_feature(coords, value, watershed_value):
        """创建矢量要素"""
        if len(coords) > 1:
            feature = ogr.Feature(dst_layer.GetLayerDefn())
            line = ogr.Geometry(ogr.wkbLineString)
            for x, y in coords:
                xp, yp = pixel2coord(x, y)
                line.AddPoint(xp, yp)
            feature.SetGeometry(line)
            feature.SetField("Value", value)

            feature.SetField("Watershed", watershed_value)
            dst_layer.CreateFeature(feature)

    # 遍历栅格数据，提取线要素
    for value in [1, 2, 3, 4]:
        for y in range(data.shape[0]):
            for x in range(data.shape[1]):
                if data[y, x] == value:
                    ws, line_coords = track_line(x, y, value)
                    create_feature(line_coords, value, ws)

    # 关闭数据源
    dst_ds.Destroy()
    src_ds = None
    dir_ds = None
    watershed_ds = None
