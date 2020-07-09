import gdal
import ogr
import sys
import numpy as np
from datetime import datetime
import dbf
import os
import os.path
import glob


input_files = ("/home/tepex/NIERSC/IEPI/20200708/radarsat/RS2_OK000000_PK0829186_DK0200708_SCNB_20200708_011434_HHHV_SGF_map_w64_gtiff.tif",
"/home/tepex/NIERSC/IEPI/20200708/radarsat/RS2_OK000000_PK0829186_DK0200708_SCNB_20200708_011514_HHHV_SGF_map_w64_gtiff.tif")


def prepare_radarsat_single(input_file, output_dir):
    scene_date = input_file.split("_")[5]
    scene_time = input_file.split("_")[6]
    # To WGS84 tif
    out_tif = "niersc_ice_class_%sT%s.tif" % (scene_date, scene_time)
    cmd = 'gdalwarp -t_srs EPSG:4326 %s %s' % (input_file, os.path.join(output_dir, out_tif))
    print cmd
    os.system(cmd)
    return os.path.join(output_dir, out_tif)


def merge_radarsat(input_file_list, out_dir):
    if len(input_file_list) > 0:
        input_file_list.sort(key=lambda x: os.path.basename(x).split("_")[3].split("T")[0])  # sorts products by their date from oldest to newest
        newest_date = os.path.basename(input_file_list[-1]).split("_")[3].split("T")[0]  # returns YYmmdd of the newest product
        out_path = os.path.join(out_dir, "niersc_ice_class_rs_%s_temp.tif" % newest_date)
        cmd = "gdal_merge.py -n 0 -a_nodata 0 -o %s" % out_path
        for input_file in input_file_list:
            cmd += " %s" % input_file
        print cmd
        os.system(cmd)
        return out_path
    else:
        return False


def prepare_radarsat_aggregation(input_file, out_dir):
    classification_date = input_file.split(".")[0].split("_")[4]
    print classification_date

    # To edge WGS84 tif
    ice_edge_wgs84_tif_name = 'niersc_ice_edge_rs_%s.tif' % classification_date
    cmd = 'gdal_calc.py -A %s --outfile=%s --NoDataValue=255 --calc="(A==0)*255 + (A==1)*255 + (A==2)*0 + (A==3)*1 + (A==4)*1 + (A==5)*0 + (A==6)*0 + (A==7)*0"' % (input_file, os.path.join(out_dir, ice_edge_wgs84_tif_name))
    print cmd
    os.system(cmd)

    # edge WGS84tif to coloured render
    fname = os.path.join(out_dir, ice_edge_wgs84_tif_name)
    cmd = 'gdaldem color-relief -alpha %s %s %s' % (fname, 'ice_edge_render.txt', "%s_render_tmp.tif" % fname[:-4])
    print cmd
    os.system(cmd)

    # render produced by gdaldem will be huge - let's deflate it:
    cmd = 'gdal_translate -co "COMPRESS=DEFLATE" %s %s' % ("%s_render_tmp.tif" % fname[:-4], "%s_render.tif" % fname[:-4])
    os.system(cmd)
    os.remove("%s_render_tmp.tif" % fname[:-4])

    # To ice_edge vector
    ice_edge_shp_temp_name = 'niersc_ice_edge_rs_' + classification_date + '_temp.shp'
    ice_edge_shp_name = 'niersc_ice_edge_rs_' + classification_date + '.shp'
    cmd = 'gdal_polygonize.py %s -f "ESRI Shapefile" %s DN' % (os.path.join(out_dir, ice_edge_wgs84_tif_name), os.path.join(out_dir, ice_edge_shp_temp_name))
    print cmd
    os.system(cmd)

    # deleting water polygons, leaving only ice ones:
    cmd = 'ogr2ogr -f "ESRI Shapefile" -where "DN=1" %s %s' % (os.path.join(out_dir, ice_edge_shp_name), os.path.join(out_dir, ice_edge_shp_temp_name))
    print cmd
    os.system(cmd)

    print "Adding DATE_TIME attribute into ice edge..."
    ice_edge_dbf_name = os.path.join(out_dir, ice_edge_shp_name.split('.')[0] + '.dbf')
    ice_edge_db = dbf.Table(ice_edge_dbf_name)
    with ice_edge_db:
        # ice_edge_db.add_fields('date_time C(100)')
        ice_edge_db.add_fields('DATE_TIME D')
        for record in ice_edge_db:
            # dbf.write(record, date_time=classification_date.strftime('%Y%m%d'))
            # classification_date = datetime.strptime(classification_date, "%Y%m%d")
            # print type(classification_date)
            class_date = datetime.strptime(classification_date, "%Y%m%d")
            dbf.write(record, date_time=class_date)

    ##########################################
    # Emulating pseudo ice classification tif:
    input_file = os.path.join(out_dir, ice_edge_wgs84_tif_name)
    ice_class_tif_name = 'niersc_ice_class_rs_%s.tif' % classification_date
    output_file = os.path.join(out_dir, ice_class_tif_name)
    cmd = 'gdal_calc.py -A %s --outfile=%s --NoDataValue=255 --calc="(A==0)*3 + (A==1)*10"' % (input_file, output_file)
    print cmd
    os.system(cmd)

    # ice class WGS84tif to coloured render
    input_file = os.path.join(out_dir, ice_class_tif_name)
    cmd = 'gdaldem color-relief -alpha %s %s %s' % (input_file, 'ice_class_summer_render.txt', "%s_render_tmp.tif" % input_file[:-4])
    print cmd
    os.system(cmd)
    cmd = 'gdal_translate -co "COMPRESS=DEFLATE" %s %s' % ("%s_render_tmp.tif" % input_file[:-4], "%s_render.tif" % input_file[:-4])
    os.system(cmd)
    os.remove("%s_render_tmp.tif" % input_file[:-4])

    # converting emulated ice class raster into vector
    input_file = os.path.join(out_dir, ice_class_tif_name)
    ice_class_shp_name = 'niersc_ice_class_rs_' + classification_date + '.shp'
    output_file = os.path.join(out_dir, ice_class_shp_name)
    cmd = 'gdal_polygonize.py %s -f "ESRI Shapefile" %s CLASS' % (input_file, output_file)
    print cmd
    os.system(cmd)

    print "Adding DATE_TIME attribute into ice class..."
    ice_class_dbf_name = os.path.join(out_dir, ice_class_shp_name.split('.')[0] + '.dbf')
    ice_class_db = dbf.Table(ice_class_dbf_name)
    with ice_class_db:
        ice_class_db.add_fields('DATE_TIME D')
        for record in ice_class_db:
            class_date = datetime.strptime(classification_date, "%Y%m%d")
            dbf.write(record, date_time=class_date)

    # clean-up:
    print 'deleting temp'
    temp_files = glob.glob(os.path.join(out_dir, '*temp*'))
    for temp_file in temp_files:
        os.remove(temp_file)


if __name__ == "__main__":

    single_rs_list = []
    for input_file in input_files:
        output_dir = os.path.dirname(input_file)
        single_rs_list.append(prepare_radarsat_single(input_file, output_dir))

    output_dir = "/home/tepex/NIERSC/IEPI/20200708/radarsat/"
    mosaic = merge_radarsat(single_rs_list, output_dir)
    prepare_radarsat_aggregation(mosaic, output_dir)
