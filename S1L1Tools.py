# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 13:15:03 2017

@author: silent
@author: antonv
"""
import gdal
import zipfile
import os
import json
import xml.etree.ElementTree as etree
from datetime import datetime
from scipy import interpolate, ndimage
import numpy as np
# import shapely.wkt
import geojson
import glob
import matplotlib.pyplot as plt


class S1L1Tools:
    def __init__(self, datasource_path):
        self.datasource_path = datasource_path
        self.datasource_body = os.path.basename(datasource_path).split('.')[0]
        self.measurements = []
        self.radiometric_correction_LUT = []
        self.noise_correction_LUT = []

        self.metadata = {}
        self.__get_metadata()

        for polarisation in self.metadata['polarisations']:
            for archived_file in self.file_list:
                if (archived_file.filename.lower().find(polarisation.lower()) != -1) and (
                        archived_file.filename.find('.tif') != -1):
                    self.measurements.append({'polarisation': polarisation, 'measurement': gdal.Open(
                        os.path.join('/vsizip/' + self.datasource_path, archived_file.filename))})

                if (archived_file.filename.lower().find(polarisation.lower()) != -1) and (
                        archived_file.filename.find('calibration/calibration') != -1) and (
                        archived_file.filename.find('.xml') != -1):
                    print (polarisation)
                    print (archived_file.filename)
                    self.radiometric_correction_LUT.append(
                        {'polarisation': polarisation, 'LUT': self.archive.read(archived_file)})

                if (archived_file.filename.lower().find(polarisation.lower()) != -1) and (
                        archived_file.filename.find('noise') != -1) and (archived_file.filename.find('.xml') != -1):
                    self.noise_correction_LUT.append(
                        {'polarisation': polarisation, 'LUT': self.archive.read(archived_file)})

    def export_to_l2(self, output_directory, polarisations=None, projection='+proj=longlat +datum=WGS84 +no_defs',
                     preview_mode=False, x_scale=5, y_scale=5):
        created_files = []
        if not polarisations:
            polarisations = []
            for measurement in self.measurements:
                polarisations.append(measurement['polarisation'])

        for measurement in self.measurements:
            if measurement['polarisation'] in polarisations:
                if preview_mode:
                    current_name = os.path.join(output_directory,
                                                self.datasource_body + '_' + measurement['polarisation'] + '.jpeg')
                else:
                    current_name = os.path.join(output_directory,
                                                self.datasource_body + '_' + measurement['polarisation'] + '.tif')

                current_measurement = self.__gcp_raster_to_projected(measurement['measurement'])
                current_measurement = self.__reproject_raster_to_projection(current_measurement, projection)

                if preview_mode:
                    pixel_size = self.__get_raster_pixel_size(current_measurement)
                    current_measurement = self.__set_raster_resolution(current_measurement,
                                                                       pixel_size['xSize'] * x_scale,
                                                                       pixel_size['ySize'] * y_scale)
                    self.__save_raster_to_jpeg(current_measurement, current_name)
                else:
                    self.__save_raster_to_gtiff(current_measurement, current_name)
                created_files.append(current_name)

        return created_files

    def render(self, output_directory, polarisations=None, projection='EPSG:3857', preserve_source_file_name=False, max_size=10):
        input_fnames = self.export_to_l2(output_directory, polarisations, projection)
        for fname in input_fnames:
            try:
                ds = gdal.Open(fname)
                band = ds.GetRasterBand(1)
                (band_min, band_max, mean, stddev) = band.GetStatistics(0, 1)
                # compute Mean+-2StdDev
                # min2sigma = mean - 2 * stddev  # will be negative, don't use it
                max2sigma = mean + 2 * stddev
                # scale values to 8-bit range and set NoDataValue to 0:
                tmp1_fname = fname[:-4] + '_tmp1.tif'
                # self.__scale_values_to_byte(dataset, min2sigma, max2sigma, tmp1_fname)
                self.__scale_values_to_byte(ds, 0, max2sigma, tmp1_fname)
                # coarse the resolution:
                tmp2_fname = fname[:-4] + '_tmp2.tif'
                # we make -dstnodata none instead of 0 to compute the alpha band:
                cmd = 'gdalwarp -tr 600 1200 -srcnodata 0 -dstnodata none %s %s' % (tmp1_fname, tmp2_fname)
                print cmd
                os.system(cmd)
                # calc alpha band from nodata
                tmp_alpha_name = fname[:-4] + '_tmpAlpha.tif'
                cmd = 'gdal_calc.py -A %s --outfile=%s --calc="(A==0)*0 + (A>0)*255"' % (tmp2_fname, tmp_alpha_name)
                print cmd
                os.system(cmd)
                # make the 4-band render and DEFLATE it
                print cmd
                render_fname = "/".join(fname.split("/")[:-1]) + "/niersc_s1_" + fname.split("/")[-1].split("_")[4] + "_render.tif"
                if preserve_source_file_name:
                    render_fname = "/".join(fname.split("/")[:-1]) + "/" + fname.split("/")[-1][:-7] + ".tif"
                    print "RENDER FILE NAME: %s" % render_fname
                cmd = 'gdal_merge.py -separate -n 0 -a_nodata 0 -o %(out)s %(in)s %(in)s %(in)s %(alpha)s -co PHOTOMETRIC=RGB -co COMPRESS=DEFLATE' % {"out": render_fname, "in": tmp2_fname, "alpha": tmp_alpha_name}
                print cmd
                os.system(cmd)
                del ds
                self.optimize_size(render_fname, max_size=max_size)
                # comment this clean-up for debugging:
                os.remove(tmp1_fname)
                os.remove(tmp2_fname)
                os.remove(tmp_alpha_name)
                os.remove(fname)
                aux_xmls = glob.glob(os.path.join(output_directory, "*.aux.xml"))
                print os.path.join(output_directory, "*.aux.xml")
                for xml in aux_xmls:
                    print xml
                    os.remove(xml)
            except Exception as e:
                print type(e)

    def optimize_size(self, raster_file, max_size=10, max_iter=5):
        iter_count = 0
        while os.path.getsize(raster_file) > max_size * 1024 * 1024 and iter_count < max_iter:
            x_res, y_res = self.__get_pixel_size(raster_file)
            x_res = 1.4 * x_res
            y_res = 1.4 * y_res
            in_file = "%s_old.tif" % raster_file[:-4]
            os.rename(raster_file, in_file)
            out_file = raster_file
            cmd = 'gdalwarp -co COMPRESS=DEFLATE -tr %(x_res)s %(y_res)s %(in)s %(out)s' % {"in": in_file,
                                                                                            "out": out_file,
                                                                                            "x_res": x_res,
                                                                                            "y_res": y_res}
            os.system(cmd)
            os.remove(in_file)
            raster_file = out_file
            iter_count += 1
        else:
            x_res, y_res = self.__get_pixel_size(raster_file)
            print "After %d iterations raster resolution is %.0fx%.0f (m per pix)" % (iter_count, x_res, y_res)

    @staticmethod
    def __get_pixel_size(raster_file):
        ds = gdal.Open(raster_file)
        gt = ds.GetGeoTransform()
        return gt[1], -gt[5]

    def perform_radiometric_calibration(self, parameter='sigma', polarisations=None):
        if parameter not in ['sigma', 'beta', 'gamma']:
            print ('Invalid parameter')
            return 0
        else:
            parameter = parameter + 'Nought'

        if not polarisations:
            polarisations = []
            for measurement in self.measurements:
                polarisations.append(measurement['polarisation'])

        for measurement in self.measurements:
            cols = measurement['measurement'].RasterXSize
            rows = measurement['measurement'].RasterYSize
            if measurement['polarisation'] in polarisations:
                # print measurement['polarisation']
                lut_xml = \
                self.__dict_search(self.radiometric_correction_LUT, 'polarisation', measurement['polarisation'])[0][
                    'LUT']
                # print lut_xml
                full_lut = self.__get_full_coefficients_array(lut_xml, 'radiometric', parameter, cols, rows)

                measurement_array = np.array(
                    measurement['measurement'].GetRasterBand(1).ReadAsArray().astype(np.float32))
                # np.save('/home/silent/ma.npy',measurement_array)
                # np.save('/home/silent/lut.npy',full_lut)
                calibrated_array = (measurement_array * measurement_array) / (full_lut * full_lut)
                self.measurements.append({'polarisation': measurement['polarisation'] + '_' + parameter,
                                          'measurement': self.__create_mem_raster_based_on_existing(
                                              measurement['measurement'], [calibrated_array], gdal.GDT_Float32)})

    def __get_metadata(self):
        namespaces = {'safe': '{http://www.esa.int/safe/sentinel-1.0}'}
        self.metadata_raw = ''
        self.archive = zipfile.ZipFile(self.datasource_path)
        self.file_list = self.archive.filelist
        for archived_file in self.file_list:
            if archived_file.filename.find('manifest.safe') != -1:
                self.metadata_raw = self.archive.read(archived_file)

        if not self.metadata_raw:
            print ('error while reading metadata!')
            return 0

        self.metadata_xml_raw = etree.ElementTree(etree.fromstring(self.metadata_raw)).getroot()

        metadata_section = self.metadata_xml_raw.find('metadataSection')
        for metadata_object in metadata_section.findall('metadataObject'):
            if metadata_object.attrib['ID'] == 'platform':
                self.metadata['mode'] = metadata_object.find('metadataWrap').find('xmlData').find(
                    '{http://www.esa.int/safe/sentinel-1.0}platform').find(
                    '{http://www.esa.int/safe/sentinel-1.0}instrument').find(
                    '{http://www.esa.int/safe/sentinel-1.0}extension').find(
                    '{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}instrumentMode').find(
                    '{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}mode').text
                self.metadata['satellite_family'] = metadata_object.find('metadataWrap').find('xmlData').find(
                    '{http://www.esa.int/safe/sentinel-1.0}platform').find(
                    '{http://www.esa.int/safe/sentinel-1.0}familyName').text
                self.metadata['satellite_name'] = metadata_object.find('metadataWrap').find('xmlData').find(
                    '{http://www.esa.int/safe/sentinel-1.0}platform').find(
                    '{http://www.esa.int/safe/sentinel-1.0}number').text

            if metadata_object.attrib['ID'] == 'generalProductInformation':
                self.metadata['instrument_configuration_id'] = metadata_object.find('metadataWrap').find(
                    'xmlData').find(
                    '{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}standAloneProductInformation').find(
                    '{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}instrumentConfigurationID').text

                self.metadata['polarisations'] = []
                for polarisation_node in metadata_object.find('metadataWrap').find('xmlData').find(
                        '{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}standAloneProductInformation').findall(
                        '{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}transmitterReceiverPolarisation'):
                    self.metadata['polarisations'].append(polarisation_node.text)

                self.metadata['product_type'] = metadata_object.find('metadataWrap').find('xmlData').find(
                    '{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}standAloneProductInformation').find(
                    '{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}productType').text

            if metadata_object.attrib['ID'] == 'acquisitionPeriod':
                self.metadata['start_time'] = metadata_object.find('metadataWrap').find('xmlData').find(
                    '{http://www.esa.int/safe/sentinel-1.0}acquisitionPeriod').find(
                    '{http://www.esa.int/safe/sentinel-1.0}startTime').text
                self.metadata['start_time'] = datetime.strptime(self.metadata['start_time'], '%Y-%m-%dT%H:%M:%S.%f')
                self.metadata['stop_time'] = metadata_object.find('metadataWrap').find('xmlData').find(
                    '{http://www.esa.int/safe/sentinel-1.0}acquisitionPeriod').find(
                    '{http://www.esa.int/safe/sentinel-1.0}stopTime').text
                self.metadata['stop_time'] = datetime.strptime(self.metadata['stop_time'], '%Y-%m-%dT%H:%M:%S.%f')

            if metadata_object.attrib['ID'] == 'measurementFrameSet':
                self.metadata['foot_print'] = metadata_object.find('metadataWrap').find('xmlData').find(
                    '{http://www.esa.int/safe/sentinel-1.0}frameSet').find(
                    '{http://www.esa.int/safe/sentinel-1.0}frame').find(
                    '{http://www.esa.int/safe/sentinel-1.0}footPrint').find(
                    '{http://www.opengis.net/gml}coordinates').text
                self.metadata['foot_print'] = self.__gml_polygon_to_wkt(self.metadata['foot_print'])

    def __gml_polygon_to_wkt(self, gml_coordinates):
        pairs = gml_coordinates.split(' ')
        wkt_pairs = []
        for pair in pairs:
            coords = pair.split(',')
            wkt_pair = coords[1] + ' ' + coords[0]
            wkt_pairs.append(wkt_pair)

        wkt_pairs.append(wkt_pairs[0])
        return 'Polygon((' + ','.join(wkt_pairs) + '))'

    def __save_raster_to_gtiff(self, raster, gtiff_path):
        driver = gdal.GetDriverByName("GTiff")
        data_type = raster.GetRasterBand(1).DataType
        dataset = driver.Create(gtiff_path, raster.RasterXSize, raster.RasterYSize, raster.RasterCount, data_type)
        dataset.SetProjection(raster.GetProjection())
        dataset.SetGeoTransform(raster.GetGeoTransform())
        i = 1
        while i <= raster.RasterCount:
            dataset.GetRasterBand(i).WriteArray(raster.GetRasterBand(i).ReadAsArray())
            i += 1
        del dataset

    def __create_mem_raster_based_on_existing(self, raster, new_arrays, new_type=None):
        driver = gdal.GetDriverByName("MEM")
        if not new_type:
            data_type = raster.GetRasterBand(1).DataType
        else:
            data_type = new_type

        dataset = driver.Create('', raster.RasterXSize, raster.RasterYSize, raster.RasterCount, data_type)
        # print ('===')
        # print (raster.GetGCPs())
        # print (raster.GetGCPProjection())
        if len(raster.GetGCPs()) > 0:
            dataset.SetGCPs(raster.GetGCPs(), raster.GetGCPProjection())

        else:
            dataset.SetProjection(raster.GetProjection())
            dataset.SetGeoTransform(raster.GetGeoTransform())

        i = 1
        while i <= raster.RasterCount:
            dataset.GetRasterBand(i).WriteArray(new_arrays[i - 1])
            i += 1
        return dataset

    def __scale_values_to_byte(self, input_dataset, src_min, src_max, output_fname):
        # gdal.Translate(destName, srcDS, **kwargs)
        output_dataset = gdal.Translate(output_fname, input_dataset, format='GTiff', outputType=gdal.GDT_Byte,
                                        scaleParams=[[src_min, src_max]], noData=0)
        return output_dataset

    def __save_raster_to_jpeg(self, raster, jpeg_path):
        gdal.Translate(jpeg_path, raster, format='JPEG')

    def __reproject_raster_to_projection(self, raster, dest_projection):
        source_projection = self.__get_projection(raster)
        output_raster = gdal.Warp('', raster, srcSRS=source_projection, dstSRS=dest_projection, format='MEM')
        return output_raster

    def __get_projection(self, raster):
        return raster.GetProjection()

    def __set_raster_resolution(self, raster, xRes, yRes):
        outraster = gdal.Warp('', raster, format='MEM', xRes=xRes, yRes=yRes)
        return outraster

    def __gcp_raster_to_projected(self, raster):
        output_raster = gdal.Warp('', raster, format='MEM')
        return output_raster

    def __get_raster_pixel_size(self, raster):
        geotransform = raster.GetGeoTransform()
        return {'xSize': geotransform[1], 'ySize': geotransform[5]}

    def __get_full_coefficients_array(self, xml_text, mode, parameter, cols, rows):
        if mode == 'radiometric':
            xml_element_name = 'calibrationVectorList'
        elif mode == 'noise':
            xml_element_name = 'noiseRangeVectorList'

        coefficients_rows = []

        e = etree.ElementTree(etree.fromstring(xml_text)).getroot()
        for noiseVectorList in e.findall(xml_element_name):
            for child in noiseVectorList:
                for param in child:
                    if param.tag == 'pixel':
                        current_pixels = str(param.text).split()
                    if param.tag == parameter:
                        current_values = str(param.text).split()

                i = 0
                current_row = np.empty([1, cols])
                current_row[:] = np.nan
                while i < len(current_pixels):
                    current_row[0, int(current_pixels[i])] = float(current_values[i])
                    i += 1

                current_row = self.__fill_nan(current_row)

                coefficients_rows.append(current_row[0])

        zoom_x = float(cols) / len(coefficients_rows[0])
        zoom_y = float(rows) / len(coefficients_rows)
        return ndimage.zoom(coefficients_rows, [zoom_y, zoom_x])

    def __dict_search(self, dictionary_list, key, value):
        return [element for element in dictionary_list if element[key] == value]

    @staticmethod
    def __fill_nan(A):
        B = A
        ok = ~np.isnan(B)
        xp = ok.ravel().nonzero()[0]
        fp = B[~np.isnan(B)]
        x = np.isnan(B).ravel().nonzero()[0]
        B[np.isnan(B)] = np.interp(x, xp, fp)
        return B

    def __dump_footprint_to_geojson(self, footprint, out_dir):
        result = False
        try:
            geojson_file_path = os.path.join(out_dir, "%s.geojson" % self.datasource_body)
            g1 = shapely.wkt.loads(footprint)
            g2 = geojson.Feature(geometry=g1, properties={})
            outfile = open(geojson_file_path, "w")
            geojson.dump(g2, outfile)
            outfile.close()
            result = geojson_file_path
        except Exception as e:
            print e
            result = False
        finally:
            return result

    def create_footprint(self, out_dir):
        result = False
        try:
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            footprint = self.metadata["foot_print"]
            created_file_path = self.__dump_footprint_to_geojson(footprint, out_dir)
            result = created_file_path
        except Exception as e:
            print e
            result = False
        finally:
            return result

    def perform_noise_correction_ESA(self):
        """
        Thermal noise removal based on LUTs from metadata
        [UNDER CONSTRUCTION]
        @return: None
        """
        measurements = self.measurements[:]  # self.measurements will be modified, so don't iterate over it
        for measurement in measurements:
            cols = measurement['measurement'].RasterXSize
            rows = measurement['measurement'].RasterYSize
            print "%sx%s" % (cols, rows)

            noise = self.__get_full_coefficients_array(self.noise_correction_LUT[0]["LUT"], "noise", "noiseRangeLut", cols, rows)
            # self.__imsave(noise, title="range_noise_%s" % measurement['polarisation'], show=False)

            orig_intensity = measurement['measurement'].GetRasterBand(1).ReadAsArray() ** 2
            denoised_intensity = orig_intensity - noise
            denoised_intensity[denoised_intensity < 0] = 0
            denoised_dn = np.sqrt(denoised_intensity)
            # self.__imsave(denoised_dn, title="denoised_%s" % measurement['polarisation'], show=False)

            self.measurements.append({'polarisation': "%s_denoised" % measurement['polarisation'],
                                      'measurement': self.__create_mem_raster_based_on_existing(
                                          measurement['measurement'], [denoised_dn], gdal.GDT_Float32)})

    def __imsave(self, array, out_dir=None, title=None, show=False):
        if out_dir is None:
            out_dir = os.path.dirname(self.datasource_path)
        plt.imsave(os.path.join(out_dir, "%s.png" % title), array)
        if show:
            plt.imshow(array)
            plt.title(title)
            plt.colorbar()
            plt.show()
            plt.clf()


if __name__ == "__main__":
    s = S1L1Tools("/home/tepex/data/s1/202007/S1B_EW_GRDM_1SDH_20200715T030607_20200715T030626_022476_02AA87_BFE2.zip")
    s.perform_noise_correction_ESA()
    s.render("/home/tepex/data/s1/202007/", polarisations=["HH", "HV"], preserve_source_file_name=True)
    s.render("/home/tepex/data/s1/202007/", polarisations=["HH_denoised", "HV_denoised"], preserve_source_file_name=True)
