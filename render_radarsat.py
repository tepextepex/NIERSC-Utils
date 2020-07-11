import gdal
import os
import os.path


def __scale_values_to_byte(input_dataset, src_min, src_max, output_fname):
	output_dataset = gdal.Translate(output_fname, input_dataset, format='GTiff', outputType=gdal.GDT_Byte, scaleParams=[[src_min, src_max]], noData=0)
	return output_dataset


def render_radarsat(source_zip):
	input_tif = os.path.join(os.path.dirname(source_zip), "imagery_HH.tif")
	fname = input_tif
	try:
		dataset = gdal.Open("/vsizip/%s/imagery_HH.tif" % source_zip)
		band = dataset.GetRasterBand(1)
		(band_min, band_max, mean, stddev) = band.GetStatistics(0, 1)
		max2sigma = mean + 2 * stddev

		# scale values to 8-bit range and set NoDataValue to 0:
		tmp1_fname = fname[:-4] + '_tmp1.tif'
		__scale_values_to_byte(dataset, 0, max2sigma, tmp1_fname)

		# coarse the resolution:
		tmp2_fname = fname[:-4] + '_tmp2.tif'
		cmd = 'gdalwarp -overwrite -t_srs EPSG:3857 -r near -wm 1000 -multi -tr 600 600 -srcnodata 0 -dstnodata none %s %s' % (tmp1_fname, tmp2_fname)  # -dstnodata none instead of 0 is for calculating alpha band
		print cmd
		os.system(cmd)

		# calc alpha band from nodata
		tmp_alpha_fname = fname[:-4] + '_tmpAlpha.tif'
		cmd = 'gdal_calc.py -A %s --outfile=%s --calc="(A==0)*0 + (A>0)*255"' % (tmp2_fname, tmp_alpha_fname)
		print cmd
		os.system(cmd)
		
		# make the 4-band render and DEFLATE it
		render_fname = "%s.tif" % source_zip[:-4]
		cmd = 'gdal_merge.py -separate -n 0 -a_nodata 0 -o %s %s %s %s %s -co PHOTOMETRIC=RGB -co COMPRESS=DEFLATE' % (render_fname, tmp2_fname, tmp2_fname, tmp2_fname, tmp_alpha_fname)
		print cmd
		os.system(cmd)

		del dataset     

		# comment this clean-up for debugging:
		os.remove(tmp1_fname)
		os.remove(tmp2_fname)
		os.remove(tmp_alpha_fname)
		for f in os.listdir(os.path.dirname(source_zip)):
			if f[:-8] == ".aux.xml":
				os.remove(os.path.join(os.path.dirname(source_zip), f))
			
	except Exception as e:
		print type(e)


if __name__ == "__main__":
	test_zip = "/home/tepex/NIERSC/IEPI/scripts/radarsat_render/RS2_OK000000_PK0829435_DK0200709_SCNB_20200709_022702_HHHV_SGF.zip"
	render_radarsat(test_zip)
