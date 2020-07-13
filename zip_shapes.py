import os
import os.path
import zipfile


def zip_shape(shp_path):
	"""
	Creates a zip file with a set of ESRI Shape files.
	DOES NOT validate the output.
	@param shp_path: absolute path to a file with .shp extension
	@return: an absolute path of created zip file or False if execution failed
	"""
	if os.path.isfile(shp_path):
		source_dir = os.path.dirname(shp_path)
		base = os.path.basename(shp_path).split(".")[0]  # without extension
		shx_path = os.path.join(source_dir, "%s.shx" % base)
		dbf_path = os.path.join(source_dir, "%s.dbf" % base)
		prj_path = os.path.join(source_dir, "%s.prj" % base)
		zip_path = os.path.join(source_dir, "%s.zip" % base)

		def write(zip_file, file_path):
			if os.path.isfile(file_path):
				zip_file.write(path, os.path.basename(file_path))
				os.remove(file_path)
				return True
			else:
				return False

		with zipfile.ZipFile(zip_path, "w") as z:
			for path in (shp_path, shx_path, dbf_path, prj_path):
				write(z, path)

		return zip_path
	else:
		return False


if __name__ == "__main__":
	zip_shape("/home/tepex/NIERSC/IEPI/20200712/test/niersc_ice_class_3d_20200712.shp")
