import os
import os.path
# import requests
import json
from datetime import datetime, timedelta
import urllib
from zip_shapes import zip_shape


class AMSR2Downloader:

    # base_url = "https://seaice.uni-bremen.de/data/amsr2/bootstrap_daygrid"
    base_url = "https://seaice.uni-bremen.de/data/amsr2/asi_daygrid_swath"  # without trailing slash

    def __init__(self):
        self.current_dir = os.path.dirname(os.path.realpath(__file__))
        self.downloaded = []
        config_file = open(os.path.join(self.current_dir, 'config_amsr.json'))
        config_content = config_file.read()
        config = json.loads(config_content)

        if config["start_date"] == "":
            self.start_date = datetime.now()  # - timedelta(days=1)  # today's AMSR product may be not created yet
        else:
            self.start_date = datetime.strptime(config["start_date"], '%Y%m%d')

        if config["end_date"] == "":
            self.end_date = datetime.now()
        else:
            self.end_date = datetime.strptime(config["end_date"], '%Y%m%d')

        self.region = config["region"]
        if self.region not in ["Arctic", "Antarctic"]:
            self.log("Incorrect region! Setting default region (Arctic)")
            self.region = "Arctic"

        if config["output_dir"] == "" or not os.path.isdir(config["output_dir"]):
            self.output_dir = self.current_dir
        else:
            self.output_dir = config["output_dir"]

    def download(self):
        self.downloaded = []
        if self.region == "Antarctic":
            region_dir = "s6250"
            # place_dir = "Antarctic3125"
            place_dir = "Antarctic"
        else:
            region_dir = "n6250"
            # place_dir = "Arctic3125"
            place_dir = "Arctic"

        dates_difference = (self.end_date - self.start_date).days

        for i in range(0, dates_difference + 1):
            current_date = self.start_date + timedelta(days=i)
            month_abv = current_date.strftime("%b").lower()
            year_abv = current_date.strftime("%Y")
            date_abv = current_date.strftime("%Y%m%d")
            # file_name = 'bootstrap-AMSR2-%s-%s.hdf' % (region_dir, date_abv)
            file_name = "asi-AMSR2-%s-%s-v5.4.tif" % (region_dir, date_abv)
            current_file_url = "%s/%s/%s/%s/%s/%s" % (self.base_url, region_dir, year_abv, month_abv, place_dir, file_name)

            if os.path.exists(os.path.join(self.output_dir, file_name)):
                continue

            self.log("%s: Downloading %s..." % (datetime.now(), current_file_url))

            status = urllib.urlopen(current_file_url)
            if status.getcode() == 200:
                try:
                    local_path = os.path.join(self.output_dir, file_name)
                    urllib.urlretrieve(current_file_url, local_path)
                    self.log(" - OK!")
                    self.downloaded.append(local_path)
                except Exception as e:
                    self.log(" - exception occurred while downloading: %s" % e)
            else:
                self.log(" - Error! Status %s" % status.getcode())

        return self.downloaded

    def log(self, msg):
        with open(os.path.join(self.current_dir, "log_amsr.txt"), "a") as log_file:
            log_file.write("\n%s" % msg)

    def post_process(self):
        for asi_path in self.downloaded:
            self.extract_vector_edge(asi_path)
        self.downloaded = []

    def extract_vector_edge(self, asi_full_path, ice_conc_threshold=15):
        threshold = self.__validate_threshold(ice_conc_threshold)
        in_file = asi_full_path
        out_file = "%s_recalc.tif" % asi_full_path[:-4]
        cmd = 'gdal_calc.py -A %s --outfile=%s --NoDataValue=0 --calc="0+logical_and(A>=%s,A<=100)*1"' % (in_file, out_file, threshold)
        os.system(cmd)

        in_file = out_file
        out_file = "%s_crop.tif" % asi_full_path[:-4]
        cmd = 'gdalwarp -te 0 -610000 3750000 1330000 -tr 3125 3125 -r near -dstalpha %s %s' % (in_file, out_file)
        os.system(cmd)

        in_file = out_file
        out_file = "%s-edge15.shp" % asi_full_path[:-4]
        cmd = 'gdal_polygonize.py %s -f "ESRI Shapefile" %s ICE' % (in_file, out_file)
        os.system(cmd)

        zip_shape(out_file)

    @staticmethod
    def __validate_threshold(value):
        """
        Checks if input value is numeric and is between 0 and 100
        Otherwise returns default value of 15%
        @return: value between 0-100
        """
        try:
            float(value)
            if (value > 0) and (value < 100):
                return float(value)
            else:
                raise ValueError
        except ValueError:
            return 15  # default value (percent)


if __name__ == "__main__":
    a = AMSR2Downloader()
    # a.download()
    # a.post_process()
    # or if you already have the needed raster file:
    a.extract_vector_edge("/home/tepex/NIERSC/IEPI/AMSR/asi-AMSR2-n6250-20200721-v5.4.tif")
